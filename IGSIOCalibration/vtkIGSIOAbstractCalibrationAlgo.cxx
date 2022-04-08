/*=IGSIO=header=begin======================================================
  Program: IGSIO
  Copyright (c) Laboratory for Percutaneous Surgery. All rights reserved.
  See License.txt for details.
=========================================================IGSIO=header=end*/

// IGSIO includes
#include "igsioConfigure.h"
#include "igsioMath.h"
#include "vtkIGSIOAbstractCalibrationAlgo.h"
#include "vtkIGSIOTransformRepository.h"
#include <vtkIGSIOAccurateTimer.h>

// VTK includes
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkTransform.h"
#include "vtkXMLUtilities.h"
#include "vtksys/SystemTools.hxx"

static const double PARALLEL_ANGLE_THRESHOLD_DEGREES = 20.0;
// Note: If the needle orientation protocol changes, only the definitions of shaftAxis and secondaryAxes need to be changed
// Define the shaft axis and the secondary shaft axis
// Current needle orientation protocol dictates: shaft axis -z, orthogonal axis +x
// If StylusX is parallel to ShaftAxis then: shaft axis -z, orthogonal axis +y
static const double SHAFT_AXIS[3] = { 0, 0, -1 };
static const double ORTHOGONAL_AXIS[3] = { 1, 0, 0 };
static const double BACKUP_AXIS[3] = { 0, 1, 0 };

//-----------------------------------------------------------------------------
vtkIGSIOAbstractCalibrationAlgo::vtkIGSIOAbstractCalibrationAlgo()
{
  this->PivotPointToMarkerTransformMatrix = NULL;

  this->PivotPointPosition_Reference[0] = 0.0;
  this->PivotPointPosition_Reference[1] = 0.0;
  this->PivotPointPosition_Reference[2] = 0.0;
  this->PivotPointPosition_Reference[3] = 1.0;

  this->PreviousMarkerToReferenceTransformMatrix = vtkMatrix4x4::New();

  this->ObjectMarkerCoordinateFrame = NULL;
  this->ReferenceCoordinateFrame = NULL;

  this->ErrorCode = CALIBRATION_NOT_STARTED;

  this->MinimumOrientationDifferenceDeg = 15.0;
  this->PositionDifferenceThresholdMm = 0.0;
  this->OrientationDifferenceThresholdDegrees = 0.0;

  this->PoseBucketSize = -1;
  this->MaximumNumberOfPoseBuckets = 1;
  this->MaximumPoseBucketError = 3.0;
}

//-----------------------------------------------------------------------------
vtkIGSIOAbstractCalibrationAlgo::~vtkIGSIOAbstractCalibrationAlgo()
{
  this->PreviousMarkerToReferenceTransformMatrix->Delete();
  this->RemoveAllCalibrationPoints();
}

//-----------------------------------------------------------------------------
void vtkIGSIOAbstractCalibrationAlgo::RemoveAllCalibrationPoints()
{
  this->MarkerToReferenceTransformMatrixBuckets.clear();
  this->SetErrorCode(CALIBRATION_NOT_STARTED);
  this->OutlierIndices.clear();
}

//----------------------------------------------------------------------------
igsioStatus vtkIGSIOAbstractCalibrationAlgo::InsertNextCalibrationPoint(vtkMatrix4x4* aMarkerToReferenceTransformMatrix)
{
  double positionDifferenceMm = this->PositionDifferenceThresholdMm;
  double orientationDifferenceDegrees = this->OrientationDifferenceThresholdDegrees;
  if (this->MarkerToReferenceTransformMatrixBuckets.size() >= 1 && this->MarkerToReferenceTransformMatrixBuckets[0].MarkerToReferenceCalibrationPoints.size() > 0)
  {
    // Compute position and orientation difference of current and previous positions
    positionDifferenceMm = igsioMath::GetPositionDifference(aMarkerToReferenceTransformMatrix, this->PreviousMarkerToReferenceTransformMatrix);
    orientationDifferenceDegrees = igsioMath::GetOrientationDifference(aMarkerToReferenceTransformMatrix, this->PreviousMarkerToReferenceTransformMatrix);
  }
  if (positionDifferenceMm < this->PositionDifferenceThresholdMm || orientationDifferenceDegrees < this->OrientationDifferenceThresholdDegrees)
  {
    LOG_DEBUG("Acquired position is too close to the previous - it is skipped");
    return IGSIO_FAIL;
  }
  this->PreviousMarkerToReferenceTransformMatrix->DeepCopy(aMarkerToReferenceTransformMatrix);

  MarkerToReferenceTransformMatrixBucket* currentBucket = NULL;
  if (this->MarkerToReferenceTransformMatrixBuckets.size() > 0)
  {
    // Last bucket is the current one
    currentBucket = &this->MarkerToReferenceTransformMatrixBuckets[this->MarkerToReferenceTransformMatrixBuckets.size() - 1];
  }

  if (this->PoseBucketSize > 0 && currentBucket && currentBucket->MarkerToReferenceCalibrationPoints.size() >= this->PoseBucketSize)
  {
    // If the current bucket is full, then we will need to create a new one.
    currentBucket = NULL;
  }

  if (!currentBucket)
  {
    // Create a new bucket
    this->MarkerToReferenceTransformMatrixBuckets.push_back(MarkerToReferenceTransformMatrixBucket());

    if (this->MaximumNumberOfPoseBuckets > 0 && this->MarkerToReferenceTransformMatrixBuckets.size() > this->MaximumNumberOfPoseBuckets)
    {
      // Maximum number of buckets reached. Remove the oldest one.
      this->MarkerToReferenceTransformMatrixBuckets.pop_front();
    }

    // Latest bucket is the current one
    currentBucket = &this->MarkerToReferenceTransformMatrixBuckets[this->MarkerToReferenceTransformMatrixBuckets.size() - 1];
  }

  vtkNew<vtkMatrix4x4> markerToReferenceTransformMatrixCopy;
  markerToReferenceTransformMatrixCopy->DeepCopy(aMarkerToReferenceTransformMatrix);
  currentBucket->MarkerToReferenceCalibrationPoints.push_back(markerToReferenceTransformMatrixCopy);
  this->InvokeEvent(InputTransformAddedEvent);

  if (this->PoseBucketSize >= 0 && currentBucket->MarkerToReferenceCalibrationPoints.size() == this->PoseBucketSize)
  {
    // Bucket is filled. Clean the input buffer
    this->CleanInputBuffer();
  }

  return IGSIO_SUCCESS;
}

//----------------------------------------------------------------------------
igsioStatus vtkIGSIOAbstractCalibrationAlgo::CleanInputBuffer()
{
  std::vector<vtkMatrix4x4*> latestBucket;
  this->GetLatestBucketMarkerToReferenceMatrices(&latestBucket);

  vtkNew<vtkMatrix4x4> pivotPointToMarkerTransformMatrix;
  double pivotPoint_Marker[4] = { 0, 0, 0, 1 };
  double pivotPoint_Reference[4] = { 0, 0, 0, 1 };
  std::set<unsigned int> outlierIndices;

  double poseBucketError = VTK_DOUBLE_MAX;
  igsioStatus status = this->DoCalibrationInternal(&latestBucket, poseBucketError);

  if (status != IGSIO_SUCCESS || poseBucketError > this->MaximumPoseBucketError)
  {
    // The latest pivot calibration bucket does not contain pivoting.
    // Discard all poses.
    this->RemoveAllCalibrationPoints();
    this->SetErrorCode(CALIBRATION_HIGH_ERROR);
  }

  return IGSIO_SUCCESS;
}

//----------------------------------------------------------------------------
std::vector<vtkMatrix4x4*> vtkIGSIOAbstractCalibrationAlgo::GetAllMarkerToReferenceMatrices()
{
  std::vector<vtkMatrix4x4*> markerToTransformMatrixArray;
  for (int i = 0; i < this->MarkerToReferenceTransformMatrixBuckets.size(); ++i)
  {
    this->GetBucketMarkerToReferenceMatrices(i, &markerToTransformMatrixArray);
  }
  return markerToTransformMatrixArray;
}

//----------------------------------------------------------------------------
void vtkIGSIOAbstractCalibrationAlgo::GetAllMarkerToReferenceMatrices(std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray)
{
  for (int i = 0; i < this->MarkerToReferenceTransformMatrixBuckets.size(); ++i)
  {
    this->GetBucketMarkerToReferenceMatrices(i, markerToTransformMatrixArray);
  }
}

//----------------------------------------------------------------------------
std::vector<vtkMatrix4x4*> vtkIGSIOAbstractCalibrationAlgo::GetBucketMarkerToReferenceMatrices(int bucketIndex)
{
  std::vector<vtkMatrix4x4*> markerToTransformMatrixArray;
  if (this->MarkerToReferenceTransformMatrixBuckets.size() > bucketIndex)
  {
    this->GetBucketMarkerToReferenceMatrices(bucketIndex, &markerToTransformMatrixArray);
  }
  return markerToTransformMatrixArray;
}

//----------------------------------------------------------------------------
void vtkIGSIOAbstractCalibrationAlgo::GetBucketMarkerToReferenceMatrices(int bucketIndex, std::vector<vtkMatrix4x4*>* matrixArray)
{
  if (this->MarkerToReferenceTransformMatrixBuckets.size() <= bucketIndex)
  {
    return;
  }

  if (!matrixArray)
  {
    return;
  }

  MarkerToReferenceTransformMatrixBucket* bucket = &this->MarkerToReferenceTransformMatrixBuckets[bucketIndex];
  std::vector<vtkSmartPointer<vtkMatrix4x4> >* poses = &bucket->MarkerToReferenceCalibrationPoints;
  for (std::vector<vtkSmartPointer<vtkMatrix4x4> >::iterator poseIt = poses->begin(); poseIt != poses->end(); ++poseIt)
  {
    matrixArray->push_back(*poseIt);
  }
}

//----------------------------------------------------------------------------
std::vector<vtkMatrix4x4*> vtkIGSIOAbstractCalibrationAlgo::GetLatestBucketMarkerToReferenceMatrices()
{
  this->GetBucketMarkerToReferenceMatrices(this->MarkerToReferenceTransformMatrixBuckets.size() - 1);
}

//----------------------------------------------------------------------------
void vtkIGSIOAbstractCalibrationAlgo::GetLatestBucketMarkerToReferenceMatrices(std::vector<vtkMatrix4x4*>* latestBucket)
{
  this->GetBucketMarkerToReferenceMatrices(this->MarkerToReferenceTransformMatrixBuckets.size() - 1, latestBucket);
}

//----------------------------------------------------------------------------
int vtkIGSIOAbstractCalibrationAlgo::GetNumberOfCalibrationPoints()
{
  int numberOfCalibrationPoints = 0;
  for (std::deque<vtkIGSIOAbstractCalibrationAlgo::MarkerToReferenceTransformMatrixBucket>::iterator bucketIt = this->MarkerToReferenceTransformMatrixBuckets.begin();
    bucketIt != this->MarkerToReferenceTransformMatrixBuckets.end(); ++bucketIt)
  {
    numberOfCalibrationPoints += bucketIt->MarkerToReferenceCalibrationPoints.size();
  }
  return numberOfCalibrationPoints;
}

//----------------------------------------------------------------------------
double vtkIGSIOAbstractCalibrationAlgo::GetOrientationDifferenceDeg(vtkMatrix4x4* aMatrix, vtkMatrix4x4* bMatrix)
{
  vtkSmartPointer<vtkMatrix4x4> diffMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
  vtkSmartPointer<vtkMatrix4x4> invBmatrix = vtkSmartPointer<vtkMatrix4x4>::New();

  vtkMatrix4x4::Invert(bMatrix, invBmatrix);

  vtkMatrix4x4::Multiply4x4(aMatrix, invBmatrix, diffMatrix);

  vtkSmartPointer<vtkTransform> diffTransform = vtkSmartPointer<vtkTransform>::New();
  diffTransform->SetMatrix(diffMatrix);

  double angleDiff_rad = vtkMath::RadiansFromDegrees(diffTransform->GetOrientationWXYZ()[0]);

  double normalizedAngleDiff_rad = atan2(sin(angleDiff_rad), cos(angleDiff_rad)); // normalize angle to domain -pi, pi

  return vtkMath::DegreesFromRadians(normalizedAngleDiff_rad);
}

//---------------------------------------------------------------------------
double vtkIGSIOAbstractCalibrationAlgo::GetMaximumToolOrientationDifferenceDeg()
{
  // this will store the maximum difference in orientation between the first transform and all the other transforms
  double maximumOrientationDifferenceDeg = 0;

  std::vector<vtkMatrix4x4*> toolToReferenceMatrices;
  this->GetAllMarkerToReferenceMatrices(&toolToReferenceMatrices);

  std::vector<vtkMatrix4x4*>::const_iterator matricesEnd = toolToReferenceMatrices.end();
  vtkMatrix4x4* referenceOrientationMatrix = toolToReferenceMatrices.front();
  std::vector<vtkMatrix4x4*>::const_iterator it;
  for (it = toolToReferenceMatrices.begin(); it != matricesEnd; it++)
  {
    double orientationDifferenceDeg = fabs(this->GetOrientationDifferenceDeg(referenceOrientationMatrix, (*it)));
    if (maximumOrientationDifferenceDeg < orientationDifferenceDeg)
    {
      maximumOrientationDifferenceDeg = orientationDifferenceDeg;
    }
  }

  return maximumOrientationDifferenceDeg;
}

//-----------------------------------------------------------------------------
int vtkIGSIOAbstractCalibrationAlgo::GetNumberOfDetectedOutliers()
{
  return this->OutlierIndices.size();
}

//---------------------------------------------------------------------------
void vtkIGSIOAbstractCalibrationAlgo::UpdateShaftDirection(vtkMatrix4x4* toolTipToToolMatrix)
{
  // We need to verify that the ToolTipToTool vector in the Shaft coordinate system is in the opposite direction of the shaft
  vtkSmartPointer< vtkMatrix4x4 > rotationMatrix = vtkSmartPointer< vtkMatrix4x4 >::New();
  this->GetToolTipToToolRotation(toolTipToToolMatrix, rotationMatrix);
  rotationMatrix->Invert();

  double toolTipToToolTranslation_ToolTip[4] = { 0, 0, 0, 0 }; // This is a vector, not a point, so the last element is 0
  toolTipToToolTranslation_ToolTip[0] = toolTipToToolMatrix->GetElement(0, 3);
  toolTipToToolTranslation_ToolTip[1] = toolTipToToolMatrix->GetElement(1, 3);
  toolTipToToolTranslation_ToolTip[2] = toolTipToToolMatrix->GetElement(2, 3);

  double toolTipToToolTranslation_Shaft[4] = { 0, 0, 0, 0 }; // This is a vector, not a point, so the last element is 0
  rotationMatrix->MultiplyPoint(toolTipToToolTranslation_ToolTip, toolTipToToolTranslation_Shaft);
  double toolTipToToolTranslation3_Shaft[3] = { toolTipToToolTranslation_Shaft[0], toolTipToToolTranslation_Shaft[1], toolTipToToolTranslation_Shaft[2] };

  // Check if it is parallel or opposite to shaft direction
  if (vtkMath::Dot(SHAFT_AXIS, toolTipToToolTranslation3_Shaft) > 0)
  {
    this->FlipShaftDirection(toolTipToToolMatrix);
  }
}

//---------------------------------------------------------------------------
void vtkIGSIOAbstractCalibrationAlgo::FlipShaftDirection(vtkMatrix4x4* toolTipToToolMatrix)
{
  // Need to rotate around the orthogonal axis
  vtkSmartPointer< vtkMatrix4x4 > rotationMatrix = vtkSmartPointer< vtkMatrix4x4 >::New();
  this->GetToolTipToToolRotation(toolTipToToolMatrix, rotationMatrix);

  double shaftAxis_Shaft[4] = { SHAFT_AXIS[0], SHAFT_AXIS[1], SHAFT_AXIS[2], 0 }; // This is a vector, not a point, so the last element is 0
  double shaftAxis_ToolTip[4] = { 0, 0, 0, 0 };
  rotationMatrix->MultiplyPoint(shaftAxis_Shaft, shaftAxis_ToolTip);

  vnl_vector< double > orthogonalAxis_Shaft = this->ComputeSecondaryAxis(vnl_vector< double >(3, 3, shaftAxis_ToolTip));

  vtkSmartPointer< vtkTransform > flipTransform = vtkSmartPointer< vtkTransform >::New();
  flipTransform->RotateWXYZ(180, orthogonalAxis_Shaft.get(0), orthogonalAxis_Shaft.get(1), orthogonalAxis_Shaft.get(2));
  vtkSmartPointer< vtkTransform > originalTransform = vtkSmartPointer< vtkTransform >::New();
  originalTransform->SetMatrix(toolTipToToolMatrix);
  originalTransform->PreMultiply();
  originalTransform->Concatenate(flipTransform);
  originalTransform->GetMatrix(toolTipToToolMatrix);
}

//---------------------------------------------------------------------------
vnl_vector< double > vtkIGSIOAbstractCalibrationAlgo::ComputeSecondaryAxis(vnl_vector< double > shaftAxis_ToolTip)
{
  // If the secondary axis 1 is parallel to the shaft axis in the tooltip frame, then use secondary axis 2
  vnl_vector< double > orthogonalAxis_Shaft(3, 3, ORTHOGONAL_AXIS);
  double angle = acos(dot_product(shaftAxis_ToolTip, orthogonalAxis_Shaft));
  // Force angle to be between -pi/2 and +pi/2
  if (angle > vtkMath::Pi() / 2)
  {
    angle -= vtkMath::Pi();
  }
  if (angle < -vtkMath::Pi() / 2)
  {
    angle += vtkMath::Pi();
  }

  if (fabs(angle) < vtkMath::RadiansFromDegrees(PARALLEL_ANGLE_THRESHOLD_DEGREES)) // If shaft axis and orthogonal axis are not parallel
  {
    return vnl_vector< double >(3, 3, BACKUP_AXIS);
  }
  return orthogonalAxis_Shaft;
}

//---------------------------------------------------------------------------
void vtkIGSIOAbstractCalibrationAlgo::GetToolTipToToolRotation(vtkMatrix4x4* toolTipToToolMatrix, vtkMatrix4x4* rotationMatrix)
{
  rotationMatrix->Identity();
  rotationMatrix->SetElement(0, 0, toolTipToToolMatrix->GetElement(0, 0));
  rotationMatrix->SetElement(0, 1, toolTipToToolMatrix->GetElement(0, 1));
  rotationMatrix->SetElement(0, 2, toolTipToToolMatrix->GetElement(0, 2));
  rotationMatrix->SetElement(1, 0, toolTipToToolMatrix->GetElement(1, 0));
  rotationMatrix->SetElement(1, 1, toolTipToToolMatrix->GetElement(1, 1));
  rotationMatrix->SetElement(1, 2, toolTipToToolMatrix->GetElement(1, 2));
  rotationMatrix->SetElement(2, 0, toolTipToToolMatrix->GetElement(2, 0));
  rotationMatrix->SetElement(2, 1, toolTipToToolMatrix->GetElement(2, 1));
  rotationMatrix->SetElement(2, 2, toolTipToToolMatrix->GetElement(2, 2));
}
