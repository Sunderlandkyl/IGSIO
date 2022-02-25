/*=IGSIO=header=begin======================================================
  Program: IGSIO
  Copyright (c) Laboratory for Percutaneous Surgery. All rights reserved.
  See License.txt for details.
=========================================================IGSIO=header=end*/

// IGSIO includes
#include "igsioConfigure.h"
#include "igsioMath.h"
#include "vtkIGSIOPivotCalibrationAlgo.h"
#include "vtkIGSIOTransformRepository.h"
#include <vtkIGSIOAccurateTimer.h>

// VTK includes
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkTransform.h"
#include "vtkXMLUtilities.h"
#include "vtksys/SystemTools.hxx"

// ITK includes
#include <vnl/algo/vnl_determinant.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_vector.h>

static const double PARALLEL_ANGLE_THRESHOLD_DEGREES = 20.0;
// Note: If the needle orientation protocol changes, only the definitions of shaftAxis and secondaryAxes need to be changed
// Define the shaft axis and the secondary shaft axis
// Current needle orientation protocol dictates: shaft axis -z, orthogonal axis +x
// If StylusX is parallel to ShaftAxis then: shaft axis -z, orthogonal axis +y
static const double SHAFT_AXIS[3] = { 0, 0, -1 };
static const double ORTHOGONAL_AXIS[3] = { 1, 0, 0 };
static const double BACKUP_AXIS[3] = { 0, 1, 0 };

vtkStandardNewMacro(vtkIGSIOPivotCalibrationAlgo);

//-----------------------------------------------------------------------------
vtkIGSIOPivotCalibrationAlgo::vtkIGSIOPivotCalibrationAlgo()
{
  this->PivotPointToMarkerTransformMatrix = NULL;
  this->PreviousMarkerToReferenceTransformMatrix = vtkMatrix4x4::New();

  this->PivotCalibrationErrorMm = -1.0;
  this->SpinCalibrationErrorMm = -1.0;

  this->ObjectMarkerCoordinateFrame = NULL;
  this->ReferenceCoordinateFrame = NULL;
  this->ObjectPivotPointCoordinateFrame = NULL;

  this->PivotPointPosition_Reference[0] = 0.0;
  this->PivotPointPosition_Reference[1] = 0.0;
  this->PivotPointPosition_Reference[2] = 0.0;
  this->PivotPointPosition_Reference[3] = 1.0;

  this->ErrorCode = CALIBRATION_NOT_STARTED;

  this->MinimumOrientationDifferenceDeg = 15.0;
  this->PositionDifferenceThresholdMm = 0.0;
  this->OrientationDifferenceThresholdDegrees = 0.0;
  this->AutoCalibrationEnabled = false;
  this->AutoCalibrationBucketSize = 10;
  this->AutoCalibrationNumberOfPoints = 50;
  this->AutoCalibrationTargetError = 2.0;
  this->AutoCalibrationMaximumBucketError = 3.0;
  this->AutoCalibrationMaximumNumberOfBuckets = 5;
  this->AutoCalibrationMode = PIVOT_CALIBRATION;
  this->AutoCalibrationAutoOrient = true;
  this->AutoCalibrationSnapRotation = false;
}

//-----------------------------------------------------------------------------
vtkIGSIOPivotCalibrationAlgo::~vtkIGSIOPivotCalibrationAlgo()
{
  this->SetPivotPointToMarkerTransformMatrix(NULL);
  this->PreviousMarkerToReferenceTransformMatrix->Delete();
  this->RemoveAllCalibrationPoints();
}

//-----------------------------------------------------------------------------
void vtkIGSIOPivotCalibrationAlgo::RemoveAllCalibrationPoints()
{
  this->MarkerToReferenceTransformMatrixBuckets.clear();
  this->SetErrorCode(CALIBRATION_NOT_STARTED);
  this->OutlierIndices.clear();
}

//----------------------------------------------------------------------------
igsioStatus vtkIGSIOPivotCalibrationAlgo::InsertNextCalibrationPoint(vtkMatrix4x4* aMarkerToReferenceTransformMatrix)
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

  if (this->AutoCalibrationEnabled && currentBucket && this->AutoCalibrationBucketSize > 0 && currentBucket->MarkerToReferenceCalibrationPoints.size() >= this->AutoCalibrationBucketSize)
  {
    // If the current bucket is full, then we will need to create a new one.
    currentBucket = NULL;
  }

  if (!currentBucket)
  {
    // Create a new bucket
    this->MarkerToReferenceTransformMatrixBuckets.push_back(MarkerToReferenceTransformMatrixBucket());

    if (this->AutoCalibrationEnabled && this->AutoCalibrationMaximumNumberOfBuckets > 0 && this->MarkerToReferenceTransformMatrixBuckets.size() > this->AutoCalibrationMaximumNumberOfBuckets)
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

  if (this->AutoCalibrationEnabled && this->AutoCalibrationBucketSize > 0 && currentBucket->MarkerToReferenceCalibrationPoints.size() >= this->AutoCalibrationBucketSize)
  {
    // We have filled up the current bucket. Attempt to perform calibration.
    this->AutoCalibrate();
  }
  return IGSIO_SUCCESS;
}

//----------------------------------------------------------------------------
igsioStatus vtkIGSIOPivotCalibrationAlgo::AutoCalibrate()
{
  // Check to make sure the data is good and remove all buckets if it is not.
  this->CleanInputBuffer();

  if (this->GetNumberOfCalibrationPoints() < this->AutoCalibrationNumberOfPoints)
  {
    return IGSIO_FAIL;
  }

  igsioStatus status = IGSIO_FAIL;
  double error = VTK_DOUBLE_MAX;

  if (this->GetAutoCalibrationMode() == PIVOT_CALIBRATION)
  {

    status = this->DoPivotCalibration(NULL, this->AutoCalibrationAutoOrient);
    error = this->PivotCalibrationErrorMm;
  }
  else if (this->GetAutoCalibrationMode() == SPIN_CALIBRATION)
  {
    status = this->DoSpinCalibration(NULL, this->AutoCalibrationSnapRotation, this->AutoCalibrationAutoOrient);
    error = this->SpinCalibrationErrorMm;
  }

  if (error > this->AutoCalibrationTargetError)
  {
    this->SetErrorCode(CALIBRATION_HIGH_ERROR);
  }

  if (status != IGSIO_SUCCESS || this->ErrorCode != CALIBRATION_SUCCESS)
  {
    return IGSIO_FAIL;
  }

  this->InvokeEvent(AutoCalibrationCompleteEvent);
  return IGSIO_SUCCESS;
}


//----------------------------------------------------------------------------
igsioStatus vtkIGSIOPivotCalibrationAlgo::CleanInputBuffer()
{
  std::vector<vtkMatrix4x4*> latestBucket;
  this->GetMarkerToReferenceTransformMatrixArray(this->MarkerToReferenceTransformMatrixBuckets.size() - 1, &latestBucket);

  vtkNew<vtkMatrix4x4> pivotPointToMarkerTransformMatrix;
  double pivotPoint_Marker[4] = { 0, 0, 0, 1 };
  double pivotPoint_Reference[4] = { 0, 0, 0, 1 };
  std::set<unsigned int> outlierIndices;

  igsioStatus status = IGSIO_FAIL;
  double meanError = VTK_DOUBLE_MAX;

  if (this->AutoCalibrationMode == PIVOT_CALIBRATION)
  {
    status = this->DoPivotCalibrationInternal(&latestBucket, true, &outlierIndices, pivotPoint_Marker, pivotPoint_Reference, pivotPointToMarkerTransformMatrix);
    meanError = this->ComputePivotCalibrationError(&latestBucket, &outlierIndices, pivotPoint_Reference, pivotPointToMarkerTransformMatrix);
  }
  else if (this->AutoCalibrationMode == SPIN_CALIBRATION)
  {
    status = this->DoSpinCalibrationInternal(&latestBucket, true, true, pivotPointToMarkerTransformMatrix, meanError);
  }

  if (status != IGSIO_SUCCESS || meanError > this->AutoCalibrationMaximumBucketError)
  {
    // The latest pivot calibration bucket does not contain pivoting.
    // Discard it and all earlier buckets.
    this->MarkerToReferenceTransformMatrixBuckets.clear();
    this->SetErrorCode(CALIBRATION_HIGH_ERROR);
  }

  return IGSIO_SUCCESS;
}

//----------------------------------------------------------------------------
/*
In homogeneous coordinates:
 PivotPoint_Reference = MarkerToReferenceTransformMatrix * PivotPoint_Marker

MarkerToReferenceTransformMatrix decomosed to rotation matrix and translation vector:
 PivotPoint_Reference = MarkerToReferenceTransformRotationMatrix * PivotPoint_Marker + MarkerToReferenceTransformTranslationVector
rearranged:
 MarkerToReferenceTransformRotationMatrix * PivotPoint_Marker - PivotPoint_Reference = -MarkerToReferenceTransformTranslationVector
in a matrix form:
 [ MarkerToReferenceTransformRotationMatrix | -Identity3x3 ] * [ PivotPoint_Marker    ] = [ -MarkerToReferenceTransformTranslationVector ]
                                                               [ PivotPoint_Reference ]

It's an Ax=b linear problem that can be solved with robust LSQR:
 Ai = [ MarkerToReferenceTransformRotationMatrix | -Identity3x3 ]
 xi = [ PivotPoint_Marker    ]
      [ PivotPoint_Reference ]
 bi = [ -MarkerToReferenceTransformTranslationVector ]
*/
igsioStatus vtkIGSIOPivotCalibrationAlgo::GetPivotPointPosition(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, std::set<unsigned int>* outlierIndices, double* pivotPoint_Marker, double* pivotPoint_Reference)
{
  std::vector<vnl_vector<double> > aMatrix;
  std::vector<double> bVector;
  vnl_vector<double> xVector(6, 0);   // result vector

  vnl_vector<double> aMatrixRow(6);
  for (std::vector<vtkMatrix4x4*>::const_iterator markerToReferenceTransformIt = markerToTransformMatrixArray->begin();
    markerToReferenceTransformIt != markerToTransformMatrixArray->end(); ++markerToReferenceTransformIt)
  {
    for (int i = 0; i < 3; i++)
    {
      aMatrixRow(0) = (*markerToReferenceTransformIt)->Element[i][0];
      aMatrixRow(1) = (*markerToReferenceTransformIt)->Element[i][1];
      aMatrixRow(2) = (*markerToReferenceTransformIt)->Element[i][2];
      aMatrixRow(3) = (i == 0 ? -1 : 0);
      aMatrixRow(4) = (i == 1 ? -1 : 0);
      aMatrixRow(5) = (i == 2 ? -1 : 0);
      aMatrix.push_back(aMatrixRow);
    }
    bVector.push_back(-(*markerToReferenceTransformIt)->Element[0][3]);
    bVector.push_back(-(*markerToReferenceTransformIt)->Element[1][3]);
    bVector.push_back(-(*markerToReferenceTransformIt)->Element[2][3]);
  }

  double mean = 0;
  double stdev = 0;
  vnl_vector<unsigned int> notOutliersIndices;
  notOutliersIndices.clear();
  notOutliersIndices.set_size(bVector.size());
  for (unsigned int i = 0; i < bVector.size(); ++i)
  {
    notOutliersIndices.put(i, i);
  }
  if (igsioMath::LSQRMinimize(aMatrix, bVector, xVector, &mean, &stdev, &notOutliersIndices) != IGSIO_SUCCESS)
  {
    LOG_ERROR("vtkIGSIOPivotCalibrationAlgo failed: LSQRMinimize error");
    return IGSIO_FAIL;
  }

  // Note: Outliers are detected and rejected for each row (each coordinate axis). Although most frequently
  // an outlier sample's every component is an outlier, it may be possible that only certain components of an
  // outlier sample are removed, which may be desirable for some cases (e.g., when the point is an outlier because
  // temporarily it was flipped around one axis) and not desirable for others (when the point is completely randomly
  // corrupted), but there would be no measurable difference anyway if the only a few percent of the points are
  // outliers.

  outlierIndices->clear();
  unsigned int processFromRowIndex = 0;
  for (unsigned int i = 0; i < notOutliersIndices.size(); i++)
  {
    unsigned int nextNotOutlierRowIndex = notOutliersIndices[i];
    if (nextNotOutlierRowIndex > processFromRowIndex)
    {
      // samples were missed, so they are outliers
      for (unsigned int outlierRowIndex = processFromRowIndex; outlierRowIndex < nextNotOutlierRowIndex; outlierRowIndex++)
      {
        int sampleIndex = outlierRowIndex / 3; // 3 rows are generated per sample
        outlierIndices->insert(sampleIndex);
      }
    }
    processFromRowIndex = nextNotOutlierRowIndex + 1;
  }

  pivotPoint_Marker[0] = xVector[0];
  pivotPoint_Marker[1] = xVector[1];
  pivotPoint_Marker[2] = xVector[2];

  pivotPoint_Reference[0] = xVector[3];
  pivotPoint_Reference[1] = xVector[4];
  pivotPoint_Reference[2] = xVector[5];

  return IGSIO_SUCCESS;
}

//----------------------------------------------------------------------------
std::vector<vtkMatrix4x4*> vtkIGSIOPivotCalibrationAlgo::GetMarkerToReferenceTransformMatrixArray()
{
  std::vector<vtkMatrix4x4*> markerToTransformMatrixArray;
  for (int i = 0; i < this->MarkerToReferenceTransformMatrixBuckets.size(); ++i)
  {
    this->GetMarkerToReferenceTransformMatrixArray(i, &markerToTransformMatrixArray);
  }
  return markerToTransformMatrixArray;
}

//----------------------------------------------------------------------------
void vtkIGSIOPivotCalibrationAlgo::GetMarkerToReferenceTransformMatrixArray(std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray)
{
  for (int i = 0; i < this->MarkerToReferenceTransformMatrixBuckets.size(); ++i)
  {
    this->GetMarkerToReferenceTransformMatrixArray(i, markerToTransformMatrixArray);
  }
}

//----------------------------------------------------------------------------
std::vector<vtkMatrix4x4*> vtkIGSIOPivotCalibrationAlgo::GetMarkerToReferenceTransformMatrixArray(int bucketIndex)
{
  std::vector<vtkMatrix4x4*> markerToTransformMatrixArray;
  if (this->MarkerToReferenceTransformMatrixBuckets.size() > bucketIndex)
  {
    this->GetMarkerToReferenceTransformMatrixArray(bucketIndex, &markerToTransformMatrixArray);
  }
  return markerToTransformMatrixArray;
}

//----------------------------------------------------------------------------
void vtkIGSIOPivotCalibrationAlgo::GetMarkerToReferenceTransformMatrixArray(int bucketIndex, std::vector<vtkMatrix4x4*>* matrixArray)
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
int vtkIGSIOPivotCalibrationAlgo::GetNumberOfCalibrationPoints()
{
  int numberOfCalibrationPoints = 0;
  for (std::deque<vtkIGSIOPivotCalibrationAlgo::MarkerToReferenceTransformMatrixBucket>::iterator bucketIt = this->MarkerToReferenceTransformMatrixBuckets.begin();
    bucketIt != this->MarkerToReferenceTransformMatrixBuckets.end(); ++bucketIt)
  {
    numberOfCalibrationPoints += bucketIt->MarkerToReferenceCalibrationPoints.size();
  }
  return numberOfCalibrationPoints;
}

//----------------------------------------------------------------------------
double vtkIGSIOPivotCalibrationAlgo::GetOrientationDifferenceDeg(vtkMatrix4x4* aMatrix, vtkMatrix4x4* bMatrix)
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
double vtkIGSIOPivotCalibrationAlgo::GetMaximumToolOrientationDifferenceDeg()
{
  // this will store the maximum difference in orientation between the first transform and all the other transforms
  double maximumOrientationDifferenceDeg = 0;

  std::vector<vtkMatrix4x4*> toolToReferenceMatrices;
  this->GetMarkerToReferenceTransformMatrixArray(&toolToReferenceMatrices);
  std::vector<vtkMatrix4x4*>::const_iterator matricesEnd = toolToReferenceMatrices.end();
  vtkMatrix4x4* referenceOrientationMatrix = toolToReferenceMatrices.front();
  std::vector<vtkMatrix4x4*>::const_iterator it;
  for (it = toolToReferenceMatrices.begin(); it != matricesEnd; it++)
  {
    double orientationDifferenceDeg = this->GetOrientationDifferenceDeg(referenceOrientationMatrix, (*it));
    if (maximumOrientationDifferenceDeg < orientationDifferenceDeg)
    {
      maximumOrientationDifferenceDeg = orientationDifferenceDeg;
    }
  }

  return maximumOrientationDifferenceDeg;
}

//----------------------------------------------------------------------------
igsioStatus vtkIGSIOPivotCalibrationAlgo::DoPivotCalibration(vtkIGSIOTransformRepository* aTransformRepository/* = NULL*/, bool autoOrient)
{
  if (this->MarkerToReferenceTransformMatrixBuckets.empty())
  {
    this->SetErrorCode(CALIBRATION_NOT_ENOUGH_POINTS);
    return IGSIO_FAIL;
  }

  if (this->GetMaximumToolOrientationDifferenceDeg() < this->MinimumOrientationDifferenceDeg)
  {
    this->SetErrorCode(CALIBRATION_NOT_ENOUGH_VARIATION);
    return IGSIO_FAIL;
  }

  std::vector<vtkMatrix4x4*> markerToTransformMatrixArray = this->GetMarkerToReferenceTransformMatrixArray();

  double pivotPoint_Marker[4] = { 0.0, 0.0, 0.0, 1.0 };
  double pivotPoint_Reference[4] = { 0.0, 0.0, 0.0, 1.0 };
  vtkNew<vtkMatrix4x4> pivotPointToMarkerTransformMatrix;
  igsioStatus status = this->DoPivotCalibrationInternal(&markerToTransformMatrixArray, autoOrient, &this->OutlierIndices, pivotPoint_Marker, pivotPoint_Reference, pivotPointToMarkerTransformMatrix);
  if (status == IGSIO_SUCCESS)
  {
    this->SetErrorCode(CALIBRATION_SUCCESS);
  }
  else
  {
    this->SetErrorCode(CALIBRATION_FAIL);
  }

  this->SetPivotPointToMarkerTransformMatrix(pivotPointToMarkerTransformMatrix);
  this->ComputePivotCalibrationError();

  // Save result
  if (aTransformRepository)
  {
    igsioTransformName pivotPointToMarkerTransformName(this->ObjectPivotPointCoordinateFrame, this->ObjectMarkerCoordinateFrame);
    aTransformRepository->SetTransform(pivotPointToMarkerTransformName, this->PivotPointToMarkerTransformMatrix);
    aTransformRepository->SetTransformPersistent(pivotPointToMarkerTransformName, true);
    aTransformRepository->SetTransformDate(pivotPointToMarkerTransformName, vtkIGSIOAccurateTimer::GetInstance()->GetDateAndTimeString().c_str());
    aTransformRepository->SetTransformError(pivotPointToMarkerTransformName, this->PivotCalibrationErrorMm);
  }
  else
  {
    LOG_DEBUG("Transform repository object is NULL, cannot save results into it");
  }

  this->PivotPointPosition_Reference[0] = pivotPoint_Reference[0];
  this->PivotPointPosition_Reference[1] = pivotPoint_Reference[1];
  this->PivotPointPosition_Reference[2] = pivotPoint_Reference[2];

  return status;
}

//----------------------------------------------------------------------------
igsioStatus vtkIGSIOPivotCalibrationAlgo::DoPivotCalibrationInternal(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, bool autoOrient, std::set<unsigned int>* outlierIndices, double pivotPoint_Marker[4], double pivotPoint_Reference[4], vtkMatrix4x4* pivotPointToMarkerTransformMatrix)
{
  if (!markerToTransformMatrixArray || markerToTransformMatrixArray->size() == 0)
  {
    return IGSIO_FAIL;
  }

  pivotPoint_Marker[0] = 0.0;
  pivotPoint_Marker[1] = 0.0;
  pivotPoint_Marker[2] = 0.0;
  pivotPoint_Marker[3] = 1.0;

  pivotPoint_Reference[0] = 0.0;
  pivotPoint_Reference[1] = 0.0;
  pivotPoint_Reference[2] = 0.0;
  pivotPoint_Reference[3] = 1.0;
  if (this->GetPivotPointPosition(markerToTransformMatrixArray, outlierIndices, pivotPoint_Marker, pivotPoint_Reference) != IGSIO_SUCCESS)
  {
    return IGSIO_FAIL;
  }

  // Get the result (tooltip to tool transform)
  double x = pivotPoint_Marker[0];
  double y = pivotPoint_Marker[1];
  double z = pivotPoint_Marker[2];

  pivotPointToMarkerTransformMatrix->SetElement(0, 3, x);
  pivotPointToMarkerTransformMatrix->SetElement(1, 3, y);
  pivotPointToMarkerTransformMatrix->SetElement(2, 3, z);

  // Compute tool orientation
  // Z axis: from the pivot point to the marker is on the -Z axis of the tool
  double pivotPointToMarkerTransformZ[3] = { x, y, z };
  vtkMath::Normalize(pivotPointToMarkerTransformZ);
  pivotPointToMarkerTransformMatrix->SetElement(0, 2, pivotPointToMarkerTransformZ[0]);
  pivotPointToMarkerTransformMatrix->SetElement(1, 2, pivotPointToMarkerTransformZ[1]);
  pivotPointToMarkerTransformMatrix->SetElement(2, 2, pivotPointToMarkerTransformZ[2]);

  // Y axis: orthogonal to tool's Z axis and the marker's X axis
  double pivotPointToMarkerTransformY[3] = { 0, 0, 0 };
  // Use the unitX vector as pivotPointToMarkerTransformX vector, unless unitX is parallel to pivotPointToMarkerTransformZ.
  // If unitX is parallel to pivotPointToMarkerTransformZ then use the unitY vector as pivotPointToMarkerTransformX.

  double unitX[3] = { 1, 0, 0 };
  double angle = acos(vtkMath::Dot(pivotPointToMarkerTransformZ, unitX));
  // Normalize between -pi/2 .. +pi/2
  if (angle > vtkMath::Pi() / 2)
  {
    angle -= vtkMath::Pi();
  }
  else if (angle < -vtkMath::Pi() / 2)
  {
    angle += vtkMath::Pi();
  }
  if (fabs(angle) * 180.0 / vtkMath::Pi() > 20.0)
  {
    // unitX is not parallel to pivotPointToMarkerTransformZ
    vtkMath::Cross(pivotPointToMarkerTransformZ, unitX, pivotPointToMarkerTransformY);
    LOG_DEBUG("Use unitX");
  }
  else
  {
    // unitX is parallel to pivotPointToMarkerTransformZ
    // use the unitY instead
    double unitY[3] = { 0, 1, 0 };
    vtkMath::Cross(pivotPointToMarkerTransformZ, unitY, pivotPointToMarkerTransformY);
    LOG_DEBUG("Use unitY");
  }
  vtkMath::Normalize(pivotPointToMarkerTransformY);
  pivotPointToMarkerTransformMatrix->SetElement(0, 1, pivotPointToMarkerTransformY[0]);
  pivotPointToMarkerTransformMatrix->SetElement(1, 1, pivotPointToMarkerTransformY[1]);
  pivotPointToMarkerTransformMatrix->SetElement(2, 1, pivotPointToMarkerTransformY[2]);

  // X axis: orthogonal to tool's Y axis and Z axis
  double pivotPointToMarkerTransformX[3] = { 0, 0, 0 };
  vtkMath::Cross(pivotPointToMarkerTransformY, pivotPointToMarkerTransformZ, pivotPointToMarkerTransformX);
  vtkMath::Normalize(pivotPointToMarkerTransformX);
  pivotPointToMarkerTransformMatrix->SetElement(0, 0, pivotPointToMarkerTransformX[0]);
  pivotPointToMarkerTransformMatrix->SetElement(1, 0, pivotPointToMarkerTransformX[1]);
  pivotPointToMarkerTransformMatrix->SetElement(2, 0, pivotPointToMarkerTransformX[2]);

  if (autoOrient)
  {
    this->UpdateShaftDirection(pivotPointToMarkerTransformMatrix); // Flip it if necessary
  }

  return IGSIO_SUCCESS;
}


//----------------------------------------------------------------------------
igsioStatus vtkIGSIOPivotCalibrationAlgo::DoSpinCalibration(vtkIGSIOTransformRepository* aTransformRepository/* = NULL*/, bool snapRotation/*=false*/, bool autoOrient/*=true*/)
{
  std::vector<vtkMatrix4x4*> markerToTransformMatrixArray = this->GetMarkerToReferenceTransformMatrixArray();
  vtkNew<vtkMatrix4x4> pivotPointToMarkerTransformMatrix;
  igsioStatus status = this->DoSpinCalibrationInternal(&markerToTransformMatrixArray, snapRotation, autoOrient, pivotPointToMarkerTransformMatrix, this->SpinCalibrationErrorMm);
  this->SetPivotPointToMarkerTransformMatrix(pivotPointToMarkerTransformMatrix);
  return status;
}

//----------------------------------------------------------------------------
igsioStatus vtkIGSIOPivotCalibrationAlgo::DoSpinCalibrationInternal(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, bool snapRotation, bool autoOrient, vtkMatrix4x4* toolTipToToolMatrix, double& error)
{
  if (!markerToTransformMatrixArray || markerToTransformMatrixArray->size() < 10)
  {
    this->SetErrorCode(CALIBRATION_NOT_ENOUGH_POINTS);
    return IGSIO_FAIL;
  }

  if (this->GetMaximumToolOrientationDifferenceDeg() < this->MinimumOrientationDifferenceDeg)
  {
    this->SetErrorCode(CALIBRATION_NOT_ENOUGH_VARIATION);
    return IGSIO_FAIL;
  }

  // Setup our system to find the axis of rotation
  unsigned int rows = 3, columns = 3;

  vnl_matrix<double> A(rows, columns, 0);

  vnl_matrix<double> I(3, 3, 0);
  I.set_identity();

  vnl_matrix<double> RI(rows, columns);

  std::vector< vtkMatrix4x4* >::const_iterator previt = markerToTransformMatrixArray->end();
  for (std::vector< vtkMatrix4x4* >::const_iterator it = markerToTransformMatrixArray->begin(); it != markerToTransformMatrixArray->end(); it++)
  {
    if (previt == markerToTransformMatrixArray->end())
    {
      previt = it;
      continue; // No comparison to make for the first matrix
    }

    vtkSmartPointer< vtkMatrix4x4 > itinverse = vtkSmartPointer< vtkMatrix4x4 >::New();
    vtkMatrix4x4::Invert((*it), itinverse);

    vtkSmartPointer< vtkMatrix4x4 > instRotation = vtkSmartPointer< vtkMatrix4x4 >::New();
    vtkMatrix4x4::Multiply4x4(itinverse, (*previt), instRotation);

    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        RI(i, j) = instRotation->GetElement(i, j);
      }
    }

    RI = RI - I;
    A = A + RI.transpose() * RI;

    previt = it;
  }

  // Setup the axes
  vnl_vector<double> shaftAxis_Shaft(columns, columns, SHAFT_AXIS);
  vnl_vector<double> orthogonalAxis_Shaft(columns, columns, ORTHOGONAL_AXIS);
  vnl_vector<double> backupAxis_Shaft(columns, columns, BACKUP_AXIS);

  // Find the eigenvector associated with the smallest eigenvalue
  // This is the best axis of rotation over all instantaneous rotations
  vnl_matrix<double> eigenvectors(columns, columns, 0);
  vnl_vector<double> eigenvalues(columns, 0);
  vnl_symmetric_eigensystem_compute(A, eigenvectors, eigenvalues);
  // Note: eigenvectors are ordered in increasing eigenvalue ( 0 = smallest, end = biggest )
  vnl_vector<double> shaftAxis_ToolTip(columns, 0);
  shaftAxis_ToolTip(0) = eigenvectors(0, 0);
  shaftAxis_ToolTip(1) = eigenvectors(1, 0);
  shaftAxis_ToolTip(2) = eigenvectors(2, 0);
  shaftAxis_ToolTip.normalize();

  // Snap the direction vector to be exactly aligned with one of the coordinate axes
  // This is if the sensor is known to be parallel to one of the axis, just not which one
  if (snapRotation)
  {
    int closestCoordinateAxis = element_product(shaftAxis_ToolTip, shaftAxis_ToolTip).arg_max();
    shaftAxis_ToolTip.fill(0);
    shaftAxis_ToolTip.put(closestCoordinateAxis, 1); // Doesn't matter the direction, will be sorted out later
  }

  //set the RMSE
  error = sqrt(eigenvalues(0) / markerToTransformMatrixArray->size());
  // Note: This error is the RMS distance from the ideal axis of rotation to the axis of rotation for each instantaneous rotation
  // This RMS distance can be computed to an angle in the following way: angle = arccos( 1 - SpinRMSE^2 / 2 )
  // Here we elect to return the RMS distance because this is the quantity that was actually minimized in the calculation

  // If the secondary axis 1 is parallel to the shaft axis in the tooltip frame, then use secondary axis 2
  vnl_vector<double> orthogonalAxis_ToolTip = this->ComputeSecondaryAxis(shaftAxis_ToolTip);
  // Do the registration find the appropriate rotation
  orthogonalAxis_ToolTip = orthogonalAxis_ToolTip - dot_product(orthogonalAxis_ToolTip, shaftAxis_ToolTip) * shaftAxis_ToolTip;
  orthogonalAxis_ToolTip.normalize();

  // Register X,Y,O points in the two coordinate frames (only spherical registration - since pure rotation)
  vnl_matrix<double> ToolTipPoints(3, 3, 0.0);
  vnl_matrix<double> ShaftPoints(3, 3, 0.0);

  ToolTipPoints.put(0, 0, shaftAxis_ToolTip(0));
  ToolTipPoints.put(0, 1, shaftAxis_ToolTip(1));
  ToolTipPoints.put(0, 2, shaftAxis_ToolTip(2));
  ToolTipPoints.put(1, 0, orthogonalAxis_ToolTip(0));
  ToolTipPoints.put(1, 1, orthogonalAxis_ToolTip(1));
  ToolTipPoints.put(1, 2, orthogonalAxis_ToolTip(2));
  ToolTipPoints.put(2, 0, 0);
  ToolTipPoints.put(2, 1, 0);
  ToolTipPoints.put(2, 2, 0);

  ShaftPoints.put(0, 0, shaftAxis_Shaft(0));
  ShaftPoints.put(0, 1, shaftAxis_Shaft(1));
  ShaftPoints.put(0, 2, shaftAxis_Shaft(2));
  ShaftPoints.put(1, 0, orthogonalAxis_Shaft(0));
  ShaftPoints.put(1, 1, orthogonalAxis_Shaft(1));
  ShaftPoints.put(1, 2, orthogonalAxis_Shaft(2));
  ShaftPoints.put(2, 0, 0);
  ShaftPoints.put(2, 1, 0);
  ShaftPoints.put(2, 2, 0);

  vnl_svd<double> ShaftToToolTipRegistrator(ShaftPoints.transpose() * ToolTipPoints);
  vnl_matrix<double> V = ShaftToToolTipRegistrator.V();
  vnl_matrix<double> U = ShaftToToolTipRegistrator.U();
  vnl_matrix<double> Rotation = V * U.transpose();

  // Make sure the determinant is positve (i.e. +1)
  double determinant = vnl_determinant(Rotation);
  if (determinant < 0)
  {
    // Switch the sign of the third column of V if the determinant is not +1
    // This is the recommended approach from Huang et al. 1987
    V.put(0, 2, -V.get(0, 2));
    V.put(1, 2, -V.get(1, 2));
    V.put(2, 2, -V.get(2, 2));
    Rotation = V * U.transpose();
  }

  // Set the elements of the output matrix
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      toolTipToToolMatrix->SetElement(i, j, Rotation[i][j]);
    }
  }
  if (autoOrient)
  {
    this->UpdateShaftDirection(toolTipToToolMatrix); // Flip it if necessary
  }

  this->SetErrorCode(CALIBRATION_SUCCESS);
  return IGSIO_SUCCESS;
}

//---------------------------------------------------------------------------
void vtkIGSIOPivotCalibrationAlgo::UpdateShaftDirection(vtkMatrix4x4* toolTipToToolMatrix)
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
void vtkIGSIOPivotCalibrationAlgo::FlipShaftDirection(vtkMatrix4x4* toolTipToToolMatrix)
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
void vtkIGSIOPivotCalibrationAlgo::GetToolTipToToolRotation(vtkMatrix4x4* toolTipToToolMatrix, vtkMatrix4x4* rotationMatrix)
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

//---------------------------------------------------------------------------
vnl_vector< double > vtkIGSIOPivotCalibrationAlgo::ComputeSecondaryAxis(vnl_vector< double > shaftAxis_ToolTip)
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


//-----------------------------------------------------------------------------
std::string vtkIGSIOPivotCalibrationAlgo::GetPivotPointToMarkerTranslationString(double aPrecision/*=3*/)
{
  if (this->PivotPointToMarkerTransformMatrix == NULL)
  {
    LOG_ERROR("Tooltip to tool transform is not initialized!");
    return "";
  }

  std::ostringstream s;
  s << std::fixed << std::setprecision(aPrecision)
    << this->PivotPointToMarkerTransformMatrix->GetElement(0, 3)
    << " x " << this->PivotPointToMarkerTransformMatrix->GetElement(1, 3)
    << " x " << this->PivotPointToMarkerTransformMatrix->GetElement(2, 3)
    << std::ends;

  return s.str();
}

//-----------------------------------------------------------------------------
igsioStatus vtkIGSIOPivotCalibrationAlgo::ReadConfiguration(vtkXMLDataElement* aConfig)
{
  XML_FIND_NESTED_ELEMENT_REQUIRED(pivotCalibrationElement, aConfig, "vtkIGSIOPivotCalibrationAlgo");
  XML_READ_CSTRING_ATTRIBUTE_REQUIRED(ObjectMarkerCoordinateFrame, pivotCalibrationElement);
  XML_READ_CSTRING_ATTRIBUTE_REQUIRED(ReferenceCoordinateFrame, pivotCalibrationElement);
  XML_READ_CSTRING_ATTRIBUTE_REQUIRED(ObjectPivotPointCoordinateFrame, pivotCalibrationElement);
  return IGSIO_SUCCESS;
}

//-----------------------------------------------------------------------------
void vtkIGSIOPivotCalibrationAlgo::ComputePivotCalibrationError()
{
  const std::vector<vtkMatrix4x4*> markerToTransformMatrixArray = this->GetMarkerToReferenceTransformMatrixArray();
  this->PivotCalibrationErrorMm = this->ComputePivotCalibrationError(&markerToTransformMatrixArray, &this->OutlierIndices, this->PivotPointPosition_Reference, this->PivotPointToMarkerTransformMatrix);
}

//-----------------------------------------------------------------------------
double vtkIGSIOPivotCalibrationAlgo::ComputePivotCalibrationError(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, std::set<unsigned int>* outlierIndices, double* pivotPoint_Reference, vtkMatrix4x4* pivotPointToMarkerTransformMatrix)
{
  vtkSmartPointer<vtkMatrix4x4> pivotPointToReferenceMatrix = vtkSmartPointer<vtkMatrix4x4>::New();

  // Compute the error for each sample as distance between the mean pivot point position and the pivot point position computed from each sample
  std::vector<double> errorValues;
  double currentPivotPoint_Reference[4] = { 0, 0, 0, 1 };
  unsigned int sampleIndex = 0;
  for (std::vector< vtkMatrix4x4* >::const_iterator markerToReferenceTransformIt = markerToTransformMatrixArray->begin();
    markerToReferenceTransformIt != markerToTransformMatrixArray->end(); ++markerToReferenceTransformIt, ++sampleIndex)
  {
    if (outlierIndices->find(sampleIndex) != outlierIndices->end())
    {
      // outlier, so skip from the error computation
      continue;
    }

    vtkMatrix4x4::Multiply4x4((*markerToReferenceTransformIt), pivotPointToMarkerTransformMatrix, pivotPointToReferenceMatrix);
    for (int i = 0; i < 3; i++)
    {
      currentPivotPoint_Reference[i] = pivotPointToReferenceMatrix->Element[i][3];
    }
    double errorValue = sqrt(vtkMath::Distance2BetweenPoints(currentPivotPoint_Reference, pivotPoint_Reference));
    errorValues.push_back(errorValue);
  }

  double mean = 0;
  double stdev = 0;
  igsioMath::ComputeMeanAndStdev(errorValues, mean, stdev);
  return mean;
}

//-----------------------------------------------------------------------------
int vtkIGSIOPivotCalibrationAlgo::GetNumberOfDetectedOutliers()
{
  return this->OutlierIndices.size();
}
