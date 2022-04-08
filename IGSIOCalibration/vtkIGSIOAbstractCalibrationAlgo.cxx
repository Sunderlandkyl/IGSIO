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

//-----------------------------------------------------------------------------
vtkIGSIOAbstractCalibrationAlgo::vtkIGSIOAbstractCalibrationAlgo()
{
  /*this->PivotPointToMarkerTransformMatrix = NULL;*/
  this->PreviousMarkerToReferenceTransformMatrix = vtkMatrix4x4::New();

  this->ObjectMarkerCoordinateFrame = NULL;
  this->ReferenceCoordinateFrame = NULL;
  /*this->ObjectPivotPointCoordinateFrame = NULL;*/

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
  this->GetMarkerToReferenceTransformMatrixArray(this->MarkerToReferenceTransformMatrixBuckets.size() - 1, &latestBucket);

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
std::vector<vtkMatrix4x4*> vtkIGSIOAbstractCalibrationAlgo::GetMarkerToReferenceTransformMatrixArray()
{
  std::vector<vtkMatrix4x4*> markerToTransformMatrixArray;
  for (int i = 0; i < this->MarkerToReferenceTransformMatrixBuckets.size(); ++i)
  {
    this->GetMarkerToReferenceTransformMatrixArray(i, &markerToTransformMatrixArray);
  }
  return markerToTransformMatrixArray;
}

//----------------------------------------------------------------------------
void vtkIGSIOAbstractCalibrationAlgo::GetMarkerToReferenceTransformMatrixArray(std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray)
{
  for (int i = 0; i < this->MarkerToReferenceTransformMatrixBuckets.size(); ++i)
  {
    this->GetMarkerToReferenceTransformMatrixArray(i, markerToTransformMatrixArray);
  }
}

//----------------------------------------------------------------------------
std::vector<vtkMatrix4x4*> vtkIGSIOAbstractCalibrationAlgo::GetMarkerToReferenceTransformMatrixArray(int bucketIndex)
{
  std::vector<vtkMatrix4x4*> markerToTransformMatrixArray;
  if (this->MarkerToReferenceTransformMatrixBuckets.size() > bucketIndex)
  {
    this->GetMarkerToReferenceTransformMatrixArray(bucketIndex, &markerToTransformMatrixArray);
  }
  return markerToTransformMatrixArray;
}

//----------------------------------------------------------------------------
void vtkIGSIOAbstractCalibrationAlgo::GetMarkerToReferenceTransformMatrixArray(int bucketIndex, std::vector<vtkMatrix4x4*>* matrixArray)
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
  this->GetMarkerToReferenceTransformMatrixArray(&toolToReferenceMatrices);
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
