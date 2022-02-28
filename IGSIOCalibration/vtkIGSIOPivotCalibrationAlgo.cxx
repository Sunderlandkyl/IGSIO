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
#include <vnl/vnl_vector.h>

vtkStandardNewMacro(vtkIGSIOPivotCalibrationAlgo);

//-----------------------------------------------------------------------------
vtkIGSIOPivotCalibrationAlgo::vtkIGSIOPivotCalibrationAlgo()
{
  this->PivotPointToMarkerTransformMatrix = NULL;
  this->CalibrationError = -1.0;
  this->ObjectMarkerCoordinateFrame = NULL;
  this->ReferenceCoordinateFrame = NULL;
  this->ObjectPivotPointCoordinateFrame = NULL;

  this->PivotPointPosition_Reference[0] = 0.0;
  this->PivotPointPosition_Reference[1] = 0.0;
  this->PivotPointPosition_Reference[2] = 0.0;
  this->PivotPointPosition_Reference[3] = 1.0;

  this->MinimumOrientationDifferenceDeg = 15.0;

  this->ErrorCode = CALIBRATION_NOT_STARTED;
  this->CalibrationPoseBucketSize = -1;
  this->MaximumBucketError = 3.0;
}

//-----------------------------------------------------------------------------
vtkIGSIOPivotCalibrationAlgo::~vtkIGSIOPivotCalibrationAlgo()
{
  this->SetPivotPointToMarkerTransformMatrix(NULL);
  this->RemoveAllCalibrationPoints();
}

//-----------------------------------------------------------------------------
void vtkIGSIOPivotCalibrationAlgo::RemoveAllCalibrationPoints()
{
  this->MarkerToReferenceTransformMatrixBuckets.clear();
  this->ErrorCode = CALIBRATION_NOT_STARTED;
  this->OutlierIndices.clear();
}

//----------------------------------------------------------------------------
igsioStatus vtkIGSIOPivotCalibrationAlgo::InsertNextCalibrationPoint(vtkMatrix4x4* aMarkerToReferenceTransformMatrix)
{
  MarkerToReferenceTransformMatrixBucket* currentBucket = NULL;
  if (this->MarkerToReferenceTransformMatrixBuckets.size() > 0)
  {
    currentBucket = &this->MarkerToReferenceTransformMatrixBuckets[this->MarkerToReferenceTransformMatrixBuckets.size() - 1];
  }
  if (!currentBucket)
  {
    this->MarkerToReferenceTransformMatrixBuckets.push_back(MarkerToReferenceTransformMatrixBucket());
    currentBucket = &this->MarkerToReferenceTransformMatrixBuckets[this->MarkerToReferenceTransformMatrixBuckets.size() - 1];
  }

  vtkNew<vtkMatrix4x4> markerToReferenceTransformMatrixCopy;
  markerToReferenceTransformMatrixCopy->DeepCopy(aMarkerToReferenceTransformMatrix);
  currentBucket->MarkerToReferenceCalibrationPoints.push_back(markerToReferenceTransformMatrixCopy);
  if (this->CalibrationPoseBucketSize > 0 && currentBucket->MarkerToReferenceCalibrationPoints.size() >= this->CalibrationPoseBucketSize)
  {
    this->CleanInputBuffer();

    // The current bucket is full. We will create a new one.
    this->MarkerToReferenceTransformMatrixBuckets.push_back(MarkerToReferenceTransformMatrixBucket());
  }
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
  igsioStatus status = this->DoPivotCalibrationInternal(&latestBucket, pivotPoint_Marker, pivotPoint_Reference, pivotPointToMarkerTransformMatrix);
  double meanError = this->ComputeCalibrationError(&latestBucket, pivotPoint_Reference, pivotPointToMarkerTransformMatrix);
  if (status != IGSIO_SUCCESS || meanError > this->MaximumBucketError)
  {
    // The latest pivot calibration bucket does not contain pivoting.
    // Discard it and all earlier buckets.
    this->MarkerToReferenceTransformMatrixBuckets.clear();
    this->ErrorCode = CALIBRATION_HIGH_ERROR;
  }

  if (this->MarkerToReferenceTransformMatrixBuckets.size() > this->MaximumNumberOfBuckets)
  {
    this->MarkerToReferenceTransformMatrixBuckets.pop_front();
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
igsioStatus vtkIGSIOPivotCalibrationAlgo::GetPivotPointPosition(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, double* pivotPoint_Marker, double* pivotPoint_Reference)
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

  this->OutlierIndices.clear();
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
        this->OutlierIndices.insert(sampleIndex);
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
  unsigned int currentRow;
  vtkMatrix4x4* referenceOrientationMatrix = toolToReferenceMatrices.front();
  std::vector<vtkMatrix4x4*>::const_iterator it;
  for (currentRow = 0, it = toolToReferenceMatrices.begin(); it != matricesEnd; it++, currentRow += 3)
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
igsioStatus vtkIGSIOPivotCalibrationAlgo::DoPivotCalibration(vtkIGSIOTransformRepository* aTransformRepository/* = NULL*/)
{
  if (this->MarkerToReferenceTransformMatrixBuckets.empty())
  {
    this->ErrorCode = CALIBRATION_NOT_ENOUGH_POINTS;
    return IGSIO_FAIL;
  }

  if (this->GetMaximumToolOrientationDifferenceDeg() < this->MinimumOrientationDifferenceDeg)
  {
    this->ErrorCode = CALIBRATION_NOT_ENOUGH_VARIATION;
    return IGSIO_FAIL;
  }

  vtkNew<vtkMatrix4x4> pivotPointToMarkerTransformMatrix;
  std::vector<vtkMatrix4x4*> markerToTransformMatrixArray = this->GetMarkerToReferenceTransformMatrixArray();

  double pivotPoint_Marker[4] = { 0.0, 0.0, 0.0, 1.0 };
  double pivotPoint_Reference[4] = { 0.0, 0.0, 0.0, 1.0 };
  igsioStatus status = this->DoPivotCalibrationInternal(&markerToTransformMatrixArray, pivotPoint_Marker, pivotPoint_Reference, pivotPointToMarkerTransformMatrix);
  if (status == IGSIO_SUCCESS)
  {
    this->ErrorCode = CALIBRATION_SUCCESS;
  }
  else
  {
    this->ErrorCode = CALIBRATION_FAIL;
  }

  this->SetPivotPointToMarkerTransformMatrix(pivotPointToMarkerTransformMatrix);

  // Save result
  if (aTransformRepository)
  {
    igsioTransformName pivotPointToMarkerTransformName(this->ObjectPivotPointCoordinateFrame, this->ObjectMarkerCoordinateFrame);
    aTransformRepository->SetTransform(pivotPointToMarkerTransformName, this->PivotPointToMarkerTransformMatrix);
    aTransformRepository->SetTransformPersistent(pivotPointToMarkerTransformName, true);
    aTransformRepository->SetTransformDate(pivotPointToMarkerTransformName, vtkIGSIOAccurateTimer::GetInstance()->GetDateAndTimeString().c_str());
    aTransformRepository->SetTransformError(pivotPointToMarkerTransformName, this->CalibrationError);
  }
  else
  {
    LOG_DEBUG("Transform repository object is NULL, cannot save results into it");
  }

  this->PivotPointPosition_Reference[0] = pivotPoint_Reference[0];
  this->PivotPointPosition_Reference[1] = pivotPoint_Reference[1];
  this->PivotPointPosition_Reference[2] = pivotPoint_Reference[2];

  this->ComputeCalibrationError();

  return status;
}

//----------------------------------------------------------------------------
igsioStatus vtkIGSIOPivotCalibrationAlgo::DoPivotCalibrationInternal(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, double pivotPoint_Marker[4], double pivotPoint_Reference[4], vtkMatrix4x4* pivotPointToMarkerTransformMatrix)
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
  if (this->GetPivotPointPosition(markerToTransformMatrixArray, pivotPoint_Marker, pivotPoint_Reference) != IGSIO_SUCCESS)
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

  return IGSIO_SUCCESS;
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
void vtkIGSIOPivotCalibrationAlgo::ComputeCalibrationError()
{
  const std::vector<vtkMatrix4x4*> markerToTransformMatrixArray = this->GetMarkerToReferenceTransformMatrixArray();
  this->CalibrationError = this->ComputeCalibrationError(&markerToTransformMatrixArray, this->PivotPointPosition_Reference, this->PivotPointToMarkerTransformMatrix);
}

//-----------------------------------------------------------------------------
double vtkIGSIOPivotCalibrationAlgo::ComputeCalibrationError(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, double* pivotPoint_Reference, vtkMatrix4x4* pivotPointToMarkerTransformMatrix)
{
  vtkSmartPointer<vtkMatrix4x4> pivotPointToReferenceMatrix = vtkSmartPointer<vtkMatrix4x4>::New();

  // Compute the error for each sample as distance between the mean pivot point position and the pivot point position computed from each sample
  std::vector<double> errorValues;
  double currentPivotPoint_Reference[4] = { 0, 0, 0, 1 };
  unsigned int sampleIndex = 0;
  for (std::vector< vtkMatrix4x4* >::const_iterator markerToReferenceTransformIt = markerToTransformMatrixArray->begin();
    markerToReferenceTransformIt != markerToTransformMatrixArray->end(); ++markerToReferenceTransformIt, ++sampleIndex)
  {
    // TODO: Discard outliers
    //if (this->OutlierIndices.find(sampleIndex) != this->OutlierIndices.end())
    //{
    //  // outlier, so skip from the error computation
    //  continue;
    //}

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
