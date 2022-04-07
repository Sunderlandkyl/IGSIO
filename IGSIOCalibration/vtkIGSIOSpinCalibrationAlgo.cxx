/*=IGSIO=header=begin======================================================
  Program: IGSIO
  Copyright (c) Laboratory for Percutaneous Surgery. All rights reserved.
  See License.txt for details.
=========================================================IGSIO=header=end*/

// IGSIO includes
#include "igsioConfigure.h"
#include "igsioMath.h"
#include "vtkIGSIOSpinCalibrationAlgo.h"
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

vtkStandardNewMacro(vtkIGSIOSpinCalibrationAlgo);

//-----------------------------------------------------------------------------
vtkIGSIOSpinCalibrationAlgo::vtkIGSIOSpinCalibrationAlgo()
{
  this->SpinCalibrationErrorMm = -1.0;
}

//-----------------------------------------------------------------------------
vtkIGSIOSpinCalibrationAlgo::~vtkIGSIOSpinCalibrationAlgo()
{
  this->SetPivotPointToMarkerTransformMatrix(NULL);
}

//-----------------------------------------------------------------------------
void vtkIGSIOSpinCalibrationAlgo::RemoveAllCalibrationPoints()
{
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
igsioStatus vtkIGSIOSpinCalibrationAlgo::GetPivotPointPosition(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, std::set<unsigned int>* outlierIndices, double* pivotPoint_Marker, double* pivotPoint_Reference)
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
    LOG_ERROR("vtkIGSIOSpinCalibrationAlgo failed: LSQRMinimize error");
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
std::vector<vtkMatrix4x4*> vtkIGSIOSpinCalibrationAlgo::GetMarkerToReferenceTransformMatrixArray()
{
  std::vector<vtkMatrix4x4*> markerToTransformMatrixArray;
  for (int i = 0; i < this->MarkerToReferenceTransformMatrixBuckets.size(); ++i)
  {
    this->GetMarkerToReferenceTransformMatrixArray(i, &markerToTransformMatrixArray);
  }
  return markerToTransformMatrixArray;
}

//----------------------------------------------------------------------------
void vtkIGSIOSpinCalibrationAlgo::GetMarkerToReferenceTransformMatrixArray(std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray)
{
  for (int i = 0; i < this->MarkerToReferenceTransformMatrixBuckets.size(); ++i)
  {
    this->GetMarkerToReferenceTransformMatrixArray(i, markerToTransformMatrixArray);
  }
}

//----------------------------------------------------------------------------
std::vector<vtkMatrix4x4*> vtkIGSIOSpinCalibrationAlgo::GetMarkerToReferenceTransformMatrixArray(int bucketIndex)
{
  std::vector<vtkMatrix4x4*> markerToTransformMatrixArray;
  if (this->MarkerToReferenceTransformMatrixBuckets.size() > bucketIndex)
  {
    this->GetMarkerToReferenceTransformMatrixArray(bucketIndex, &markerToTransformMatrixArray);
  }
  return markerToTransformMatrixArray;
}

//----------------------------------------------------------------------------
void vtkIGSIOSpinCalibrationAlgo::GetMarkerToReferenceTransformMatrixArray(int bucketIndex, std::vector<vtkMatrix4x4*>* matrixArray)
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
igsioStatus vtkIGSIOSpinCalibrationAlgo::DoCalibrationInternal(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, double& error)
{
  double pivotPoint_Marker[4] = { 0.0, 0.0, 0.0, 0.0 };
  double pivotPoint_Reference[4] = { 0.0, 0.0, 0.0, 0.0 };
  vtkNew<vtkMatrix4x4> pivotPointToTool;
  bool snapRotation = false;
  bool autoOrient = true;
  return this->DoSpinCalibrationInternal(markerToTransformMatrixArray, snapRotation, autoOrient, pivotPointToTool, error);
}

//----------------------------------------------------------------------------
igsioStatus vtkIGSIOSpinCalibrationAlgo::DoSpinCalibration(vtkIGSIOTransformRepository* aTransformRepository/* = NULL*/, bool snapRotation/*=false*/, bool autoOrient/*=true*/)
{
  std::vector<vtkMatrix4x4*> markerToTransformMatrixArray = this->GetMarkerToReferenceTransformMatrixArray();
  vtkNew<vtkMatrix4x4> pivotPointToMarkerTransformMatrix;
  igsioStatus status = this->DoSpinCalibrationInternal(&markerToTransformMatrixArray, snapRotation, autoOrient, pivotPointToMarkerTransformMatrix, this->SpinCalibrationErrorMm);
  this->SetPivotPointToMarkerTransformMatrix(pivotPointToMarkerTransformMatrix);
  return status;
}

//----------------------------------------------------------------------------
igsioStatus vtkIGSIOSpinCalibrationAlgo::DoSpinCalibrationInternal(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, bool snapRotation, bool autoOrient, vtkMatrix4x4* toolTipToToolMatrix, double& error)
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
void vtkIGSIOSpinCalibrationAlgo::UpdateShaftDirection(vtkMatrix4x4* toolTipToToolMatrix)
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
void vtkIGSIOSpinCalibrationAlgo::FlipShaftDirection(vtkMatrix4x4* toolTipToToolMatrix)
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
void vtkIGSIOSpinCalibrationAlgo::GetToolTipToToolRotation(vtkMatrix4x4* toolTipToToolMatrix, vtkMatrix4x4* rotationMatrix)
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
vnl_vector< double > vtkIGSIOSpinCalibrationAlgo::ComputeSecondaryAxis(vnl_vector< double > shaftAxis_ToolTip)
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
std::string vtkIGSIOSpinCalibrationAlgo::GetPivotPointToMarkerTranslationString(double aPrecision/*=3*/)
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
igsioStatus vtkIGSIOSpinCalibrationAlgo::ReadConfiguration(vtkXMLDataElement* aConfig)
{
  XML_FIND_NESTED_ELEMENT_REQUIRED(pivotCalibrationElement, aConfig, "vtkIGSIOSpinCalibrationAlgo");
  XML_READ_CSTRING_ATTRIBUTE_REQUIRED(ObjectMarkerCoordinateFrame, pivotCalibrationElement);
  XML_READ_CSTRING_ATTRIBUTE_REQUIRED(ReferenceCoordinateFrame, pivotCalibrationElement);
  XML_READ_CSTRING_ATTRIBUTE_REQUIRED(ObjectPivotPointCoordinateFrame, pivotCalibrationElement);
  return IGSIO_SUCCESS;
}
