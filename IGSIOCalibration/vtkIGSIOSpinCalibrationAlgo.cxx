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
  this->PivotPointToMarkerTransformMatrix = NULL;
  this->SpinCalibrationErrorMm = -1.0;
}

//-----------------------------------------------------------------------------
vtkIGSIOSpinCalibrationAlgo::~vtkIGSIOSpinCalibrationAlgo()
{
  this->SetPivotPointToMarkerTransformMatrix(NULL);
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

//----------------------------------------------------------------------------
igsioStatus vtkIGSIOSpinCalibrationAlgo::DoCalibrationInternal(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, double& error)
{
  double pivotPoint_Marker[4] = { 0.0, 0.0, 0.0, 0.0 };
  double pivotPoint_Reference[4] = { 0.0, 0.0, 0.0, 0.0 };

  // For internal calibration used for bucket error detection, we don't care about orienting the tool corretcly, only the error.
  bool snapRotation = false;
  bool autoOrient = false;
  vtkNew<vtkMatrix4x4> pivotPointToTool;
  return this->DoSpinCalibrationInternal(markerToTransformMatrixArray, snapRotation, autoOrient, pivotPointToTool, error);
}

//----------------------------------------------------------------------------
igsioStatus vtkIGSIOSpinCalibrationAlgo::DoSpinCalibration(vtkIGSIOTransformRepository* aTransformRepository/* = NULL*/, bool snapRotation/*=false*/, bool autoOrient/*=true*/)
{
  if (!this->PivotPointToMarkerTransformMatrix)
  {
    vtkNew<vtkMatrix4x4> pivotPointToMarkerTransformMatrix;
    this->SetPivotPointToMarkerTransformMatrix(pivotPointToMarkerTransformMatrix);
  }

  std::vector<vtkMatrix4x4*> markerToTransformMatrixArray = this->GetAllMarkerToReferenceMatrices();
  return this->DoSpinCalibrationInternal(&markerToTransformMatrixArray, snapRotation, autoOrient, this->PivotPointToMarkerTransformMatrix, this->SpinCalibrationErrorMm);
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
