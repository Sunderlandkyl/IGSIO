/*=Plus=header=begin======================================================
  Program: Plus
  Copyright (c) Laboratory for Percutaneous Surgery. All rights reserved.
  See License.txt for details.
=========================================================Plus=header=end*/

//#include "PlusConfigure.h"

#include "igsioMath.h"

//#include "vnl/vnl_sparse_matrix.h"
//#include "vnl/vnl_sparse_matrix_linear_system.h"  
//#include "vnl/algo/vnl_lsqr.h"  
//#include "vnl/vnl_cross.h"  

#include "vtkMath.h"
#include "vtkTransform.h"

#define MINIMUM_NUMBER_OF_CALIBRATION_EQUATIONS 8

//----------------------------------------------------------------------------
igsioMath::igsioMath()
{

}

//----------------------------------------------------------------------------
igsioMath::~igsioMath()
{

}

//----------------------------------------------------------------------------
double igsioMath::GetPositionDifference(vtkMatrix4x4* aMatrix, vtkMatrix4x4* bMatrix)
{
  //**LOG_TRACE("igsioMath::GetPositionDifference"); 
  vtkSmartPointer<vtkTransform> aTransform = vtkSmartPointer<vtkTransform>::New(); 
  aTransform->SetMatrix(aMatrix); 

  vtkSmartPointer<vtkTransform> bTransform = vtkSmartPointer<vtkTransform>::New(); 
  bTransform->SetMatrix(bMatrix); 

  double ax = aTransform->GetPosition()[0]; 
  double ay = aTransform->GetPosition()[1]; 
  double az = aTransform->GetPosition()[2]; 

  double bx = bTransform->GetPosition()[0]; 
  double by = bTransform->GetPosition()[1]; 
  double bz = bTransform->GetPosition()[2]; 

  // Euclidean distance
  double distance = sqrt( pow(ax-bx,2) + pow(ay-by,2) + pow(az-bz,2) ); 

  return distance; 
}

//----------------------------------------------------------------------------
double igsioMath::GetOrientationDifference(vtkMatrix4x4* aMatrix, vtkMatrix4x4* bMatrix)
{
  vtkSmartPointer<vtkMatrix4x4> diffMatrix = vtkSmartPointer<vtkMatrix4x4>::New(); 
  vtkSmartPointer<vtkMatrix4x4> invBmatrix = vtkSmartPointer<vtkMatrix4x4>::New(); 

  vtkMatrix4x4::Invert(bMatrix, invBmatrix);  

  vtkMatrix4x4::Multiply4x4(aMatrix, invBmatrix, diffMatrix); 

  vtkSmartPointer<vtkTransform> diffTransform = vtkSmartPointer<vtkTransform>::New(); 
  diffTransform->SetMatrix(diffMatrix); 

  double angleDiff_rad= vtkMath::RadiansFromDegrees(diffTransform->GetOrientationWXYZ()[0]);

  double normalizedAngleDiff_rad = atan2( sin(angleDiff_rad), cos(angleDiff_rad) ); // normalize angle to domain -pi, pi 

  return vtkMath::DegreesFromRadians(normalizedAngleDiff_rad);
}
