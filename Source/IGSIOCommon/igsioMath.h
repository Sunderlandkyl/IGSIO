/*=Plus=header=begin======================================================
  Program: Plus
  Copyright (c) Laboratory for Percutaneous Surgery. All rights reserved.
  See License.txt for details.
=========================================================Plus=header=end*/

#ifndef __IGSIOMATH_H
#define __IGSIOMATH_H

//#include "PlusConfigure.h"
#include "vtkigsiocommon_export.h"
#include "igsioCommon.h"

#include <vector>

#include "vtkMath.h"

class vtkMatrix4x4;
class vtkTransform;

#if _MSC_VER == 1600 // VS 2010
namespace std
{
  double VTKIGSIOCOMMON_EXPORT round(double arg);
}
#endif

/*!
  \class igsioMath
  \brief A utility class that contains static functions for various useful commonly used computations
  \ingroup PlusLibCommon
*/
class VTKIGSIOCOMMON_EXPORT igsioMath
{
public:

  /*! Returns the Euclidean distance between two 4x4 homogeneous transformation matrix */
  static double GetPositionDifference(vtkMatrix4x4* aMatrix, vtkMatrix4x4* bMatrix); 

  /*! Returns the orientation difference in degrees between two 4x4 homogeneous transformation matrix, in degrees. */
  static double GetOrientationDifference(vtkMatrix4x4* aMatrix, vtkMatrix4x4* bMatrix); 

protected:
  igsioMath(); 
  ~igsioMath();

private: 
  igsioMath(igsioMath const&);
  igsioMath& operator=(igsioMath const&);
};

#endif 