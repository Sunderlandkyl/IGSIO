/*=IGSIO=header=begin======================================================
Program: IGSIO
Copyright (c) Laboratory for Percutaneous Surgery. All rights reserved.
See License.txt for details.
=========================================================IGSIO=header=end*/

#ifndef __vtkIGSIOSpinCalibrationAlgo_h
#define __vtkIGSIOSpinCalibrationAlgo_h

// Local includes
#include "igsioConfigure.h"
#include "igsioCommon.h"
#include "vtkIGSIOAbstractCalibrationAlgo.h"
#include "vtkIGSIOCalibrationExport.h"

// VTK includes
#include <vtkCommand.h>
#include <vtkObject.h>
#include <vtkMatrix4x4.h>

class vtkIGSIOTransformRepository;
class vtkXMLDataElement;

//-----------------------------------------------------------------------------

/*!
  \class vtkIGSIOSpinCalibrationAlgo
  \brief Spin calibration algorithm to calibrate a stylus. It determines the direction of the stylus shaft.

  \ingroup igsioCalibrationAlgorithm
*/
class vtkIGSIOCalibrationExport vtkIGSIOSpinCalibrationAlgo : public vtkIGSIOAbstractCalibrationAlgo
{
public:
  vtkTypeMacro(vtkIGSIOSpinCalibrationAlgo, vtkIGSIOAbstractCalibrationAlgo);
  static vtkIGSIOSpinCalibrationAlgo* New();

  /// Read configuration
  /// \param aConfig Root element of the device set configuration
  igsioStatus ReadConfiguration(vtkXMLDataElement* aConfig) override;

  /// Calibrate (call the minimizer and set the result)
  /// \param aTransformRepository Transform repository to save the results into
  igsioStatus DoSpinCalibration(vtkIGSIOTransformRepository* aTransformRepository = NULL, bool snapRotation = false, bool autoOrient = true);

  /// Mean error of the pivot calibration result in mm
  vtkGetMacro(SpinCalibrationErrorMm, double);

protected:
  vtkIGSIOSpinCalibrationAlgo();
  virtual ~vtkIGSIOSpinCalibrationAlgo();

protected:
  igsioStatus DoCalibrationInternal(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, double& error) override;
  igsioStatus DoSpinCalibrationInternal(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, bool snapRotation, bool autoOrient, vtkMatrix4x4* pivotPointToMarkerTransformMatrix, double& error);

protected:
  double SpinCalibrationErrorMm;

  class vtkInternal;
  vtkInternal* Internal;
};

#endif
