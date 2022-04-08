/*=IGSIO=header=begin======================================================
Program: IGSIO
Copyright (c) Laboratory for Percutaneous Surgery. All rights reserved.
See License.txt for details.
=========================================================IGSIO=header=end*/

#ifndef __vtkIGSIOAbstractCalibrationAlgo_h
#define __vtkIGSIOAbstractCalibrationAlgo_h

// Local includes
#include "igsioConfigure.h"
#include "igsioCommon.h"
#include "vtkIGSIOCalibrationExport.h"

// VTK includes
#include <vtkCommand.h>
#include <vtkObject.h>
#include <vtkMatrix4x4.h>

// STL includes
#include <set>
#include <deque>

class vtkIGSIOTransformRepository;
class vtkXMLDataElement;

//-----------------------------------------------------------------------------

/*!
  \class vtkIGSIOAbstractCalibrationAlgo
  \brief Pivot calibration algorithm to calibrate a stylus. It determines the pose of the stylus tip relative to the marker attached to the stylus.

  The stylus tip position is computed by robust LSQR method, which detects and ignores outliers (that have much larger reprojection error than other points).

  The stylus pose is computed assuming that the marker is attached on the center of one of the stylus axes, which is often a good approximation.
  The axis that points towards the marker is the PivotPoint coordinate system's -Z axis (so that points in front of the stylus have positive Z coordinates
  in the PivotPoint coordinate system). The X axis of the PivotPoint coordinate system is
  aligned with the marker coordinate system's X axis (unless the Z axis of the PivotPoint coordinate system is parallel with the marker coordinate
  system's X axis; in this case the X axis of the PivotPoint coordinate system is aligned with the marker coordinate system's Y axis). The Y axis
  of the PivotPoint coordinate system is chosen to be the cross product of the Z and X axes.

  The method detects outlier points (points that have larger than 3x error than the standard deviation) and ignores them when computing the pivot point
  coordinates and the calibration error.

  \ingroup igsioCalibrationAlgorithm
*/
class vtkIGSIOCalibrationExport vtkIGSIOAbstractCalibrationAlgo : public vtkObject
{
public:
  vtkTypeMacro(vtkIGSIOAbstractCalibrationAlgo, vtkObject);
  /*static vtkIGSIOAbstractCalibrationAlgo* New();*/

  enum Events
  {
    InputTransformAddedEvent = vtkCommand::UserEvent + 173,
  };

  enum CalibrationErrorCodes
  {
    CALIBRATION_FAIL = IGSIO_FAIL,
    CALIBRATION_SUCCESS = IGSIO_SUCCESS,
    CALIBRATION_NOT_STARTED,
    CALIBRATION_NOT_ENOUGH_POINTS,
    CALIBRATION_NOT_ENOUGH_VARIATION,
    CALIBRATION_HIGH_ERROR,
  };

  /*!
  * Read configuration
  * \param aConfig Root element of the device set configuration
  */
  virtual igsioStatus ReadConfiguration(vtkXMLDataElement* aConfig) = 0;

  /*!
    Remove all previously inserted calibration points.
    Call this method to get rid of previously added calibration points
    before starting a new calibration.
  */
  virtual void RemoveAllCalibrationPoints();

  /*!
    Insert acquired point to calibration point list
    \param aMarkerToReferenceTransformMatrix New calibration point (tool to reference transform)
  */
  virtual igsioStatus InsertNextCalibrationPoint(vtkMatrix4x4* aMarkerToReferenceTransformMatrix);

  /*!
    Get the number of outlier points. It is recommended to display a warning to the user
    if the percentage of outliers vs total number of points is larger than a few percent.
  */
  virtual int GetNumberOfDetectedOutliers();

  virtual int GetNumberOfCalibrationPoints();

  // Computes the maximum orientation difference in degrees between the first tool transformation
  // and all the others. Used for determining if there was enough variation in the input data.
  virtual double GetMaximumToolOrientationDifferenceDeg();

  //@{
  /// Output calibration tip to marker position
  /*vtkGetObjectMacro(PivotPointToMarkerTransformMatrix, vtkMatrix4x4);*/
  /*vtkGetVector3Macro(PivotPointPosition_Reference, double);*/
  //@}

  /// Name of the object marker coordinate frame (eg. Stylus)
  vtkGetStringMacro(ObjectMarkerCoordinateFrame);

  /// Name of the reference coordinate frame (eg. Reference)
  vtkGetStringMacro(ReferenceCoordinateFrame);

  /// Name of the object pivot point coordinate frame (eg. StylusTip)
  vtkGetStringMacro(ObjectPivotPointCoordinateFrame);

  /// Error code indicating what went wrong with the calibration.
  vtkGetMacro(ErrorCode, int);

  //@{
  /// Required minimum amount of variation within the recorded poses
  vtkSetMacro(MinimumOrientationDifferenceDeg, double);
  vtkGetMacro(MinimumOrientationDifferenceDeg, double);
  //@}

  //@{
  /// Required minimum amount of variation in position from the previous position in order for a transform to be accepted (0 degrees by default).
  vtkGetMacro(PositionDifferenceThresholdMm, double);
  vtkSetMacro(PositionDifferenceThresholdMm, double);
  //@}

  //@{
  /// Required minimum amount of variation in position from the previous position in order for a transform to be accepted (0 degrees by default).
  vtkGetMacro(OrientationDifferenceThresholdDegrees, double);
  vtkSetMacro(OrientationDifferenceThresholdDegrees, double);
  //@}

  //@{
  /// The number of input points to be stored in each bucket
  vtkGetMacro(PoseBucketSize, int);
  vtkSetMacro(PoseBucketSize, int);
  //@}

  //@{
  /// The maximum number of buckets to be stored. If the number of buckets is exceeded, the oldest will be discarded.
  vtkGetMacro(MaximumNumberOfPoseBuckets, int);
  vtkSetMacro(MaximumNumberOfPoseBuckets, int);
  //@}

  //@{
  /// The accepted amount of error within a single bucket. If the error in the latest bucket exceeds this amount, then all of the buckets will be discarded.
  vtkGetMacro(MaximumPoseBucketError, double);
  vtkSetMacro(MaximumPoseBucketError, double);
  //@}

  //@{
  /// Output calibration tip to marker position
  vtkGetObjectMacro(PivotPointToMarkerTransformMatrix, vtkMatrix4x4);
  vtkSetObjectMacro(PivotPointToMarkerTransformMatrix, vtkMatrix4x4);
  vtkGetVector3Macro(PivotPointPosition_Reference, double);
  //@}

protected:
  /*vtkSetObjectMacro(PivotPointToMarkerTransformMatrix, vtkMatrix4x4);*/
  vtkSetStringMacro(ObjectMarkerCoordinateFrame);
  vtkSetStringMacro(ReferenceCoordinateFrame);
  vtkSetStringMacro(ObjectPivotPointCoordinateFrame);
  vtkSetMacro(ErrorCode, int);

protected:
  vtkIGSIOAbstractCalibrationAlgo();
  virtual ~vtkIGSIOAbstractCalibrationAlgo();

protected:

  virtual std::vector<vtkMatrix4x4*> GetMarkerToReferenceTransformMatrixArray();
  virtual void GetMarkerToReferenceTransformMatrixArray(std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray);
  virtual std::vector<vtkMatrix4x4*> GetMarkerToReferenceTransformMatrixArray(int bucket);
  virtual void GetMarkerToReferenceTransformMatrixArray(int bucket, std::vector<vtkMatrix4x4*>* matrixArray);

  /// Removes all data if the latest bucket is not within the acceptable error threshold.
  igsioStatus CleanInputBuffer();

  virtual igsioStatus DoCalibrationInternal(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, double& error) = 0;

  /*virtual void GetToolTipToToolRotation(vtkMatrix4x4* toolTipToToolMatrix, vtkMatrix4x4* rotationMatrix);*/

  // Returns the orientation difference in degrees between two 4x4 homogeneous transformation matrix, in degrees.
  virtual double GetOrientationDifferenceDeg(vtkMatrix4x4* aMatrix, vtkMatrix4x4* bMatrix);

  // Verify whether the tool's shaft is in the same direction as the ToolTip to Tool vector.
  // Rotate the ToolTip coordinate frame by 180 degrees about the secondary axis to make the
  // shaft in the same direction as the ToolTip to Tool vector, if this is not already the case.
  void UpdateShaftDirection(vtkMatrix4x4* toolTipToToolMatrix);

  // Flip the direction of the shaft axis
  void FlipShaftDirection(vtkMatrix4x4* toolTipToToolMatrix);

  void GetToolTipToToolRotation(vtkMatrix4x4* toolTipToToolMatrix, vtkMatrix4x4* rotationMatrix);

  // Helper method to compute the secondary axis, given a shaft axis
  vnl_vector< double > ComputeSecondaryAxis(vnl_vector< double > shaftAxis_ToolTip);

protected:
  vtkMatrix4x4* PivotPointToMarkerTransformMatrix;
  double        PivotPointPosition_Reference[4];

  vtkMatrix4x4* PreviousMarkerToReferenceTransformMatrix;

  char*                     ObjectMarkerCoordinateFrame;
  char*                     ReferenceCoordinateFrame;
  char*                     ObjectPivotPointCoordinateFrame;

  /*! List of outlier sample indices */
  std::set<unsigned int>    OutlierIndices;

  int                       ErrorCode;

  double                    MinimumOrientationDifferenceDeg;

  double                    PositionDifferenceThresholdMm;
  double                    OrientationDifferenceThresholdDegrees;

  int                       PoseBucketSize;
  int                       MaximumNumberOfPoseBuckets;
  double                    MaximumPoseBucketError;

  struct MarkerToReferenceTransformMatrixBucket
  {
    MarkerToReferenceTransformMatrixBucket()
    {
      this->MarkerToReferenceCalibrationPoints = std::vector< vtkSmartPointer<vtkMatrix4x4> >();
    }
    std::vector< vtkSmartPointer<vtkMatrix4x4> > MarkerToReferenceCalibrationPoints;
  };
  std::deque<MarkerToReferenceTransformMatrixBucket>  MarkerToReferenceTransformMatrixBuckets;

  class vtkInternal;
  vtkInternal* Internal;
};

#endif
