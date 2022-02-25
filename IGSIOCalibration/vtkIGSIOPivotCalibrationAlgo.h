/*=IGSIO=header=begin======================================================
Program: IGSIO
Copyright (c) Laboratory for Percutaneous Surgery. All rights reserved.
See License.txt for details.
=========================================================IGSIO=header=end*/

#ifndef __vtkIGSIOPivotCalibrationAlgo_h
#define __vtkIGSIOPivotCalibrationAlgo_h

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
  \class vtkIGSIOPivotCalibrationAlgo
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
class vtkIGSIOCalibrationExport vtkIGSIOPivotCalibrationAlgo : public vtkObject
{
public:
  vtkTypeMacro(vtkIGSIOPivotCalibrationAlgo, vtkObject);
  static vtkIGSIOPivotCalibrationAlgo* New();

  enum Events
  {
    InputTransformAddedEvent = vtkCommand::UserEvent + 173,
    AutoCalibrationCompleteEvent
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

  enum AutoCalibrationModeTypes
  {
    PIVOT_CALIBRATION,
    SPIN_CALIBRATION
  };

  /*!
  * Read configuration
  * \param aConfig Root element of the device set configuration
  */
  virtual igsioStatus ReadConfiguration(vtkXMLDataElement* aConfig);

  /*!
    Remove all previously inserted calibration points.
    Call this method to get rid of previously added calibration points
    before starting a new calibration.
  */
  void RemoveAllCalibrationPoints();

  /*!
    Insert acquired point to calibration point list
    \param aMarkerToReferenceTransformMatrix New calibration point (tool to reference transform)
  */
  igsioStatus InsertNextCalibrationPoint(vtkMatrix4x4* aMarkerToReferenceTransformMatrix);

  /*!
    Calibrate (call the minimizer and set the result)
    \param aTransformRepository Transform repository to save the results into
  */
  igsioStatus DoPivotCalibration(vtkIGSIOTransformRepository* aTransformRepository = NULL, bool autoOrient = true);
  igsioStatus DoSpinCalibration(vtkIGSIOTransformRepository* aTransformRepository = NULL, bool snapRotation = false, bool autoOrient = true);

  /*!
    Get calibration result string to display
    \param aPrecision Number of decimals shown
    \return Calibration result (e.g. stylus tip to stylus translation) string
  */
  std::string GetPivotPointToMarkerTranslationString(double aPrecision = 3);

  /*!
    Get the number of outlier points. It is recommended to display a warning to the user
    if the percentage of outliers vs total number of points is larger than a few percent.
  */
  int GetNumberOfDetectedOutliers();

  int GetNumberOfCalibrationPoints();

  // Computes the maximum orientation difference in degrees between the first tool transformation
  // and all the others. Used for determining if there was enough variation in the input data.
  double GetMaximumToolOrientationDifferenceDeg();

  //@{
  /// Mean error of the pivot calibration result in mm
  vtkGetMacro(PivotCalibrationErrorMm, double);
  //@}

  //@{
  /// Mean error of the spin calibration result in mm
  vtkGetMacro(SpinCalibrationErrorMm, double);
  //@}

  //@{
  /// Output calibration tip to marker position
  vtkGetObjectMacro(PivotPointToMarkerTransformMatrix, vtkMatrix4x4);
  vtkGetVector3Macro(PivotPointPosition_Reference, double);
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
  /// Flag that indicates if auto calibration is enabled.
  /// While auto calibration is enabled, the algorithm will maintain a buffer of input points stored in buckets (\sa SetAutoCalibrationBucketSize, \sa SetAutoCalibrationMaximumNumberOfBuckets).
  /// If the error in the latest bucket is too high, then all input data will be discarded (\sa SetAutoCalibrationMaximumBucketError).
  /// When the required number of points are recorded (\sa SetAutoCalibrationNumberOfPoints), and the error is below the threshold (\sa SetAutoCalibrationTargetError), calibration will be run and AutoCalibrationCompleteEvent will be invoked.
  /// \sa SetAutoCalibrationMode
  vtkGetMacro(AutoCalibrationEnabled, bool);
  vtkSetMacro(AutoCalibrationEnabled, bool);
  vtkBooleanMacro(AutoCalibrationEnabled, bool);
  //@}

  //@{
  /// The number of input points to be stored in each bucket
  vtkGetMacro(AutoCalibrationBucketSize, int);
  vtkSetMacro(AutoCalibrationBucketSize, int);
  //@}

  //@{
  /// The maximum number of buckets to be stored. If the number of buckets is exceeded, the oldest will be discarded.
  vtkGetMacro(AutoCalibrationMaximumNumberOfBuckets, int);
  vtkSetMacro(AutoCalibrationMaximumNumberOfBuckets, int);
  //@}

  //@{
  /// Target threshold to reach before AutoCalibrationCompleteEvent will be invoked.
  vtkGetMacro(AutoCalibrationTargetError, double);
  vtkSetMacro(AutoCalibrationTargetError, double);

  //@{
  /// Number of points to record before attempting to perform calibration
  vtkGetMacro(AutoCalibrationNumberOfPoints, int);
  vtkSetMacro(AutoCalibrationNumberOfPoints, int);
  //@}

  //@{
  /// Calibration mode to use for auto calibration.
  /// Value should be either PIVOT_CALIBRATION or SPIN_CALIBRATION. Default is PIVOT_CALIBRATION.
  vtkGetMacro(AutoCalibrationMode, int);
  vtkSetMacro(AutoCalibrationMode, int);
  //@}

  //@{
  /// The accepted amount of error within a single bucket. If the error in the latest bucket exceeds this amount, then all of the buckets will be discarded.
  vtkGetMacro(AutoCalibrationMaximumBucketError, double);
  vtkSetMacro(AutoCalibrationMaximumBucketError, double);
  //@}

  //@{
  /// Automatically flips the shaft direction to be consistent with the needle orientation protocol.
  vtkGetMacro(AutoCalibrationAutoOrient, bool);
  vtkSetMacro(AutoCalibrationAutoOrient, bool);
  //@}

  //@{
  /// Snaps the rotation to be a 90 degree rotation about one of the coordinate axes.
  vtkGetMacro(AutoCalibrationSnapRotation, bool);
  vtkSetMacro(AutoCalibrationSnapRotation, bool);
  //@}

protected:
  vtkSetObjectMacro(PivotPointToMarkerTransformMatrix, vtkMatrix4x4);
  vtkSetStringMacro(ObjectMarkerCoordinateFrame);
  vtkSetStringMacro(ReferenceCoordinateFrame);
  vtkSetStringMacro(ObjectPivotPointCoordinateFrame);
  vtkSetMacro(ErrorCode, int);

protected:
  vtkIGSIOPivotCalibrationAlgo();
  virtual ~vtkIGSIOPivotCalibrationAlgo();

protected:
  /*! Compute the mean position error of the pivot point (in mm) */
  void ComputePivotCalibrationError();

  double ComputePivotCalibrationError(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, std::set<unsigned int>* outlierIndices, double* pivotPoint_Reference, vtkMatrix4x4* pivotPointToMarkerTransformMatrix);

  igsioStatus GetPivotPointPosition(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, std::set<unsigned int>* outlierIndices, double* pivotPoint_Marker, double* pivotPoint_Reference);

  std::vector<vtkMatrix4x4*> GetMarkerToReferenceTransformMatrixArray();
  void GetMarkerToReferenceTransformMatrixArray(std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray);
  std::vector<vtkMatrix4x4*> GetMarkerToReferenceTransformMatrixArray(int bucket);
  void GetMarkerToReferenceTransformMatrixArray(int bucket, std::vector<vtkMatrix4x4*>* matrixArray);

  /// Removes all data if the latest bucket is not within the acceptable error threshold.
  igsioStatus CleanInputBuffer();

  /// Attempt to perform automatic calibration
  /// If enough points have been recorded, and the error is less than the threshold, then AutoCalibrationCompleteEvent will be invoked.
  igsioStatus AutoCalibrate();

  igsioStatus DoPivotCalibrationInternal(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, bool autoOrient, std::set<unsigned int>* outlierIndices, double pivotPoint_Marker[4], double pivotPoint_Reference[4], vtkMatrix4x4* pivotPointToMarkerTransformMatrix);
  igsioStatus DoSpinCalibrationInternal(const std::vector<vtkMatrix4x4*>* markerToTransformMatrixArray, bool snapRotation, bool autoOrient, vtkMatrix4x4* pivotPointToMarkerTransformMatrix, double& error);

  // Verify whether the tool's shaft is in the same direction as the ToolTip to Tool vector.
  // Rotate the ToolTip coordinate frame by 180 degrees about the secondary axis to make the
  // shaft in the same direction as the ToolTip to Tool vector, if this is not already the case.
  void UpdateShaftDirection(vtkMatrix4x4* toolTipToToolMatrix);

  // Flip the direction of the shaft axis
  void FlipShaftDirection(vtkMatrix4x4* toolTipToToolMatrix);

  void GetToolTipToToolRotation(vtkMatrix4x4* toolTipToToolMatrix, vtkMatrix4x4* rotationMatrix);

  // Returns the orientation difference in degrees between two 4x4 homogeneous transformation matrix, in degrees.
  double GetOrientationDifferenceDeg(vtkMatrix4x4* aMatrix, vtkMatrix4x4* bMatrix);

  // Helper method to compute the secondary axis, given a shaft axis
  vnl_vector< double > vtkIGSIOPivotCalibrationAlgo::ComputeSecondaryAxis(vnl_vector< double > shaftAxis_ToolTip);

protected:
  vtkMatrix4x4* PivotPointToMarkerTransformMatrix;
  vtkMatrix4x4* PreviousMarkerToReferenceTransformMatrix;

  double                    PivotCalibrationErrorMm;
  double                    SpinCalibrationErrorMm;


  char*                     ObjectMarkerCoordinateFrame;
  char*                     ReferenceCoordinateFrame;
  char*                     ObjectPivotPointCoordinateFrame;

  /*! Pivot point position in the Reference coordinate system */
  double                    PivotPointPosition_Reference[4];

  /*! List of outlier sample indices */
  std::set<unsigned int>    OutlierIndices;

  int                       ErrorCode;

  double                    MinimumOrientationDifferenceDeg;

  double                    PositionDifferenceThresholdMm;
  double                    OrientationDifferenceThresholdDegrees;

  bool                      AutoCalibrationEnabled;
  int                       AutoCalibrationBucketSize;
  int                       AutoCalibrationNumberOfPoints;
  double                    AutoCalibrationTargetError;
  double                    AutoCalibrationMaximumBucketError;
  int                       AutoCalibrationMaximumNumberOfBuckets;
  int                       AutoCalibrationMode;
  bool                      AutoCalibrationAutoOrient;
  bool                      AutoCalibrationSnapRotation;

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
