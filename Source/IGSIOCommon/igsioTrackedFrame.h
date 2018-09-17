/*=Plus=header=begin======================================================
  Program: Plus
  Copyright (c) Laboratory for Percutaneous Surgery. All rights reserved.
  See License.txt for details.
=========================================================Plus=header=end*/

#ifndef __IGSIOTRACKEDFRAME_H
#define __IGSIOTRACKEDFRAME_H

#include "vtkigsiocommon_export.h"

#include "IGSIOCommon.h"
#include "igsioVideoFrame.h"

class vtkMatrix4x4;
class vtkPoints;

/*!
  \enum igsioTrackedFrameFieldStatus
  \brief Tracked frame field status
  Image field is valid if the image data is not NULL.
  Tool status is valid only if the ToolStatus is TOOL_OK.
  \ingroup PlusLibCommon
*/
enum igsioTrackedFrameFieldStatus
{
  FIELD_OK,            /*!< Field is valid */
  FIELD_INVALID       /*!< Field is invalid */
};

/*!
  \class igsioTrackedFrame
  \brief Stores tracked frame (image + pose information)
  \ingroup PlusLibCommon
*/
class VTKIGSIOCOMMON_EXPORT igsioTrackedFrame
{
public:
  static const char* FIELD_FRIENDLY_DEVICE_NAME;

public:
  static const std::string TransformPostfix;
  static const std::string TransformStatusPostfix;
  typedef std::map<std::string, std::string> FieldMapType;

public:
  igsioTrackedFrame();
  ~igsioTrackedFrame();
  igsioTrackedFrame(const igsioTrackedFrame& frame);
  igsioTrackedFrame& operator=(igsioTrackedFrame const& trackedFrame);

public:
  /*! Set image data */
  void SetImageData(const igsioVideoFrame& value);

  /*! Get image data */
  igsioVideoFrame* GetImageData() { return &(this->ImageData); };

  /*! Set timestamp */
  void SetTimestamp(double value);

  /*! Get timestamp */
  double GetTimestamp() { return this->Timestamp; };

  /*! Set custom frame field */
  void SetCustomFrameField(std::string name, std::string value);

  /*! Get custom frame field value */
  const char* GetCustomFrameField(const char* fieldName);
  const char* GetCustomFrameField(const std::string& fieldName);

  /*! Delete custom frame field */
  igsioStatus DeleteCustomFrameField(const char* fieldName);

  /*!
    Check if a custom frame field is defined or not
    \return true, if the field is defined; false, if the field is not defined
  */
  bool IsCustomFrameFieldDefined(const char* fieldName);

  /*!
    Check if a custom frame transform name field is defined or not
    \return true, if the field is defined; false, if the field is not defined
  */
  bool IsCustomFrameTransformNameDefined(const igsioTransformName& transformName);

  /*! Get custom frame transform */
  igsioStatus GetCustomFrameTransform(const igsioTransformName& frameTransformName, double transform[16]);
  /*! Get custom frame transform */
  igsioStatus GetCustomFrameTransform(const igsioTransformName& frameTransformName, vtkMatrix4x4* transformMatrix);

  /*! Get custom frame status */
  igsioStatus GetCustomFrameTransformStatus(const igsioTransformName& frameTransformName, igsioTrackedFrameFieldStatus& status);
  /*! Set custom frame status */
  igsioStatus SetCustomFrameTransformStatus(const igsioTransformName& frameTransformName, igsioTrackedFrameFieldStatus status);

  /*! Set custom frame transform */
  igsioStatus SetCustomFrameTransform(const igsioTransformName& frameTransformName, double transform[16]);

  /*! Set custom frame transform */
  igsioStatus SetCustomFrameTransform(const igsioTransformName& frameTransformName, vtkMatrix4x4* transform);

  /*! Get the list of the name of all custom frame fields */
  void GetCustomFrameFieldNameList(std::vector<std::string>& fieldNames);

  /*! Get the list of the transform name of all custom frame transforms*/
  void GetCustomFrameTransformNameList(std::vector<igsioTransformName>& transformNames);

  /*! Get tracked frame size in pixel. Returns: FrameSizeType.  */
  FrameSizeType GetFrameSize();

  /*! Set tracked frame size in pixel.*/
  void SetFrameSize(FrameSizeType frameSize);

  /*! Get tracked frame pixel size in bits (scalar size * number of scalar components) */
  int GetNumberOfBitsPerScalar();

  /*! Get number of scalar components in a pixel */
  igsioStatus GetNumberOfScalarComponents(unsigned int& numberOfScalarComponents);

  /*! Get number of bits in a pixel */
  int GetNumberOfBitsPerPixel();

  /*! Set Segmented fiducial point pixel coordinates */
  void SetFiducialPointsCoordinatePx(vtkPoints* fiducialPoints);

  /*! Get Segmented fiducial point pixel coordinates */
  vtkPoints* GetFiducialPointsCoordinatePx() { return this->FiducialPointsCoordinatePx; };

  /*! Write image with image to tracker transform to file */
  igsioStatus WriteToFile(const std::string& filename, vtkMatrix4x4* imageToTracker);

  /*! Print tracked frame human readable serialization data to XML data
      If requestedTransforms is empty, all stored CustomFrameFields are sent
  */
  igsioStatus PrintToXML(vtkXMLDataElement* xmlData, const std::vector<igsioTransformName>& requestedTransforms);

  /*! Serialize Tracked frame human readable data to xml data and return in string */
  igsioStatus GetigsioTrackedFrameInXmlData(std::string& strXmlData, const std::vector<igsioTransformName>& requestedTransforms);

  /*! Deserialize igsioTrackedFrame human readable data from xml data string */
  igsioStatus SetigsioTrackedFrameFromXmlData(const char* strXmlData);
  /*! Deserialize igsioTrackedFrame human readable data from xml data string */
  igsioStatus SetigsioTrackedFrameFromXmlData(const std::string& xmlData);

  /*! Convert from field status string to field status enum */
  static igsioTrackedFrameFieldStatus ConvertFieldStatusFromString(const char* statusStr);

  /*! Convert from field status enum to field status string */
  static std::string ConvertFieldStatusToString(igsioTrackedFrameFieldStatus status);

  /*! Return all custom fields in a map */
  const FieldMapType& GetCustomFields() { return this->CustomFrameFields; }

  /*! Returns true if the input string ends with "Transform", else false */
  static bool IsTransform(std::string str);

  /*! Returns true if the input string ends with "TransformStatus", else false */
  static bool IsTransformStatus(std::string str);

public:
  bool operator< (igsioTrackedFrame data) { return Timestamp < data.Timestamp; }
  bool operator== (const igsioTrackedFrame& data) const
  {
    return (Timestamp == data.Timestamp);
  }

protected:
  igsioVideoFrame ImageData;
  double Timestamp;

  FieldMapType CustomFrameFields;

  FrameSizeType FrameSize;

  /*! Stores segmented fiducial point pixel coordinates */
  vtkPoints* FiducialPointsCoordinatePx;
};

//----------------------------------------------------------------------------

/*!
  \enum igsioTrackedFrameValidationRequirements
  \brief If any of the requested requirement is not fulfilled then the validation fails.
  \ingroup PlusLibCommon
*/
enum igsioTrackedFrameValidationRequirements
{
  REQUIRE_UNIQUE_TIMESTAMP = 0x0001, /*!< the timestamp shall be unique */
  REQUIRE_TRACKING_OK = 0x0002, /*!<  the tracking flags shall be valid (TOOL_OK) */
  REQUIRE_CHANGED_ENCODER_POSITION = 0x0004, /*!<  the stepper encoder position shall be different from the previous ones  */
  REQUIRE_SPEED_BELOW_THRESHOLD = 0x0008, /*!<  the frame acquisition speed shall be less than a threshold */
  REQUIRE_CHANGED_TRANSFORM = 0x0010, /*!<  the transform defined by name shall be different from the previous ones  */
};

/*!
  \class igsioTrackedFrameTimestampFinder
  \brief Helper class used for validating timestamps in a tracked frame list
  \ingroup PlusLibCommon
*/
class igsioTrackedFrameTimestampFinder
{
public:
  igsioTrackedFrameTimestampFinder(igsioTrackedFrame* frame): migsioTrackedFrame(frame) {};
  bool operator()(igsioTrackedFrame* newFrame)
  {
    return newFrame->GetTimestamp() == migsioTrackedFrame->GetTimestamp();
  }
  igsioTrackedFrame* migsioTrackedFrame;
};

//----------------------------------------------------------------------------

/*!
  \class igsioTrackedFrameEncoderPositionFinder
  \brief Helper class used for validating encoder position in a tracked frame list
  \ingroup PlusLibCommon
*/
class VTKIGSIOCOMMON_EXPORT igsioTrackedFrameEncoderPositionFinder
{
public:
  igsioTrackedFrameEncoderPositionFinder(igsioTrackedFrame* frame, double minRequiredTranslationDifferenceMm, double minRequiredAngleDifferenceDeg);
  ~igsioTrackedFrameEncoderPositionFinder();

  /*! Get stepper encoder values from the tracked frame */
  static igsioStatus GetStepperEncoderValues(igsioTrackedFrame* trackedFrame, double& probePosition, double& probeRotation, double& templatePosition);

  /*!
    Predicate unary function for std::find_if to validate encoder position
    \return Returning true if the encoder position difference is less than required
  */
  bool operator()(igsioTrackedFrame* newFrame);

protected:
  igsioTrackedFrame* mTrackedFrame;
  double mMinRequiredTranslationDifferenceMm;
  double mMinRequiredAngleDifferenceDeg;
};

//----------------------------------------------------------------------------

/*!
  \class TrackedFrameTransformFinder
  \brief Helper class used for validating frame transform in a tracked frame list
  \ingroup PlusLibCommon
*/
class TrackedFrameTransformFinder
{
public:
  TrackedFrameTransformFinder(igsioTrackedFrame* frame, const igsioTransformName& frameTransformName, double minRequiredTranslationDifferenceMm, double minRequiredAngleDifferenceDeg);
  ~TrackedFrameTransformFinder();

  /*!
    Predicate unary function for std::find_if to validate transform
    \return Returning true if the transform difference is less than required
  */
  bool operator()(igsioTrackedFrame* newFrame);

protected:
  igsioTrackedFrame* mTrackedFrame;
  double mMinRequiredTranslationDifferenceMm;
  double mMinRequiredAngleDifferenceDeg;
  igsioTransformName mFrameTransformName;
};

#endif
