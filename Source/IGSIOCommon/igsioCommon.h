/*=Plus=header=begin======================================================
  Program: Plus
  Copyright (c) Laboratory for Percutaneous Surgery. All rights reserved.
  See License.txt for details.
=========================================================Plus=header=end*/

#ifndef __IGSIOCommon_h
#define __IGSIOCommon_h

// IGSIO includes
#include "vtkigsiocommon_export.h"
//#include "itkImageIOBase.h"
//#include "vtkPlusLogger.h"
//#include "vtkPlusMacro.h"
//#include "vtkPlusRecursiveCriticalSection.h"

// VTK includes
#include <vtkImageData.h>
#include <vtksys/SystemTools.hxx>

// System includes
#include <float.h>

// STL includes
#include <array>
#include <list>
#include <locale>
#include <sstream>

//class vtkPlusUsScanConvert;
//class vtkPlusTrackedFrameList;

//const int MAX_PATH = 1000;

enum igsioStatus
{
  IGSIO_FAIL = 0,
  IGSIO_SUCCESS = 1
};

enum igsioImagingMode
{
  IGSIO_UnknownMode,
  IGSIO_BMode,
  IGSIO_RfMode
};

typedef std::array<unsigned int, 3> FrameSizeType;

#define UNDEFINED_TIMESTAMP DBL_MAX

/* Define case insensitive string compare for Windows. */
#if defined( _WIN32 ) && !defined(__CYGWIN__)
#if defined(__BORLANDC__)
#define STRCASECMP stricmp
#else
#define STRCASECMP _stricmp
#endif
#else
#define STRCASECMP strcasecmp
#endif

///////////////////////////////////////////////////////////////////

/*!
\class igsioLockGuard
\brief A class for automatically unlocking objects

This class is used for locking a objects (buffers, mutexes, etc.)
and releasing the lock automatically when the guard object is deleted
(typically by getting out of scope).

Example:
\code
igsioLockGuard<vtkPlusRecursiveCriticalSection> updateMutexGuardedLock(this->UpdateMutex);
\endcode

\ingroup PlusLibCommon
*//*
template <typename T>
class igsioLockGuard
{
public:
  igsioLockGuard(T* lockableObject)
  {
    m_LockableObject = lockableObject;
    m_LockableObject->Lock();
  }
  ~igsioLockGuard()
  {
    m_LockableObject->Unlock();
    m_LockableObject = NULL;
  }
private:
  igsioLockGuard(igsioLockGuard&);
  void operator=(igsioLockGuard&);

  T* m_LockableObject;
};*/

/*!
\def DELETE_IF_NOT_NULL(Object)
\brief A macro to safely delete a VTK object (usable if the VTK object pointer is already NULL).
\ingroup PlusLibCommon
*/
#define DELETE_IF_NOT_NULL( Object ) {\
  if ( Object != NULL ) {\
    Object->Delete();\
    Object = NULL;\
  }\
}

//
// Get character string.  Creates member Get"name"()
// (e.g., char *GetFilename());
//
#define vtkGetStdStringMacro(name) \
virtual std::string Get##name () const { \
  return this->name; \
}

#define GetStdStringMacro(name) \
virtual std::string Get##name () const { \
  return this->name; \
}

//
// Set character string.  Creates member Set"name"()
// (e.g., SetFilename(char *));
//
#define vtkSetStdStringMacro(name) \
virtual void Set##name (const std::string& _arg) \
{ \
  if ( this->name.compare(_arg) == 0 ) { return;} \
  this->name = _arg; \
  this->Modified(); \
} \
virtual void Set##name (const char* _arg) \
{ \
  this->Set##name(std::string(_arg ? _arg : "")); \
}

#define SetStdStringMacro(name) \
virtual void Set##name (const std::string& _arg) \
{ \
  if ( this->name.compare(_arg) == 0 ) { return;} \
  this->name = _arg; \
} \
virtual void Set##name (const char* _arg) \
{ \
  this->Set##name(std::string(_arg ? _arg : "")); \
}

#define vtkGetMacroConst(name,type) \
virtual type Get##name () const { \
  return this->name; \
}

#define SetMacro(name,type) \
virtual void Set##name (type _arg) \
{ \
  if (this->name != _arg) \
  { \
    this->name = _arg; \
  } \
}

#define GetMacro(name,type) \
virtual void Set##name (type _arg) \
{ \
  if (this->name != _arg) \
  { \
    this->name = _arg; \
  } \
}

class vtkPlusTrackedFrameList;
class vtkXMLDataElement;

namespace igsioCommon
{
  //typedef itk::ImageIOBase::IOComponentType ITKScalarPixelType;
  typedef int VTKScalarPixelType;
  typedef int IGTLScalarPixelType;

  //----------------------------------------------------------------------------
  /*! Quick and robust string to int conversion */
  template<class T>
  igsioStatus StringToInt(const char* strPtr, T& result)
  {
    if (strPtr == NULL || strlen(strPtr) == 0)
    {
      return IGSIO_FAIL;
    }
    char* pEnd = NULL;
    result = static_cast<int>(strtol(strPtr, &pEnd, 10));
    if (pEnd != strPtr + strlen(strPtr))
    {
      return IGSIO_FAIL;
    }
    return IGSIO_SUCCESS;
  }

  //----------------------------------------------------------------------------
  /*! Quick and robust string to int conversion */
  template<class T>
  igsioStatus StringToUInt(const char* strPtr, T& result)
  {
    if (strPtr == NULL || strlen(strPtr) == 0)
    {
      return IGSIO_FAIL;
    }
    char* pEnd = NULL;
    result = static_cast<unsigned int>(strtol(strPtr, &pEnd, 10));
    if (pEnd != strPtr + strlen(strPtr))
    {
      return IGSIO_FAIL;
    }
    return IGSIO_SUCCESS;
  }

  //----------------------------------------------------------------------------
  /*! Quick and robust string to double conversion */
  template<class T>
  igsioStatus StringToDouble(const char* strPtr, T& result)
  {
    if (strPtr == NULL || strlen(strPtr) == 0)
    {
      return IGSIO_FAIL;
    }
    char* pEnd = NULL;
    result = strtod(strPtr, &pEnd);
    if (pEnd != strPtr + strlen(strPtr))
    {
      return IGSIO_FAIL;
    }
    return IGSIO_SUCCESS;
  }
  //---------------------------------------------------------------------------
  struct ImageMetaDataItem
  {
    std::string Id;  /* device name to query the IMAGE and COLORT */
    std::string Description;        /* name / description (< 64 bytes)*/
    std::string Modality;    /* modality name (< 32 bytes) */
    std::string PatientName; /* patient name (< 64 bytes) */
    std::string PatientId;   /* patient ID (MRN etc.) (< 64 bytes) */
    double TimeStampUtc;   /* scan time in UTC*/
    unsigned int Size[3];     /* entire image volume size */
    unsigned char ScalarType;  /* scalar type. see scalar_type in IMAGE message */
  };

  typedef std::list<ImageMetaDataItem> ImageMetaDataList;
  //----------------------------------------------------------------------------
  /*! Quick and robust string to int conversion */
  template<class T>
  igsioStatus StringToLong(const char* strPtr, T& result)
  {
    if (strPtr == NULL || strlen(strPtr) == 0)
    {
      return IGSIO_FAIL;
    }
    char* pEnd = NULL;
    result = strtol(strPtr, &pEnd, 10);
    if (pEnd != strPtr + strlen(strPtr))
    {
      return IGSIO_FAIL;
    }
    return IGSIO_SUCCESS;
  }

  //----------------------------------------------------------------------------
  enum LINE_STYLE
  {
    LINE_STYLE_SOLID,
    LINE_STYLE_DOTS,
  };

  enum ALPHA_BEHAVIOR
  {
    ALPHA_BEHAVIOR_SOURCE,
    ALPHA_BEHAVIOR_OPAQUE
  };

  //----------------------------------------------------------------------------
  //VTKIGSIOCOMMON_EXPORT igsioStatus DrawLine(vtkImageData& imageData,
  //  const std::array<float, 3>& colour,
  //  LINE_STYLE style,
  //  const std::array<int, 3>& startPixel,
  //  const std::array<int, 3>& endPixel,
  //  unsigned int numberOfPoints,
  //  ALPHA_BEHAVIOR alphaBehavior = ALPHA_BEHAVIOR_OPAQUE);

  //----------------------------------------------------------------------------
  //VTKIGSIOCOMMON_EXPORT igsioStatus DrawLine(vtkImageData& imageData,
  //  float greyValue,
  //  LINE_STYLE style,
  //  const std::array<int, 3>& startPixel,
  //  const std::array<int, 3>& endPixel,
  //  unsigned int numberOfPoints,
  //  ALPHA_BEHAVIOR alphaBehavior = ALPHA_BEHAVIOR_OPAQUE);

  //----------------------------------------------------------------------------
  bool VTKIGSIOCOMMON_EXPORT IsEqualInsensitive(std::string const& a, std::string const& b);
  bool VTKIGSIOCOMMON_EXPORT IsEqualInsensitive(std::wstring const& a, std::wstring const& b);
  bool VTKIGSIOCOMMON_EXPORT HasSubstrInsensitive(std::string const& a, std::string const& b);
  bool VTKIGSIOCOMMON_EXPORT HasSubstrInsensitive(std::wstring const& a, std::wstring const& b);

  //----------------------------------------------------------------------------
  typedef std::array<int, 3> PixelPoint;
  typedef std::pair<PixelPoint, PixelPoint> PixelLine;
  typedef std::vector<PixelLine> PixelLineList;
  //VTKIGSIOCOMMON_EXPORT igsioStatus DrawScanLines(int* inputImageExtent, float greyValue, const PixelLineList& scanLineEndPoints, vtkPlusTrackedFrameList* trackedFrameList);
  //VTKIGSIOCOMMON_EXPORT igsioStatus DrawScanLines(int* inputImageExtent, const std::array<float, 3>& colour, const PixelLineList& scanLineEndPoints, vtkPlusTrackedFrameList* trackedFrameList);
  //VTKIGSIOCOMMON_EXPORT igsioStatus DrawScanLines(int* inputImageExtent, float greyValue, const PixelLineList& scanLineEndPoints, vtkImageData* imageData);
  //VTKIGSIOCOMMON_EXPORT igsioStatus DrawScanLines(int* inputImageExtent, const std::array<float, 3>& colour, const PixelLineList& scanLineEndPoints, vtkImageData* imageData);

  //----------------------------------------------------------------------------
  template<typename T>
  static std::string ToString(T number)
  {
#if defined(_MSC_VER) && _MSC_VER < 1700
    // This method can be used for number to string conversion
    // until std::to_string is supported by more compilers.
    std::ostringstream ss;
    ss << number;
    return ss.str();
#else
    return std::to_string(number);
#endif
  }

  static const int NO_CLIP = -1;
  VTKIGSIOCOMMON_EXPORT bool IsClippingRequested(const std::array<int, 3>& clipOrigin, const std::array<int, 3>& clipSize);
  VTKIGSIOCOMMON_EXPORT bool IsClippingWithinExtents(const std::array<int, 3>& clipOrigin, const std::array<int, 3>& clipSize, const int extents[6]);

  VTKIGSIOCOMMON_EXPORT void SplitStringIntoTokens(const std::string& s, char delim, std::vector<std::string>& elems, bool keepEmptyParts = true);
  VTKIGSIOCOMMON_EXPORT std::vector<std::string> SplitStringIntoTokens(const std::string& s, char delim, bool keepEmptyParts = true);
  template<typename ElemType>
  VTKIGSIOCOMMON_EXPORT void JoinTokensIntoString(const std::vector<ElemType>& elems, std::string& output, char separator = ' ');

  VTKIGSIOCOMMON_EXPORT igsioStatus CreateTemporaryFilename(std::string& aString, const std::string& anOutputDirectory);

  /*! Trim whitespace characters from the left and right */
  VTKIGSIOCOMMON_EXPORT std::string& Trim(std::string& str);

  /*!
  On some systems fwrite may fail if a large chunk of data is attempted to written in one piece.
  This method writes the data in smaller chunks as long as all data is written or no data
  can be written anymore.
  */
  VTKIGSIOCOMMON_EXPORT igsioStatus RobustFwrite(FILE* fileHandle, void* data, size_t dataSize, size_t& writtenSize);

  VTKIGSIOCOMMON_EXPORT std::string GetPlusLibVersionString();

  //----------------------------------------------------------------------------
  namespace XML
  {
    /*!
    Writes an XML element to file. The output is nicer that with the built-in vtkXMLDataElement::PrintXML, as
    there are no extra lines, if there are many attributes then each of them is printed on separate line, and
    matrix elements (those that contain Matrix or Transform in the attribute name and 16 numerical elements in the attribute value)
    are printed in 4 lines.
    */
    VTKIGSIOCOMMON_EXPORT igsioStatus PrintXML(const std::string& filename, vtkXMLDataElement* elem);
    /*!
    Writes an XML element to a stream. The output is nicer that with the built-in vtkXMLDataElement::PrintXML, as
    there are no extra lines, if there are many attributes then each of them is printed on separate line, and
    matrix elements (those that contain Matrix or Transform in the attribute name and 16 numerical elements in the attribute value)
    are printed in 4 lines.
    */
    VTKIGSIOCOMMON_EXPORT igsioStatus PrintXML(ostream& os, vtkIndent indent, vtkXMLDataElement* elem);

    VTKIGSIOCOMMON_EXPORT igsioStatus SafeCheckAttributeValueInsensitive(vtkXMLDataElement& element, const std::string& attributeName, const std::string& value, bool& isEqual);
    VTKIGSIOCOMMON_EXPORT igsioStatus SafeCheckAttributeValueInsensitive(vtkXMLDataElement& element, const std::wstring& attributeName, const std::wstring& value, bool& isEqual);

    VTKIGSIOCOMMON_EXPORT igsioStatus SafeGetAttributeValueInsensitive(vtkXMLDataElement& element, const std::string& attributeName, std::string& value);
    VTKIGSIOCOMMON_EXPORT igsioStatus SafeGetAttributeValueInsensitive(vtkXMLDataElement& element, const std::wstring& attributeName, std::string& value);
    template<typename T> VTKIGSIOCOMMON_EXPORT igsioStatus SafeGetAttributeValueInsensitive(vtkXMLDataElement& element, const std::string& attributeName, T& value);
    template<typename T> VTKIGSIOCOMMON_EXPORT igsioStatus SafeGetAttributeValueInsensitive(vtkXMLDataElement& element, const std::wstring& attributeName, T& value);
  }

};

/*!
\class igsioTransformName
\brief Stores the from and to coordinate frame names for transforms

The igsioTransformName stores and generates the from and to coordinate frame names for transforms.
To enable robust serialization to/from a simple string (...To...Transform), the coordinate frame names must
start with an uppercase character and it shall not contain "To" followed by an uppercase character. E.g., valid
coordinate frame names are Tracker, TrackerBase, Tool; invalid names are tracker, trackerBase, ToImage.

Example usage:
Setting a transform name:
\code
igsioTransformName tnImageToProbe("Image", "Probe");
\endcode
or
\code
igsioTransformName tnImageToProbe;
if ( tnImageToProbe->SetTransformName("ImageToProbe") != IGSIO_SUCCESS )
{
LOG_ERROR("Failed to set transform name!");
return IGSIO_FAIL;
}
\endcode
Getting coordinate frame or transform names:
\code
std::string fromFrame = tnImageToProbe->From();
std::string toFrame = tnImageToProbe->To();
\endcode
or
\code
std::string strImageToProbe;
if ( tnImageToProbe->GetTransformName(strImageToProbe) != IGSIO_SUCCESS )
{
LOG_ERROR("Failed to get transform name!");
return IGSIO_FAIL;
}
\endcode

\ingroup PlusLibCommon
*/
class VTKIGSIOCOMMON_EXPORT igsioTransformName
{
public:
  igsioTransformName();
  ~igsioTransformName();
  igsioTransformName(std::string aFrom, std::string aTo);
  igsioTransformName(const std::string& transformName);

  /*!
  Set 'From' and 'To' coordinate frame names from a combined transform name with the following format [FrameFrom]To[FrameTo].
  The combined transform name might contain only one 'To' phrase followed by a capital letter (e.g. ImageToToProbe is not allowed)
  and the coordinate frame names should be in camel case format starting with capitalized letters.
  */
  igsioStatus SetTransformName(const char* aTransformName);
  igsioStatus SetTransformName(const std::string& aTransformName);

  /*! Return combined transform name between 'From' and 'To' coordinate frames: [From]To[To] */
  igsioStatus GetTransformName(std::string& aTransformName) const;
  std::string GetTransformName() const;

  /*! Return 'From' coordinate frame name, give a warning if it's not capitalized and capitalize it*/
  std::string From() const;

  /*! Return 'To' coordinate frame name, give a warning if it's not capitalized and capitalize it */
  std::string To() const;

  /*! Clear the 'From' and 'To' fields */
  void Clear();

  /*! Check if the current transform name is valid */
  bool IsValid() const;

  inline bool operator== (const igsioTransformName& in) const
  {
    return (in.m_From == m_From && in.m_To == m_To);
  }

  inline bool operator!= (const igsioTransformName& in) const
  {
    return !(in == *this);
  }

  friend std::ostream& operator<< (std::ostream& os, const igsioTransformName& transformName)
  {
    os << transformName.GetTransformName();
    return os;
  }

private:
  /*! Check if the input string is capitalized, if not capitalize it */
  void Capitalize(std::string& aString);
  std::string m_From; /*! From coordinate frame name */
  std::string m_To; /*! To coordinate frame name */
};


#define RETRY_UNTIL_TRUE(command_, numberOfRetryAttempts_, delayBetweenRetryAttemptsSec_) \
  { \
    bool success = false; \
    int numOfTries = 0; \
    while ( !success && numOfTries < numberOfRetryAttempts_ ) \
    { \
      success = (command_);   \
      if (success)  \
      { \
        /* command successfully completed, continue without waiting */ \
        break; \
      } \
      /* command failed, wait for some time and retry */ \
      numOfTries++;   \
      vtkPlusAccurateTimer::Delay(delayBetweenRetryAttemptsSec_); \
    } \
  }

//#include "vtkPlusConfig.h"
//#include "PlusXmlUtils.h"
#include "igsioCommon.txx"

#endif //__igsioCommon_h
