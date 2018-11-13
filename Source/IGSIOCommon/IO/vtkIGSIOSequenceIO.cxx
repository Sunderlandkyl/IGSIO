/*=Plus=header=begin======================================================
  Program: Plus
  Copyright (c) Laboratory for Percutaneous Surgery. All rights reserved.
  See License.txt for details.
=========================================================Plus=header=end*/

//#include "PlusConfigure.h"

#include "vtkIGSIOMetaImageSequenceIO.h"
#include "vtkIGSIONrrdSequenceIO.h"
#include "vtkIGSIOSequenceIO.h"
#include "vtkIGSIOTrackedFrameList.h"

/// VTK includes
#include <vtkNew.h>

#ifdef PLUS_USE_VTKVIDEOIO_MKV
  #include "vtkPlusMkvSequenceIO.h"
#endif

//----------------------------------------------------------------------------
igsioStatus vtkIGSIOSequenceIO::Write(const std::string& filename, vtkIGSIOTrackedFrameList* frameList, US_IMAGE_ORIENTATION orientationInFile/*=US_IMG_ORIENT_MF*/, bool useCompression/*=true*/, bool enableImageDataWrite/*=true*/)
{
  // Convert local filename to plus output filename
  if (vtksys::SystemTools::FileExists(filename.c_str()))
  {
    // Remove the file before replacing it
    vtksys::SystemTools::RemoveFile(filename.c_str());
  }

   // Parse sequence filename to determine if it's metafile or NRRD
   if (vtkIGSIOMetaImageSequenceIO::CanWriteFile(filename))
   {
     if (frameList->SaveToSequenceMetafile(filename, orientationInFile, useCompression, enableImageDataWrite) != IGSIO_SUCCESS)
     {
       //**LOG_ERROR("Unable to save file: " << filename << " as sequence metafile.");
       return IGSIO_FAIL;
     }
     return IGSIO_SUCCESS;
   }
  else if (vtkIGSIONrrdSequenceIO::CanWriteFile(filename))
  {
    if (frameList->SaveToNrrdFile(filename, orientationInFile, useCompression, enableImageDataWrite) != IGSIO_SUCCESS)
    {
      //**LOG_ERROR("Unable to save file: " << filename << " as Nrrd file.");
      return IGSIO_FAIL;
    }

    return IGSIO_SUCCESS;
  }
#ifdef PLUS_USE_VTKVIDEOIO_MKV
  else if (vtkPlusMkvSequenceIO::CanWriteFile(filename))
  {
    if (frameList->SaveToMatroskaFile(filename, orientationInFile, useCompression, enableImageDataWrite) != IGSIO_SUCCESS)
    {
      //**LOG_ERROR("Unable to save file: " << filename << " as MKV file.");
      return IGSIO_FAIL;
    }
  }
#endif

  //**LOG_ERROR("No writer for file: " << filename);
  return IGSIO_FAIL;
}

//----------------------------------------------------------------------------
igsioStatus vtkIGSIOSequenceIO::Write(const std::string& filename, igsioTrackedFrame* frame, US_IMAGE_ORIENTATION orientationInFile /*= US_IMG_ORIENT_MF*/, bool useCompression /*= true*/, bool EnableImageDataWrite /*= true*/)
{
  vtkNew<vtkIGSIOTrackedFrameList> list;
  list->AddTrackedFrame(frame);
  return vtkIGSIOSequenceIO::Write(filename, list.GetPointer(), orientationInFile, useCompression, EnableImageDataWrite);
}

//----------------------------------------------------------------------------
igsioStatus vtkIGSIOSequenceIO::Read(const std::string& filename, vtkIGSIOTrackedFrameList* frameList)
{
  if (!vtksys::SystemTools::FileExists(filename.c_str()))
  {
    //**LOG_ERROR("File: " << filename << " does not exist.");
    return IGSIO_FAIL;
  }

  if (vtkIGSIOMetaImageSequenceIO::CanReadFile(filename))
  {
    // Attempt metafile read
    if (frameList->ReadFromSequenceMetafile(filename) != IGSIO_SUCCESS)
    {
      //**LOG_ERROR("Failed to read video buffer from sequence metafile: " << filename);
      return IGSIO_FAIL;
    }
    return IGSIO_SUCCESS;
  }
  // Parse sequence filename to determine if it's metafile or NRRD
  else if (vtkIGSIONrrdSequenceIO::CanReadFile(filename))
  {
    // Attempt Nrrd read
    if (frameList->ReadFromNrrdFile(filename.c_str()) != IGSIO_SUCCESS)
    {
      //**LOG_ERROR("Failed to read video buffer from Nrrd file: " << filename);
      return IGSIO_FAIL;
    }

    return IGSIO_SUCCESS;
  }
#ifdef PLUS_USE_VTKVIDEOIO_MKV
  else if (vtkPlusMkvSequenceIO::CanReadFile(filename))
  {
    // Attempt MKV read
    if (frameList->ReadFromMatroskaFile(filename.c_str()) != IGSIO_SUCCESS)
    {
      //**LOG_ERROR("Failed to read video buffer from MKV file: " << filename);
      return IGSIO_FAIL;
    }

    return IGSIO_SUCCESS;
  }
#endif

  //**LOG_ERROR("No reader for file: " << filename);
  return IGSIO_FAIL;
}

//----------------------------------------------------------------------------
vtkIGSIOSequenceIOBase* vtkIGSIOSequenceIO::CreateSequenceHandlerForFile(const std::string& filename)
{
  // Parse sequence filename to determine if it's metafile or NRRD
  if (vtkIGSIOMetaImageSequenceIO::CanWriteFile(filename))
  {
    return vtkIGSIOMetaImageSequenceIO::New();
  }
  else if (vtkIGSIONrrdSequenceIO::CanWriteFile(filename))
  {
    return vtkIGSIONrrdSequenceIO::New();
  }
#ifdef PLUS_USE_VTKVIDEOIO_MKV
  else if (vtkPlusMkvSequenceIO::CanReadFile(filename))
  {
    return vtkPlusMkvSequenceIO::New();
  }
#endif

  //**LOG_ERROR("No writer for file: " << filename);
  return NULL;
}
