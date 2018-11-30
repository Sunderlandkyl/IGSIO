/*=Plus=header=begin======================================================
Program: Plus
Copyright (c) Laboratory for Percutaneous Surgery. All rights reserved.
See License.txt for details.
=========================================================Plus=header=end*/

// std includes
#include <iomanip>

// VTK includes
#include <vtkMatrix4x4.h>
#include <vtkNew.h>
#include <vtksys/CommandLineArguments.hxx>

// IGSIO includes
#include "igsioTrackedFrame.h"
#include "vtkIGSIOMkvSequenceIO.h"
#include "vtkIGSIOTrackedFrameList.h"

///////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  std::string inputImageSequenceFileName;
  std::string outputImageSequenceFileName;

  int numberOfFailures(0);

  vtksys::CommandLineArguments args;
  args.Initialize(argc, argv);

  args.AddArgument("--input-filename", vtksys::CommandLineArguments::EQUAL_ARGUMENT, &inputImageSequenceFileName, "Filename of the input image sequence.");
  args.AddArgument("--output-filename", vtksys::CommandLineArguments::EQUAL_ARGUMENT, &outputImageSequenceFileName, "Filename of the output image sequence.");

  if (!args.Parse())
  {
    std::cerr << "Problem parsing arguments" << std::endl;
    std::cout << "Help: " << args.GetHelp() << std::endl;
    exit(EXIT_FAILURE);
  }

  if (inputImageSequenceFileName.empty() && outputImageSequenceFileName.empty())
  {
    std::cerr << "--input-filename or --output-filename is required" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (!inputImageSequenceFileName.empty())
  {
    // ******************************************************************************
    // Test reading

    vtkNew<vtkIGSIOMkvSequenceIO> reader;
    reader->SetFileName(inputImageSequenceFileName.c_str());
    if (reader->Read() != IGSIO_SUCCESS)
    {
      std::cerr << "Couldn't read sequence MKV: " << inputImageSequenceFileName << std::endl;
      return EXIT_FAILURE;
    }
    vtkIGSIOTrackedFrameList* trackedFrameList = reader->GetTrackedFrameList();

    if (trackedFrameList == NULL)
    {
      std::cerr << "Unable to get trackedFrameList!" << std::endl;
      return EXIT_FAILURE;
    }

    if (trackedFrameList->GetNumberOfTrackedFrames() < 0)
    {
      std::cerr << "No frames in trackedFrameList!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (!outputImageSequenceFileName.empty())
  {
    // ******************************************************************************
    // Test writing

    int width = 10;
    int height = 10;
    int numFrames = 100;
    float fps = 10.0;
    double timeBetweenFrames = 1.0 / fps;
    double timestamp = 0.0;

    FrameSizeType frameSize = { width, height, 1 };

    vtkNew<vtkIGSIOTrackedFrameList> outputTrackedFrameList;

    igsioTrackedFrame emptyFrame;
    for (unsigned int i = outputTrackedFrameList->GetNumberOfTrackedFrames(); i < numFrames; i++)
    {
      outputTrackedFrameList->AddTrackedFrame(&emptyFrame, vtkIGSIOTrackedFrameList::ADD_INVALID_FRAME);
    }

    for (int i = 0; i < numFrames; ++i)
    {
      igsioTrackedFrame* trackedFrame = outputTrackedFrameList->GetTrackedFrame(i);
      igsioVideoFrame* videoFrame = trackedFrame->GetImageData();
      videoFrame->AllocateFrame(frameSize, VTK_UNSIGNED_CHAR, 3);

      vtkImageData* imageData = videoFrame->GetImage();
      unsigned char* imageDataScalars = (unsigned char*)imageData->GetScalarPointer();
      unsigned char blue = 255 * (i / (double)numFrames);
      for (int y = 0; y < height; ++y)
      {
        unsigned char green = 255 * (y / (double)height);
        for (int x = 0; x < width; ++x)
        {
          unsigned char red = 255 * (x / (double)width);
          imageDataScalars[0] = red;
          imageDataScalars[1] = green;
          imageDataScalars[2] = blue;
          imageDataScalars += 3;
        }
      }
      imageDataScalars = (unsigned char*)imageData->GetScalarPointer();
      timestamp += timeBetweenFrames;
      trackedFrame->SetTimestamp(timestamp);
    }

    vtkNew<vtkIGSIOMkvSequenceIO> writer;
    writer->SetFileName(outputImageSequenceFileName.c_str());
    writer->SetTrackedFrameList(outputTrackedFrameList.GetPointer());
    if (!writer->Write())
    {
      std::cerr << "Could not write trackedFrameList!" << std::endl;
      return EXIT_FAILURE;
    }
    writer->Close();
  }

  std::cout << "vtkIGSIOMkvSequenceIOTest completed successfully!" << std::endl;
  return EXIT_SUCCESS;
}
