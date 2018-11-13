/*==============================================================================

Copyright (c) Laboratory for Percutaneous Surgery (PerkLab)
Queen's University, Kingston, ON, Canada. All Rights Reserved.

See COPYRIGHT.txt
or http://www.slicer.org/copyright/copyright.txt for details.

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

This file was originally developed by Kyle Sunderland, PerkLab, Queen's University

==============================================================================*/

// libwebm includes
#include "mkvwriter.h"
#include "mkvparser.h"

// VTK includes
#include <vtkObjectFactory.h>

// vtkVideoIO includes
#include "vtkMKVWriter.h"

#include "igsioTrackedFrame.h"

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkMKVWriter);

class vtkMKVWriter::vtkInternal
{
public:
  vtkMKVWriter* External;

  double FrameRate;
  mkvmuxer::MkvWriter* MKVWriter;
  mkvparser::EBMLHeader* EBMLHeader;
  mkvmuxer::Segment* MKVWriteSegment;

  //----------------------------------------------------------------------------
  vtkInternal(vtkMKVWriter* external)
    : External(external)
    , MKVWriter(new mkvmuxer::MkvWriter())
    , EBMLHeader(NULL)
    , MKVWriteSegment(NULL)
  {
  }

  //----------------------------------------------------------------------------
  virtual ~vtkInternal()
  {
    delete this->MKVWriter;
    this->MKVWriter = NULL;
  }

};

//---------------------------------------------------------------------------
vtkMKVWriter::vtkMKVWriter()
  : Internal(new vtkInternal(this))
{
}

//---------------------------------------------------------------------------
vtkMKVWriter::~vtkMKVWriter()
{
  this->Close();
  delete this->Internal;
  this->Internal = NULL;
}

//---------------------------------------------------------------------------
bool vtkMKVWriter::WriteFile()
{
  if (!this->TrackedFrameList || this->TrackedFrameList->GetNumberOfTrackedFrames() < 1)
  {
    return false;
  }

  if (!this->WriteHeader())
  {
    return false;
    this->Close();
  }

  std::string trackName = ""; // TODO: !!!! this->TrackedFrameList->GetImageName();

  FrameSizeType frameSize;
  this->TrackedFrameList->GetFrameSize(frameSize);
  // TODO: !!!! this->TrackedFrameList->GetCompressedFrameSize();

  // TODO: !!!!
  // TODO: !!!!
  // TODO: !!!!
  // TODO: !!!!

  //igsioTrackedFrame* trackedFrame = this->TrackedFrameList->GetTrackedFrame(0);
  //std::map<std::string, std::string> customFields = trackedFrame->GetCustomFields();

  //int videoTrack = this->AddVideoTrack(trackName, this->TrackedFrameList->GetCodecFourCC(), frameSize[0], frameSize[1]);
  //if (videoTrack < 1)
  //{
  //  this->Close();
  //  return false;
  //}
  //this->Internal->MKVWriteSegment->CuesTrack(videoTrack);

  //std::map<std::string, int> metaDataTracks;
  //for (std::map<std::string, std::string>::iterator fieldIt = customFields.begin(); fieldIt != customFields.end(); ++fieldIt)
  //{
  //  metaDataTracks[fieldIt->first] = this->AddMetadataTrack(fieldIt->first);
  //}

  //double initialTimestamp = -1.0;
  //for (int i = 0; i < this->TrackedFrameList->GetNumberOfTrackedFrames(); ++i)
  //{
  //  igsioTrackedFrame* trackedFrame = this->TrackedFrameList->GetTrackedFrame(i);
  //  std::map<std::string, std::string> customFields = trackedFrame->GetCustomFields();

  //  double timestamp = trackedFrame->GetTimestamp();
  //  if (initialTimestamp < 0)
  //  {
  //    initialTimestamp = timestamp;
  //  }
  //  double currentTimestamp = timestamp - initialTimestamp;

  //  vtkUnsignedCharArray* compressedImage = trackedFrame->GetImageData()->GetCompressedFrameData();
  //  if (!compressedImage)
  //  {
  //    vtkErrorMacro("Error writing video, frame missing!");
  //    this->Close();
  //    return false;
  //  }

  //  unsigned long size = compressedImage->GetSize();
  //  bool isKeyFrame = trackedFrame->GetImageData()->IsKeyFrame();
  //  this->WriteEncodedVideoFrame((unsigned char*)compressedImage->GetPointer(0), size, isKeyFrame, videoTrack, currentTimestamp);
  //  for (std::map<std::string, std::string>::iterator fieldIt = customFields.begin(); fieldIt != customFields.end(); ++fieldIt)
  //  {
  //    this->WriteMetadata(fieldIt->second, metaDataTracks[fieldIt->first], currentTimestamp);
  //  }
  //}

  //this->Close();
  return true;
}

//---------------------------------------------------------------------------
bool vtkMKVWriter::WriteHeader()
{
  this->Close();
  
  this->Internal->MKVWriter = new mkvmuxer::MkvWriter();
  if (!this->Internal->MKVWriter->Open(this->Filename.c_str()))
  {
    vtkErrorMacro("Could not open file " << this->Filename << "!");
    return false;
  }

  if (!this->Internal->EBMLHeader)
  {
    this->Internal->EBMLHeader = new mkvparser::EBMLHeader();
  }

  if (!this->Internal->MKVWriteSegment)
  {
    this->Internal->MKVWriteSegment = new mkvmuxer::Segment();
  }

  if (!this->Internal->MKVWriteSegment->Init(this->Internal->MKVWriter))
  {
    vtkErrorMacro("Could not initialize MKV file segment!");
    return false;
  }

  return true;
}

//----------------------------------------------------------------------------
int vtkMKVWriter::AddVideoTrack(std::string name, std::string encodingFourCC, int width, int height, std::string language/*="und"*/)
{
  int trackNumber = this->Internal->MKVWriteSegment->AddVideoTrack(width, height, 0);
  if (trackNumber <= 0)
  {
    vtkErrorMacro("Could not create video track!");
    return false;
  }

  mkvmuxer::VideoTrack* videoTrack = (mkvmuxer::VideoTrack*)this->Internal->MKVWriteSegment->GetTrackByNumber(trackNumber);
  if (!videoTrack)
  {
    vtkErrorMacro("Could not find video track: " << trackNumber << "!");
    return false;
  }

  std::string codecId = vtkMKVUtil::FourCCToCodecId(encodingFourCC);
  videoTrack->set_codec_id(codecId.c_str());
  videoTrack->set_name(name.c_str());
  videoTrack->set_language(language.c_str());

  if (codecId == VTKVIDEOIO_MKV_UNCOMPRESSED_CODECID)
  {
    videoTrack->set_colour_space(encodingFourCC.c_str());
  }

  this->Internal->MKVWriteSegment->CuesTrack(trackNumber);

  return trackNumber;
}

//----------------------------------------------------------------------------
int vtkMKVWriter::AddMetadataTrack(std::string name, std::string language/*="und"*/)
{
  mkvmuxer::Track* const track = this->Internal->MKVWriteSegment->AddTrack(0);
  if (!track)
  {
    vtkErrorMacro("Could not create metadata track!");
    return false;
  }

  track->set_name(name.c_str());
  track->set_type(mkvparser::Track::kSubtitle);
  track->set_codec_id("S_TEXT/UTF8");
  track->set_language(language.c_str());
  return track->number();
}

//-----------------------------------------------------------------------------
bool vtkMKVWriter::WriteEncodedVideoFrame(unsigned char* encodedFrame, uint64_t size, bool isKeyFrame, int trackNumber, double timestampSeconds)
{
  if (!this->Internal->MKVWriter)
  {
    vtkErrorMacro("Header not initialized.")
    return false;
  }

  uint64_t timestampNanoSeconds = std::floor(NANOSECONDS_IN_SECOND * timestampSeconds);
  if (!this->Internal->MKVWriteSegment->AddFrame(encodedFrame, size, trackNumber, timestampNanoSeconds, isKeyFrame))
  {
    vtkErrorMacro("Error writing frame to file!");
    return false;
  }
  this->Internal->MKVWriteSegment->AddCuePoint(timestampNanoSeconds, trackNumber);
  return true;
}

//---------------------------------------------------------------------------
bool vtkMKVWriter::WriteMetadata(std::string metadata, int trackNumber, double timestampSeconds, double durationSeconds/*=1.0/NANOSECONDS_IN_SECOND*/)
{
  uint64_t timestampNanoSeconds = std::floor(NANOSECONDS_IN_SECOND * timestampSeconds);
  uint64_t durationNanoSeconds = std::floor(NANOSECONDS_IN_SECOND * durationSeconds);
  if (!this->Internal->MKVWriteSegment->AddMetadata((uint8_t*)metadata.c_str(), sizeof(char) * (metadata.length() + 1), trackNumber, timestampNanoSeconds, durationNanoSeconds))
  {
    vtkErrorMacro("Error writing metadata to file!");
    return false;
  }
  return true;
}

//---------------------------------------------------------------------------
void vtkMKVWriter::Close()
{
  if (this->Internal->MKVWriteSegment)
  {
    this->Internal->MKVWriteSegment->Finalize();
  }

  if (this->Internal->MKVWriter)
  {
    this->Internal->MKVWriter->Close();
    delete this->Internal->MKVWriter;
    this->Internal->MKVWriter = NULL;
  }

  if (this->Internal->MKVWriteSegment)
  {
    delete this->Internal->MKVWriteSegment;
    this->Internal->MKVWriteSegment = NULL;
  }
}

//---------------------------------------------------------------------------
bool vtkMKVWriter::CanWriteFile(std::string filename)
{
  std::string extension = vtksys::SystemTools::GetFilenameExtension(filename);
  extension = vtksys::SystemTools::LowerCase(extension);

  if (extension == ".mkv")
  {
    return true;
  }
  else if (extension == ".webm")
  {
    return true;
  }

  return false;
}

//---------------------------------------------------------------------------
void vtkMKVWriter::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os, indent);
}
