PROJECT(IGSIO)

# --------------------------------------------------------------------------
# Configure library
SET(${PROJECT_NAME}_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${IGSIOZLIB_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/IGSIOCommon
  ${CMAKE_CURRENT_SOURCE_DIR}/SequenceIO
  CACHE INTERNAL ""
  )

SET(${PROJECT_NAME}Common_SRCS
  IGSIOCommon/igsioCommon.cxx
  IGSIOCommon/igsioMath.cxx
  IGSIOCommon/igsioVideoFrame.cxx
  IGSIOCommon/igsioTrackedFrame.cxx
  IGSIOCommon/vtkIGSIOTrackedFrameList.cxx
  IGSIOCommon/vtkIGSIOTransformRepository.cxx
  IGSIOCommon/vtkIGSIORecursiveCriticalSection.cxx
  )

SET(${PROJECT_NAME}Common_HDRS
  IGSIOCommon/igsioCommon.h
  IGSIOCommon/igsioCommon.txx
  IGSIOCommon/igsioXMLUtils.h
  IGSIOCommon/igsioMath.h
  IGSIOCommon/igsioVideoFrame.h
  IGSIOCommon/igsioVideoFrame.txx
  IGSIOCommon/igsioTrackedFrame.h
  IGSIOCommon/vtkIGSIOTrackedFrameList.h
  IGSIOCommon/vtkIGSIOTransformRepository.h
  IGSIOCommon/vtkIGSIORecursiveCriticalSection.h
  )

OPTION(ENABLE_MKV_IO "Enable MKV reading/writing" ON)
IF (ENABLE_MKV_IO)
  FIND_PACKAGE(libwebm REQUIRED)
ENDIF()
  
SET(${PROJECT_NAME}SequenceIO_SRCS
  SequenceIO/vtkIGSIOSequenceIO.cxx
  SequenceIO/vtkIGSIOSequenceIOBase.cxx
  SequenceIO/vtkIGSIOMetaImageSequenceIO.cxx
  SequenceIO/vtkIGSIONrrdSequenceIO.cxx
  )
  
SET(${PROJECT_NAME}SequenceIO_HDRS
  SequenceIO/vtkIGSIOSequenceIO.h
  SequenceIO/vtkIGSIOSequenceIOBase.h
  SequenceIO/vtkIGSIOMetaImageSequenceIO.h
  SequenceIO/vtkIGSIONrrdSequenceIO.h
  )

IF (ENABLE_MKV_IO)
  ADD_DEFINITIONS(-DENABLE_MKV_IO)
  LIST(APPEND ${PROJECT_NAME}_INCLUDE_DIRS
    ${LIBWEBM_INCLUDE_DIRS}
    )
  LIST(APPEND ${PROJECT_NAME}SequenceIO_SRCS
    SequenceIO/vtkIGSIOMkvSequenceIO.cxx
    )
  LIST(APPEND ${PROJECT_NAME}SequenceIO_HDRS
    SequenceIO/vtkIGSIOMkvSequenceIO.h
    )
ENDIF()
  
SET(${PROJECT_NAME}_SRCS
  ${${PROJECT_NAME}Common_SRCS}
  ${${PROJECT_NAME}SequenceIO_SRCS}
  )

SET(${PROJECT_NAME}_HDRS
  ${${PROJECT_NAME}Common_HDRS}
  ${${PROJECT_NAME}SequenceIO_HDRS}
  )
  
SET(${PROJECT_NAME}_LIBS
  vtkCommonCore
  vtksys
  vtkIOImage
  vtkIOXMLParser
  vtkImagingCore
  ${IGSIOZLIB_LIBRARY}
  )

#-----------------------------------------------------------------------------
# Create the library
IF(MSVC OR ${CMAKE_GENERATOR} MATCHES "Xcode")
  ADD_LIBRARY(${PROJECT_NAME} ${${PROJECT_NAME}_SRCS} ${${PROJECT_NAME}_HDRS})
ELSE()
  ADD_LIBRARY(${PROJECT_NAME} ${${PROJECT_NAME}_SRCS})
ENDIF()
TARGET_LINK_LIBRARIES(${PROJECT_NAME} PUBLIC ${${PROJECT_NAME}_LIBS})

INSTALL(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}
  RUNTIME DESTINATION ${IGSIO_INSTALL_BIN_DIR} COMPONENT RuntimeLibraries
  LIBRARY DESTINATION ${IGSIO_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${IGSIO_INSTALL_LIB_DIR} COMPONENT Development)

FOREACH(p IN LISTS ${PROJECT_NAME}_INCLUDE_DIRS)
  TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${p}>)
ENDFOREACH()
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}>)

GENERATE_EXPORT_HEADER(${PROJECT_NAME})

SET_PROPERTY(GLOBAL APPEND PROPERTY IGSIO_TARGETS ${PROJECT_NAME})

SET(${PROJECT_NAME}_INCLUDE_DIRS ${${PROJECT_NAME}_INCLUDE_DIRS} CACHE INTERNAL "")

# --------------------------------------------------------------------------
# Testing
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(IGSIOCommon/Testing)
  ADD_SUBDIRECTORY(SequenceIO/Testing)
ENDIF()
