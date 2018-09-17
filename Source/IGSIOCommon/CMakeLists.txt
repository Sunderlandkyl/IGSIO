PROJECT(vtkIGSIOCommon)

# --------------------------------------------------------------------------
# Configure library
SET(${PROJECT_NAME}_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  CACHE INTERNAL ""
  )

SET(${PROJECT_NAME}_SRCS
  igsioCommon.cxx
  igsioVideoFrame.cxx
  igsioTrackedFrame.cxx
  vtkTrackedFrameList.cxx
  )

SET(${PROJECT_NAME}_HDRS
  igsioCommon.h
  igsioCommon.txx
  igsioVideoFrame.h
  igsioVideoFrame.txx
  igsioTrackedFrame.h
  vtkTrackedFrameList.h
  )

SET(${PROJECT_NAME}_LIBS
  vtkCommonCore
  vtksys
  vtkIOImage
  vtkImagingCore
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
  #ADD_SUBDIRECTORY(Testing)
ENDIF()
