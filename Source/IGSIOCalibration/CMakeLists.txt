PROJECT(IGSIOCalibration)

# --------------------------------------------------------------------------
# Sources
SET(${PROJECT_NAME}_SRCS
  vtkPivotCalibrationAlgo/vtkIGSIOPivotCalibrationAlgo.cxx
  vtkLandmarkDetectionAlgo/vtkIGSIOLandmarkDetectionAlgo.cxx
  )

SET(${PROJECT_NAME}_HDRS
  vtkPivotCalibrationAlgo/vtkIGSIOPivotCalibrationAlgo.h
  vtkLandmarkDetectionAlgo/vtkIGSIOLandmarkDetectionAlgo.h
  )

SET(${PROJECT_NAME}_INCLUDE_DIRS
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/vtkPivotCalibrationAlgo
  ${CMAKE_CURRENT_SOURCE_DIR}/vtkLandmarkDetectionAlgo
  CACHE INTERNAL "" FORCE)

# --------------------------------------------------------------------------
# Build the library
SET(${PROJECT_NAME}_LIBS
  ITKCommon
  ITKIOImageBase
  ITKOptimizers
  ITKTransform
  vtkIGSIOCommon
  )
IF(TARGET MGHIO)
  LIST(APPEND ${PROJECT_NAME}_LIBS
    MGHIO
    )
ENDIF()
IF(TARGET ITKIODCMTK)
  LIST(APPEND ${PROJECT_NAME}_LIBS
    ITKIODCMTK
    )
ENDIF()

IF(${ITK_VERSION_MAJOR} GREATER 4)
  LIST(APPEND ${PROJECT_NAME}_LIBS
    ITKIOBruker
    ITKIOJPEG2000
    ITKIOMINC
    )
ENDIF()
IF(IGSIO_RENDERING_ENABLED)
  LIST(APPEND ${PROJECT_NAME}_LIBS
    vtkIGSIORendering
    )
ENDIF()

GENERATE_EXPORT_DIRECTIVE_FILE(vtk${PROJECT_NAME})

ADD_LIBRARY(vtk${PROJECT_NAME} ${${PROJECT_NAME}_SRCS} ${${PROJECT_NAME}_HDRS})
FOREACH(p IN LISTS ${PROJECT_NAME}_INCLUDE_DIRS)
  target_include_directories(vtk${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${p}>)
ENDFOREACH()
target_include_directories(vtk${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:${IGSIOLIB_INCLUDE_INSTALL}>)
TARGET_LINK_LIBRARIES(vtk${PROJECT_NAME} PUBLIC ${${PROJECT_NAME}_LIBS})
#IGSIOAddVersionInfo(vtk${PROJECT_NAME} "Library containing various calibration algorithms. Part of the IGSIO toolkit." vtk${PROJECT_NAME} vtk${PROJECT_NAME}) # TODO

# --------------------------------------------------------------------------
# Set external MS project
SET(IGSIOLIB_DEPENDENCIES ${IGSIOLIB_DEPENDENCIES} vtk${PROJECT_NAME} CACHE INTERNAL "" FORCE)
LIST(REMOVE_DUPLICATES IGSIOLIB_DEPENDENCIES)
# Add this variable to UseIGSIOLib.cmake.in INCLUDE_IGSIOLIB_MS_PROJECTS macro
SET(vcProj_vtk${PROJECT_NAME} vtk${PROJECT_NAME};${IGSIOLib_BINARY_DIR}/src/${PROJECT_NAME}/vtk${PROJECT_NAME}.vcxproj;vtkIGSIOCommon CACHE INTERNAL "" FORCE)

IF(IGSIOBUILD_BUILD_IGSIOLib_TOOLS)
  #--------------------------------------------------------------------------------------------
  ADD_EXECUTABLE( ConvertFcsvToXml
    Tools/ConvertFcsvToXml.cxx
    )
  SET_TARGET_PROPERTIES(ConvertFcsvToXml PROPERTIES FOLDER Tools)
  # Link the executable to the algo library.
  TARGET_LINK_LIBRARIES( ConvertFcsvToXml
    vtk${PROJECT_NAME}
    ITKCommon
    vtkIGSIODataCollection
    )

  #--------------------------------------------------------------------------------------------
  ADD_EXECUTABLE( ConvertXmlToFcsv
    Tools/ConvertXmlToFcsv.cxx
    )
  SET_TARGET_PROPERTIES(ConvertXmlToFcsv PROPERTIES FOLDER Tools)
  # Link the executable to the algo library.
  TARGET_LINK_LIBRARIES( ConvertXmlToFcsv
    vtk${PROJECT_NAME}
    ITKCommon
    vtkIGSIODataCollection
    )

  #--------------------------------------------------------------------------------------------
  ADD_EXECUTABLE(TemporalCalibration Tools/TemporalCalibration.cxx)
  SET_TARGET_PROPERTIES(TemporalCalibration PROPERTIES FOLDER Tools)
  TARGET_LINK_LIBRARIES(TemporalCalibration vtkIGSIOCommon vtk${PROJECT_NAME} vtkIGSIODataCollection)
  GENERATE_HELP_DOC(TemporalCalibration)

  #--------------------------------------------------------------------------------------------
  ADD_EXECUTABLE(ProbeCalibration Tools/ProbeCalibration.cxx)
  SET_TARGET_PROPERTIES(ProbeCalibration PROPERTIES FOLDER Tools)
  TARGET_LINK_LIBRARIES(ProbeCalibration itkvnl itkvnl_algo vtk${PROJECT_NAME} vtkIGSIODataCollection )
  GENERATE_HELP_DOC(ProbeCalibration)

  INSTALL(
    TARGETS
      ConvertFcsvToXml
      ConvertXmlToFcsv
      TemporalCalibration
      ProbeCalibration
    EXPORT IGSIOLib
    DESTINATION "${IGSIOLIB_BINARY_INSTALL}"
    COMPONENT RuntimeExecutables
    )
ENDIF()

# --------------------------------------------------------------------------
# Install
#

IGSIOInstallLibrary(vtk${PROJECT_NAME} ${PROJECT_NAME})
  
# --------------------------------------------------------------------------
# Testing
#
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY (Testing)
ENDIF()