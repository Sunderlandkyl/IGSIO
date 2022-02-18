SET( TestDataDir ${IGSIOLIB_DATA_DIR}/TestImages )
SET( ConfigFilesDir ${IGSIOLIB_DATA_DIR}/ConfigFiles )

#--------------------------------------------------------------------------------------------
#ADD_EXECUTABLE(vtkStylusCalibrationTest vtkStylusCalibrationTest.cxx)
#SET_TARGET_PROPERTIES(vtkStylusCalibrationTest PROPERTIES FOLDER Tests)
#TARGET_LINK_LIBRARIES(vtkStylusCalibrationTest itkvnl itkvnl_algo vtkIGSIOCalibration vtkIGSIODataCollection )

#ADD_TEST(vtkStylusCalibrationTest
#  ${IGSIO_EXECUTABLE_OUTPUT_PATH}/vtkStylusCalibrationTest
#  --config-file=${ConfigFilesDir}/IGSIODeviceSet_fCal_Sim_PivotCalibration.xml
#  --baseline-file=${TestDataDir}/StylusCalibration.results.xml 
#  --outlier-generation-probability=0.05
#  )
#SET_TESTS_PROPERTIES(vtkStylusCalibrationTest PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR;WARNING")
# TODO: Kyle