SET( TestDataDir ${CMAKE_CURRENT_SOURCE_DIR}/Data )

#*************************** vtkMetaImageSequenceIOTest  ***************************
ADD_EXECUTABLE(vtkMetaImageSequenceIOTest vtkMetaImageSequenceIOTest.cxx )
SET_TARGET_PROPERTIES(vtkMetaImageSequenceIOTest PROPERTIES FOLDER Tests)
TARGET_LINK_LIBRARIES(vtkMetaImageSequenceIOTest vtkIGSIOCommon)
ADD_TEST(vtkMetaImageSequenceIOTest 
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/vtkMetaImageSequenceIOTest
  --img-seq-file=${TestDataDir}/MetaImageSequenceIOTest1.igs.mhd
  --output-img-seq-file=MetaImageSequenceIOTest1Output.igs.mha
  )
SET_TESTS_PROPERTIES(vtkMetaImageSequenceIOTest PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR;WARNING")
