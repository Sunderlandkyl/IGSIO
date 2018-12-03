
#--------------------------------------------------------------------------------------------
#ADD_EXECUTABLE(igsioCommonTest igsioCommonTest.cxx )
#SET_TARGET_PROPERTIES(igsioCommonTest PROPERTIES FOLDER Tests)
#TARGET_LINK_LIBRARIES(igsioCommonTest vtkPlusCommon vtkIGSIOCommon )

#ADD_TEST(igsioCommonTest
#  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/igsioCommonTest
#  )

#--------------------------------------------------------------------------------------------
#ADD_EXECUTABLE(vtkTransformRepositoryTest vtkTransformRepositoryTest.cxx )
#SET_TARGET_PROPERTIES(vtkTransformRepositoryTest PROPERTIES FOLDER Tests)
#TARGET_LINK_LIBRARIES(vtkTransformRepositoryTest vtkPlusCommon vtkIGSIOCommon )

#ADD_TEST(vtkTransformRepositoryTest
#  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/vtkTransformRepositoryTest
#  --verbose=3
#  )
# This test prints some errors when testing error cases, therefore the output is not
# checked for the presence of ERROR or WARNING string