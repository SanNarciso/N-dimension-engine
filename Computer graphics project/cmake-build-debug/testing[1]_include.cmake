if(EXISTS "C:/University/IV/Computer graphics project/cmake-build-debug/testing[1]_tests.cmake")
  include("C:/University/IV/Computer graphics project/cmake-build-debug/testing[1]_tests.cmake")
else()
  add_test(testing_NOT_BUILT testing_NOT_BUILT)
endif()
