add_test( MatrixTest.TestMatrixAddition [==[C:/University/IV/Computer graphics project/cmake-build-debug/testing.exe]==] [==[--gtest_filter=MatrixTest.TestMatrixAddition]==] --gtest_also_run_disabled_tests)
set_tests_properties( MatrixTest.TestMatrixAddition PROPERTIES WORKING_DIRECTORY [==[C:/University/IV/Computer graphics project/cmake-build-debug]==] SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
add_test( MatrixTest.TestMatrixOutOfBounds [==[C:/University/IV/Computer graphics project/cmake-build-debug/testing.exe]==] [==[--gtest_filter=MatrixTest.TestMatrixOutOfBounds]==] --gtest_also_run_disabled_tests)
set_tests_properties( MatrixTest.TestMatrixOutOfBounds PROPERTIES WORKING_DIRECTORY [==[C:/University/IV/Computer graphics project/cmake-build-debug]==] SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
set( testing_TESTS MatrixTest.TestMatrixAddition MatrixTest.TestMatrixOutOfBounds)