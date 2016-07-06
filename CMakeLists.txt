cmake_minimum_required(VERSION 3.0)
project(cheb)
set(CMAKE_CXX_FLAGS "--std=c++14 ${CMAKE_CXX_FLAGS}")
set(SRC src/minimal.cpp)
set(CHEBDIR include/)
add_library(cheb ${SRC})
target_include_directories(cheb PUBLIC ${CHEBDIR})
set(GTEST_ROOT ${CMAKE_CURRENT_LIST_DIR}/../googletest)
add_subdirectory(${GTEST_ROOT} ${CMAKE_CURRENT_BINARY_DIR}/gtest)


set(TEST_LIB gtest_main cheb)
set(TESTS
  tests/points_test.cpp
)

foreach(test ${TESTS})
get_filename_component(filepart ${test} NAME)
string(REPLACE ".cpp" "" exec ${filepart})
add_executable(${exec} ${test})
target_link_libraries(${exec} ${TEST_LIB})
endforeach()


