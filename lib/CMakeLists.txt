cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(matplotlib_cpp LANGUAGES CXX)

set(PACKAGE_NAME matplotlib_cpp)

# Library target
add_library(${PACKAGE_NAME} INTERFACE)
add_library(sub::libplot ALIAS ${PACKAGE_NAME})

target_compile_features(${PACKAGE_NAME} INTERFACE cxx_std_11)

set(CMAKE_FIND_FRAMEWORK NEVER)

#set(Python3_ROOT_DIR ../../../venv)
#set(Python3_FIND_VIRTUALENV FIRST)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
target_link_libraries(${PACKAGE_NAME} INTERFACE Python3::Python Python3::Module)

message("Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
message("Python3_EXECUTABLE: ${Python3_EXECUTABLE}")
message("Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")

find_package(Python3 COMPONENTS NumPy)
if(Python3_NumPy_FOUND)
    target_link_libraries(${PACKAGE_NAME} INTERFACE Python3::NumPy)
else()
    target_compile_definitions(${PACKAGE_NAME} INTERFACE WITHOUT_NUMPY)
endif()
