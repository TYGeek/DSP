cmake_minimum_required(VERSION 3.22)
project(DSP)

set(CMAKE_CXX_STANDARD 14)
set(TARGET_NAME DSP)
set(SOURCES main.cpp
        src/signalUnit.cpp
        src/vectorUnit.cpp
        src/utilsUnit.cpp)

add_subdirectory(lib)

add_executable(${TARGET_NAME} ${SOURCES})
target_link_libraries(${TARGET_NAME} PRIVATE sub::libplot)
target_include_directories(${TARGET_NAME} PRIVATE ./lib ./inc)
