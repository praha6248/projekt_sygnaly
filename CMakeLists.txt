cmake_minimum_required(VERSION 3.4...3.18)
project(cmake_example)
add_subdirectory(pybind11)
pybind11_add_module(cmake_example src/main.cpp)
add_subdirectory(matplotplusplus)
add_subdirectory(AudioFile)

target_link_libraries(cmake_example PUBLIC matplot)
target_link_libraries(cmake_example PRIVATE AudioFile)

target_compile_definitions(cmake_example
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
