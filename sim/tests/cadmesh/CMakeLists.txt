cmake_minimum_required(VERSION 3.16...3.31)
project(cadmesh_test LANGUAGES CXX)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_executable(cadmesh_test cadmesh_test.cc)
# Header-only vendored CADMesh; no ASSIMP/tetgen macros defined (STL via built-in reader).
target_include_directories(cadmesh_test PRIVATE
  ${PROJECT_SOURCE_DIR}/../../external/cadmesh)
target_link_libraries(cadmesh_test ${Geant4_LIBRARIES})
