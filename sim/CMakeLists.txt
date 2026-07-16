# sim/CMakeLists.txt --- Geant4 application `mcs_sim` (S3).
# MCS kink-angle transport engine. Built against Geant4 11.3.2 (conda env
# g4highland). Output via G4AnalysisManager (ROOT ntuple, MT-merged).

cmake_minimum_required(VERSION 3.16...3.31)
project(mcs_sim LANGUAGES CXX)

# Geant4 (with all UI/vis components available in the env).
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${PROJECT_SOURCE_DIR}/include)
# Header-only vendored CADMesh (ASCII STL -> G4TessellatedSolid); no ASSIMP/tetgen.
include_directories(${PROJECT_SOURCE_DIR}/external/cadmesh)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(mcs_sim ${sources} ${headers})
target_link_libraries(mcs_sim ${Geant4_LIBRARIES})

# Copy macros next to the binary for convenience.
file(COPY ${PROJECT_SOURCE_DIR}/macros DESTINATION ${PROJECT_BINARY_DIR})

install(TARGETS mcs_sim DESTINATION bin)
