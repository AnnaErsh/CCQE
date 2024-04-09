# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(CCQE)
set(CMAKE_CXX_STANDARD 17)

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
#find_package(ROOT REQUIRED
#  COMPONENTS Geom Physics Matrix MathCore Tree RIO EG)
#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
#include(${ROOT_USE_FILE})

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/edep-io)

execute_process(COMMAND root-config --incdir OUTPUT_VARIABLE ROOT_INCLUDE_DIR)
execute_process(COMMAND root-config --nonew --libs COMMAND tr -d \\n OUTPUT_VARIABLE ROOT_LD_DIR)

include_directories(${ROOT_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h*)

add_executable(CCQE CCQE.cc ${sources} ${headers})
target_link_libraries(CCQE ${ROOT_LD_DIR})

install(TARGETS CCQE DESTINATION bin)

# set(source
#   src/CCQE.cc
#   src/FormFactor.cc
#   src/nuCCQE.cc)

# set(include
#   include/nuCCQE.hh
#   include/FormFactor.hh
#   include/vec.h
#   include/vect.h
#   include/jednostki.h)

# ROOT_GENERATE_DICTIONARY(G__edep ${include} LINKDEF edep-io/edepsim_io_LinkDef.h)

# #---Create a shared library with geneated dictionary
# add_library(edep SHARED ${source} G__edep.cxx)
# target_link_libraries(edep ${ROOT_LIBRARIES})

# #---Create  a main program using the library
# add_executable(analysis_example analysis_example.cxx)
# target_link_libraries(analysis_example edep)