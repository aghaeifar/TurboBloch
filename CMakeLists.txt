cmake_minimum_required(VERSION 3.7)

project(bloch_simulator)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
#set(CMAKE_SUPPRESS_REGENERATION true)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
else()
    message("-- CMake run for msvc")
endif()

############################################################
# Create a library
############################################################

#Generate the shared library from the library sources
add_library(bloch_sim_shared SHARED bloch.cpp) 
set_target_properties(bloch_sim_shared PROPERTIES OUTPUT_NAME bloch)
add_compile_definitions(__SINGLE_PRECISION__)
#add_compile_definitions(__SEQUENTIAL__)  
#add_compile_definitions(__FASTER__)

find_package(Matlab)

if (UNIX)
    target_link_libraries(bloch_sim_shared pthread tbb dl)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64" )
    install(TARGETS bloch_sim_shared LIBRARY DESTINATION lib)    
endif()  

if (WIN32) 
    add_compile_definitions(__EXPORT_CLASS_BLOCH__) 
endif() 

if(Matlab_FOUND)
	matlab_add_mex(NAME bloch_mex LINK_TO bloch_sim_shared SRC ./MATLAB/bloch_mex.cpp R2018a) 
endif()


