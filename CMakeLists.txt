cmake_minimum_required(VERSION 3.7)

project(bloch_simulator)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
#set(CMAKE_SUPPRESS_REGENERATION true)

set(CMAKE_CXX_STANDARD 17)

############################################################
# Create a library
############################################################

#Generate the shared library from the library sources
add_library(bloch_sim_shared SHARED bloch.cpp)

add_compile_definitions(__MEASURE_ELAPSED_TIME__)
#add_compile_definitions(__SEQUENTIAL__)

find_package(Matlab)

if (UNIX)
    target_include_directories(bloch_sim_shared
	PRIVATE
        /usr/include/mkl/
    )
    
    target_link_libraries(bloch_sim_shared 
        libmkl_intel_lp64.a 
        libmkl_intel_thread.a
        libmkl_core.a 
        iomp5
        pthread 
        dl
    )

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64" )

    install(TARGETS bloch_sim_shared LIBRARY DESTINATION lib)    

    if(Matlab_FOUND)
        matlab_add_mex(NAME bloch_mex LINK_TO bloch_sim_shared tbb SRC ./MATLAB/bloch_mex.cpp R2018a) 
        add_custom_command(TARGET bloch_mex POST_BUILD  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:bloch_mex> ${CMAKE_CURRENT_LIST_DIR}/lib)   
    endif()

endif()  
if (WIN32) 
    add_compile_definitions(__EXPORT_CLASS_BLOCH__)

    message("MKL is installed in " "$ENV{MKLROOT}")
    target_include_directories(bloch_sim_shared
	PRIVATE
        $ENV{MKLROOT}/include
    ) 

    target_link_libraries(bloch_sim_shared 
        PRIVATE 
        $ENV{MKLROOT}/lib/intel64/mkl_intel_lp64.lib 
        $ENV{MKLROOT}/lib/intel64/mkl_intel_thread.lib
        $ENV{MKLROOT}/lib/intel64/mkl_core.lib
        $ENV{ONEAPI_ROOT}/compiler/latest/windows/compiler/lib/intel64_win/libiomp5md.lib
    ) 

    if(Matlab_FOUND)  
        matlab_add_mex(NAME bloch_mex LINK_TO bloch_sim_shared SRC ./MATLAB/bloch_mex.cpp R2018a) 
        add_custom_command(TARGET bloch_mex POST_BUILD  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:bloch_mex> ${CMAKE_CURRENT_LIST_DIR}/lib)   
    endif()

endif() 

file(MAKE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib)

#set_target_properties(bloch_sim_shared PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_LIST_DIR}/lib)
set_target_properties(bloch_sim_shared PROPERTIES OUTPUT_NAME bloch_sim)
add_custom_command(TARGET bloch_sim_shared POST_BUILD  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:bloch_sim_shared> ${CMAKE_CURRENT_LIST_DIR}/lib)   
