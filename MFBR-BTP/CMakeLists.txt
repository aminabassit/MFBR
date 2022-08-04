cmake_minimum_required (VERSION 3.5.1)



project(demo CXX)
set(CMAKE_CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED ON)
option( BUILD_STATIC "Set to ON to include static versions of the library" OFF)

find_package(Palisade)

set( CMAKE_CXX_FLAGS ${PALISADE_CXX_FLAGS} )

include_directories( ${OPENMP_INCLUDES} )
include_directories( ${PALISADE_INCLUDE} )
include_directories( ${PALISADE_INCLUDE}/third-party/include )
include_directories( ${PALISADE_INCLUDE}/core )
include_directories( ${PALISADE_INCLUDE}/pke )

link_directories( ${PALISADE_LIBDIR} )
link_directories( ${OPENMP_LIBRARIES} )
if(BUILD_STATIC)
    set( CMAKE_EXE_LINKER_FLAGS "${PALISADE_EXE_LINKER_FLAGS} -static")
    link_libraries( ${PALISADE_STATIC_LIBRARIES} )
else()
    set( CMAKE_EXE_LINKER_FLAGS ${PALISADE_EXE_LINKER_FLAGS} )
    link_libraries( ${PALISADE_SHARED_LIBRARIES} )
endif()



# Ensure we have pthreads
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)



include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
set( SOURCESMFBR ${CMAKE_CURRENT_SOURCE_DIR}/src/mfbr.cpp ${CMAKE_CURRENT_SOURCE_DIR}/include/mfbr.hpp)






# Baseline

add_executable( expIPBaselineClearComp expIPBaselineClearComp.cpp ${SOURCESMFBR}) 
add_executable( expIPBaselineEncComp expIPBaselineEncComp.cpp ${SOURCESMFBR}) 

add_executable( expSEDBaselineClearComp expSEDBaselineClearComp.cpp ${SOURCESMFBR}) 


# MFBR

add_executable( expMFIPClearComp expMFIPClearComp.cpp ${SOURCESMFBR}) 
add_executable( expMFIPEncComp expMFIPEncComp.cpp ${SOURCESMFBR}) 

add_executable( expMFSEDClearComp expMFSEDClearComp.cpp ${SOURCESMFBR}) 



# Boddeti

add_executable( expIPBoddetiClearComp expIPBoddetiClearComp.cpp ${SOURCESMFBR}) 
add_executable( expIPBoddetiEncComp expIPBoddetiEncComp.cpp ${SOURCESMFBR}) 


# Parameters

add_executable( serialTemplateMFBR serialTemplateMFBR.cpp ${SOURCESMFBR}) 
add_executable( serialProbeMFBR serialProbeMFBR.cpp ${SOURCESMFBR}) 
add_executable( serialProbeAndTemplateBaseline serialProbeAndTemplateBaseline.cpp ${SOURCESMFBR}) 



