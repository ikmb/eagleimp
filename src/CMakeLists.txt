cmake_minimum_required(VERSION 3.9)

# Import local CMakeModules folder
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

project(eagleimp C CXX)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED 17)

#############################################################################
# C++ source files and target definition
#############################################################################
set(cpp_sources
    main.cpp
    version.cpp
    Args.cpp
    hybridsys/Buffer.h
    hybridsys/BufferAllocator.h
    hybridsys/BufferFactory.h
    hybridsys/Device.h
    hybridsys/DeviceCategory.cpp
    hybridsys/FPGA.cpp
    hybridsys/GPU.cpp
    hybridsys/Hybridsys.cpp
    hybridsys/Spinlock.h
    hybridsys/ThreadUtils.cpp
    utils.h
    MyMalloc.h
    Stopwatch.cpp
    MinMaxHeap.h
    StatusFile.h
    Datatypes.cpp
    RingBuffer.cpp
    History.h
    BooleanHistory.cpp
    HaplotypeHistory.h
    MapInterpolater.cpp
    VCFData.cpp
    PBWTInterval.h
    PBWT.cpp
    HaplotypePath.cpp
    BeamPath.cpp
    Beam.cpp
    Target.cpp
    TargetImp.cpp
    FPGAConfigurationEagleImp.cpp
    FPGAHandler.cpp
    GPUHandler.cpp
    PBWTPhaser.cpp
    PBWTImputer.cpp
    # GPU-only
    GPUKernels.cu
    GPUEngine.cu
)


#############################################################################
# Set up Git version extraction
#############################################################################
execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file(
    ${CMAKE_SOURCE_DIR}/version.h.template
    ${CMAKE_BINARY_DIR}/generated/version.h
    )

include_directories("${CMAKE_BINARY_DIR}/generated")


#############################################################################
# Configure libraries
#############################################################################

## NVIDIA
#cmake_policy(SET CMP0077 NEW)
#set(CUDA_PROPAGATE_HOST_FLAGS OFF)
#set(CUDA_TOOLKIT_ROOT_DIR "/opt/cuda")
#set(CUDA_BIN_PATH "/opt/cuda/bin")
#set(GPU_DEPLOYMENT_KIT_ROOT_DIR "/opt/cuda")
#find_package(NVML) #optional!
#find_package(CUDA 10.0) #optional!
#set(CUDA_HOST_COMPILER "/usr/bin/gcc-9")
#if(NVML_FOUND AND CUDA_FOUND)
##    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-Xcompiler;-Wall;-arch=sm_60;-gencode=arch=compute_35,code=sm_35;-gencode=arch=compute_50,code=sm_50;-gencode=arch=compute_60,code=sm_60;-gencode=arch=compute_70,code=sm_70")
#    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-Xcompiler;-Wall;-arch=sm_60;-gencode=arch=compute_60,code=sm_60;-gencode=arch=compute_70,code=sm_70")
#    set(CUDA_NVCC_FLAGS_DEBUG "${CUDA_NVCC_FLAGS_DEBUG};-O0;-g;-G")
#    set(CUDA_NVCC_FLAGS_RELEASE "${CUDA_NVCC_FLAGS_RELEASE};-O3")
#    set(CUDA_NVCC_FLAGS_RELWITHDEBINFO "${CUDA_NVCC_FLAGS_RELWITHDEBINFO};-O2;-g;-G")
#    cuda_add_executable(${PROJECT_NAME} ${cpp_sources})
##    target_compile_definitions(${PROJECT_NAME} PUBLIC USE_CUDA_GPU _FORCE_INLINES HAVE_SM_60)
#    target_compile_definitions(${PROJECT_NAME} PUBLIC USE_CUDA_GPU _FORCE_INLINES)
#    include_directories(${NVML_INCLUDE_DIR})
#    target_link_libraries(${PROJECT_NAME} ${NVML_LIBRARY})
#    message("Found CUDA ${CUDA_VERSION_STRING}. Enabling GPU support!")
#else()
#    message("CUDA or NVML not found. Compiling without GPU support!")
    add_definitions(-x c++)
    add_executable(${PROJECT_NAME} ${cpp_sources})
#endif()

# pthreads
set(CMAKE_THREAD_PREFER_PTHREAD ON)
find_package(Threads REQUIRED)
target_link_libraries(${PROJECT_NAME} ${CMAKE_THREAD_LIBS_INIT})

# OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)

# Zlib
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${ZLIB_LIBRARIES})

# HTSlib
find_package(HTSlib REQUIRED)
include_directories(${HTSlib_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${HTSlib_LIBRARIES})

# Alpha Data API
find_package(ADMXRC3) # optional!
if(ADMXRC3_FOUND)
    target_link_libraries(${PROJECT_NAME} ${ADMXRC3_LIBRARY})
    include_directories(${ADMXRC3_INCLUDE_DIR})
    target_compile_definitions(${PROJECT_NAME} PUBLIC USE_AD_FPGA)
    message("Found ADMXRC3. Enabling FPGA support!")
else()
    message("Alpha Data API not found. Compiling without FPGA support!")
endif()

# Boost
find_package(Boost 1.65.1 REQUIRED
    COMPONENTS
    filesystem
    program_options
)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})

# Intel TBB
#find_package(TBB REQUIRED)
#include_directories(${TBB_INCLUDE_DIRS})
#target_link_libraries(${PROJECT_NAME} ${TBB_LIBRARIES})
target_link_libraries(${PROJECT_NAME} tbb)
#target_compile_definitions(${PROJECT_NAME} PUBLIC ${TBB_DEFINITIONS} TBB_USE_EXCEPTIONS)

target_link_libraries(${PROJECT_NAME} rt)
target_link_libraries(${PROJECT_NAME} curl)

# additional CXX options
set(CXX_WARNINGS "-Wall -Wextra -Wshadow -Wno-deprecated-declarations -Wmissing-include-dirs -Wstrict-aliasing=2 -Wfloat-equal -Wcast-align")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_OPTIONS} ${cxx_compiler_opts} -std=c++11 ${CXX_WARNINGS} -fopenmp -march=native")
