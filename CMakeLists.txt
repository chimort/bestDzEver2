cmake_minimum_required(VERSION 3.10)
project(bestDz)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(bestDz fd.cpp)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_compile_options(bestDz PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(bestDz PRIVATE ${OpenMP_CXX_LIBRARIES})
else()
    message(STATUS "OpenMP not found. Attempting to fetch OpenMP.")

    include(FetchContent)
    FetchContent_Declare(
        openmp
        GIT_REPOSITORY https://github.com/llvm/llvm-project.git
        GIT_TAG release/12.x
    )

    FetchContent_MakeAvailable(openmp)
    
    target_compile_options(bestDz PRIVATE -fopenmp)
    target_link_libraries(bestDz PRIVATE OpenMP::OpenMP_CXX)
endif()

target_link_libraries(bestDz PRIVATE pthread)
