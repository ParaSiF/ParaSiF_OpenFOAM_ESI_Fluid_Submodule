cmake_minimum_required(VERSION 3.0)
project(dummyOF VERSION 1.0 DESCRIPTION "Multiscale Universal Interface dummyOF" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED)
find_package(MUI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
if(MUI_BASE_DIR)
  include_directories(${MUI_INCLUDE_DIR})
endif()

add_executable (PUSHER_FETCHER_1 PUSHER_FETCHER_1.cpp)
target_link_libraries(PUSHER_FETCHER_1 ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES})
