
cmake_minimum_required(VERSION 3.8)

project(GadgetMC)

find_package(ROOT REQUIRED COMPONENTS Core RIO Net)

set(sources 
    simulate.cpp
    GadgetUtils.cpp
)

set(include 
    GadgetUtils.hpp
)

add_executable(simulate ${sources} ${include})

target_include_directories(simulate PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" "${ROOT_INCLUDE_DIRS}")

target_link_libraries(simulate PUBLIC ROOT::Core "${ROOT_LIBRARIES}")
