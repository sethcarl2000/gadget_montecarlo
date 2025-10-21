
cmake_minimum_required(VERSION 3.8)

project(GadgetMC)

find_package(ROOT REQUIRED COMPONENTS Core RIO Net)

set(sources 
    simulate.cpp
    EventType.cpp
    GadgetUtils.cpp
    LatticeFreePath.cpp
)

set(include 
    GadgetUtils.hpp
    Vec3.hpp
    Neutron.hpp
    EventType.hpp
    LatticeFreePath.hpp
    EnergyDependentCS.hpp
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

add_executable(simulate ${sources} ${include})

root_generate_dictionary(G__GadgetMC ${include} LINKDEF LinkDef.h)

#---------------------------------------------------------
#   Add Libraries
add_library(GadgetMC SHARED ${include} ${sources} G__GadgetMC.cxx)
target_include_directories(GadgetMC PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(GadgetMC PUBLIC ROOT::Core)


target_include_directories(simulate PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" "${ROOT_INCLUDE_DIRS}")

target_link_libraries(simulate PUBLIC ROOT::Core "${ROOT_LIBRARIES}")


