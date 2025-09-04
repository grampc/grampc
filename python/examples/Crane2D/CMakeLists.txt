cmake_minimum_required(VERSION 3.15)
project(crane_problem_project LANGUAGES CXX C)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Eigen 3.4 is required
find_package(Eigen3 3.4 REQUIRED NO_MODULE)

# Path to <grampc_root>
set(GRAMPC_ROOT
    ../../../
)

include_directories(
    ${GRAMPC_ROOT}include # include directory of GRAMPC
    ${GRAMPC_ROOT}python/include # include directory of PyGRAMPC
)

pybind11_add_module(crane_problem MODULE Crane2D.cpp)
target_link_libraries(crane_problem PRIVATE Eigen3::Eigen) # linking against Eigen is necessary
install(TARGETS crane_problem DESTINATION .) # copy the .pyd (Windows) or .so (Linux) to the correct installation folder where Python finds the extension
