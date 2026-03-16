cmake_minimum_required(VERSION 3.15)
project(crane_problem_project LANGUAGES CXX C)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# find Eigen and pygrampc
find_package(Eigen3 3.4 REQUIRED CONFIG)
find_package(pygrampc REQUIRED CONFIG)

pybind11_add_module(crane_problem MODULE Crane2D.cpp)
target_link_libraries(crane_problem PRIVATE Eigen3::Eigen pygrampc) # link against Eigen and pygrampc so the header files are found
install(TARGETS crane_problem DESTINATION .) # copy the .pyd (Windows) or .so (Linux) to the correct installation folder where Python finds the extension
