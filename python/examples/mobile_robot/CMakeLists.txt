cmake_minimum_required(VERSION 3.15)
project(grampc_mobile_robot LANGUAGES CXX C)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
find_package(Eigen3 3.4 REQUIRED CONFIG)
find_package(pygrampc REQUIRED CONFIG)

pybind11_add_module(${PROJECT_NAME} MODULE
    mobile_robot_problem_description.cpp
    binding.cpp
)
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen pygrampc) # link against Eigen and pygrampc so the header files are found
install(TARGETS ${PROJECT_NAME} DESTINATION .) # copy the .pyd (Windows) or .so (Linux) to the correct installation folder where Python finds the extension
