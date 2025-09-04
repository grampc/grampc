find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Eigen 3.4 is required
find_package(Eigen3 3.4 REQUIRED NO_MODULE)

include_directories(include)

pybind11_add_module(_core MODULE)
add_subdirectory(src)

target_link_libraries(_core PRIVATE Eigen3::Eigen grampc)

install(TARGETS _core DESTINATION pygrampc)