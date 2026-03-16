find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

include_directories(include)

pybind11_add_module(_core MODULE)
add_subdirectory(src)

target_link_libraries(_core PRIVATE Eigen3::Eigen grampc)

install(TARGETS _core DESTINATION .)

# install header files and export cmake target

# add pygrampc as an interface
add_library(pygrampc INTERFACE)
target_include_directories(pygrampc INTERFACE
    $<INSTALL_INTERFACE:include/pygrampc> # specifies the location of the header files after installation
)

install(TARGETS pygrampc EXPORT pygrampc DESTINATION .)

# export target for CMake
install(
  EXPORT pygrampc
  DESTINATION lib/cmake/pygrampc
)
# install the config file so CMake later finds the lib and header files
install(FILES
    pygrampc-config.cmake
    DESTINATION lib/cmake/pygrampc
)

# install the header files of PyGRAMPC
install(
    DIRECTORY include/
    DESTINATION include/pygrampc
)

# install the header files of GRAMPC
install(
    DIRECTORY ../include/
    DESTINATION include/pygrampc
)
