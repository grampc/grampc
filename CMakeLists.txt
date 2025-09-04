cmake_minimum_required(VERSION 3.15)
project(grampc)

# link_directories should treat paths as relative to the source dir
cmake_policy(SET CMP0015 NEW)

# Libraries should be created with position independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# This flag is required to compile the examples using MSVC,
# as the function fopen may lead to a warning or even an error.
if(MSVC)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()

# Check whether math library needs to be linked again (Yes with GCC, no with MSVC)
find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
    set(grampc_link_math m)
else()
    set(grampc_link_math "")
endif()

# Include header files
include_directories(
    include
    cpp/include
)

# Create list of source files for GRAMPC
list(APPEND grampc_source_files
  src/discrete.c
  src/grampc_alloc.c
  src/grampc_configreader.c
  src/grampc_erk.c
  src/grampc_fixedsize.c
  src/grampc_init.c
  src/grampc_mess.c
  src/grampc_run.c
  src/grampc_setopt.c
  src/grampc_setparam.c
  src/grampc_util.c
  src/rodas.c
  src/ruku45.c
  src/simpson.c
  src/timing.c
  src/trapezoidal.c
  src/finite_diff.c
)

# Convert relative to absolute paths for ease of reuse
list(TRANSFORM grampc_source_files PREPEND ${CMAKE_SOURCE_DIR}/)
# Compile GRAMPC as library
add_library(grampc ${grampc_source_files})

if(PYGRAMPC) # True when invoked via pip
    add_subdirectory(python)
else()
    # Add examples
    add_subdirectory(examples)

    # Add C++ interface
    add_subdirectory(cpp)
endif()

# Display warning message to notify users about the change in probfct signatures
message(WARNING 
    "Probfct signatures changed with GRAMPC 2.3. 
    The ordering of dfd*_vec and dHdxdt changed to (out, t, x, u, p, vec, param, userparam). 
    Make sure, your implementation uses the correct new signatures."
)
