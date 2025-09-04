# Include header files
include_directories(
	include
)

# Compile GRAMPC with C++ interface
add_library(grampc_cpp
	${grampc_source_files}
	src/grampc.cpp
	src/problem_description.cpp
)

# Add examples
add_subdirectory(examples)
