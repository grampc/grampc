# Compile GRAMPC with C++ interface
add_library(grampc_cpp
	${grampc_source_files}
	src/grampc.cpp
	src/problem_description.cpp
)
target_include_directories(grampc_cpp PUBLIC 
    $<BUILD_INTERFACE:${grampc_SOURCE_DIR}/cpp/include>
    $<INSTALL_INTERFACE:${grampc_include_dest}>
)

# installation
install(DIRECTORY include/ DESTINATION ${grampc_include_dest})
install(TARGETS grampc_cpp EXPORT grampc DESTINATION ${grampc_lib_dest})

# Add examples
add_subdirectory(examples)
