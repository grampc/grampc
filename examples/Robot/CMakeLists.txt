add_executable(Robot main_OCP_ROBOT.c 
    probfct_ROBOT.c)
target_link_libraries(Robot grampc ${grampc_link_math})

add_custom_command(
    TARGET Robot
    POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:Robot>/resOCP"
)
