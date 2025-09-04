add_executable(DAE_Integrator main_DAE_Integrator.c 
    probfct_DAE_Integrator.c)
target_link_libraries(DAE_Integrator grampc ${grampc_link_math})

add_custom_command(
    TARGET DAE_Integrator
    POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:DAE_Integrator>"
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/config_DAE_Integrator.cfg" "$<TARGET_FILE_DIR:DAE_Integrator>"
)
