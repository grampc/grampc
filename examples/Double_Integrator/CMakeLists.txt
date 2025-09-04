add_executable(Double_Integrator main_DOUBLE_INTEGRATOR.c 
    probfct_DOUBLE_INTEGRATOR.c)
target_link_libraries(Double_Integrator grampc ${grampc_link_math})

add_executable(Double_Integrator_OCP main_OCP_DOUBLE_INTEGRATOR.c 
    probfct_DOUBLE_INTEGRATOR_OCP.c)
target_link_libraries(Double_Integrator_OCP grampc ${grampc_link_math})

add_custom_command(
    TARGET Double_Integrator
    POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:Double_Integrator>/resMPC"
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/config_Double_Integrator.cfg" "$<TARGET_FILE_DIR:Double_Integrator>"
)

add_custom_command(
    TARGET Double_Integrator_OCP
    POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:Double_Integrator_OCP>/resOCP"
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/config_Double_Integrator.cfg" "$<TARGET_FILE_DIR:Double_Integrator_OCP>"
)
