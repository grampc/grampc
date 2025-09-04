add_executable(NLChain_4 main_NLChain_4.c 
    probfct_NLCHAIN_4.c)
target_link_libraries(NLChain_4 grampc ${grampc_link_math})

add_executable(NLChain_6 main_NLChain_6.c 
    probfct_NLCHAIN_6.c)
target_link_libraries(NLChain_6 grampc ${grampc_link_math})

add_executable(NLChain_8 main_NLChain_8.c 
    probfct_NLCHAIN_8.c)
target_link_libraries(NLChain_8 grampc ${grampc_link_math})

add_executable(NLChain_10 main_NLChain_10.c 
    probfct_NLCHAIN_10.c)
target_link_libraries(NLChain_10 grampc ${grampc_link_math})

add_executable(NLChain_12 main_NLChain_12.c 
    probfct_NLCHAIN_12.c)
target_link_libraries(NLChain_12 grampc ${grampc_link_math})

add_custom_command(
    TARGET NLChain_4
    POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:NLChain_4>/resNLC_4"
)

add_custom_command(
    TARGET NLChain_6
    POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:NLChain_6>/resNLC_6"
)

add_custom_command(
    TARGET NLChain_8
    POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:NLChain_8>/resNLC_8"
)

add_custom_command(
    TARGET NLChain_10
    POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:NLChain_10>/resNLC_10"
)

add_custom_command(
    TARGET NLChain_12
    POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:NLChain_12>/resNLC_12"
)