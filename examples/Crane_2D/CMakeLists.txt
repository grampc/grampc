add_executable(Crane_2D
    main_CRANE_2D.c
    probfct_CRANE_2D.c
)
target_link_libraries(Crane_2D grampc ${grampc_link_math})

add_executable(Crane_2D_fixedsize
    ${grampc_source_files}
    main_CRANE_2D.c
    probfct_CRANE_2D.c
)
target_include_directories(Crane_2D_fixedsize PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(Crane_2D_fixedsize PUBLIC FIXEDSIZE)
target_link_libraries(Crane_2D_fixedsize ${grampc_link_math})

add_custom_command(
    TARGET Crane_2D
    POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:Crane_2D>"
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/config_CRANE_2D.cfg" "$<TARGET_FILE_DIR:Crane_2D>"
)

add_custom_command(
    TARGET Crane_2D_fixedsize
    POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:Crane_2D_fixedsize>"
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/config_CRANE_2D.cfg" "$<TARGET_FILE_DIR:Crane_2D_fixedsize>"
)