add_executable(Crane_3D
    main_CRANE_3D.c
    probfct_CRANE_3D.c
)
target_link_libraries(Crane_3D grampc ${grampc_link_math})

add_executable(Crane_3D_fixedsize
    ${grampc_source_files}
    main_CRANE_3D.c
    probfct_CRANE_3D.c
)
target_include_directories(Crane_3D_fixedsize PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(Crane_3D_fixedsize PUBLIC FIXEDSIZE)
target_link_libraries(Crane_3D_fixedsize ${grampc_link_math})

add_custom_command(
    TARGET Crane_3D
    POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:Crane_3D>"
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/config_CRANE_3D.cfg" "$<TARGET_FILE_DIR:Crane_3D>"
)

add_custom_command(
    TARGET Crane_3D_fixedsize
    POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:Crane_3D_fixedsize>"
    COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/config_CRANE_3D.cfg" "$<TARGET_FILE_DIR:Crane_3D_fixedsize>"
)
