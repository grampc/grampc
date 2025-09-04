add_executable(VTOL
    main_VTOL.c
    probfct_VTOL.c
)
target_link_libraries(VTOL grampc ${grampc_link_math})

add_executable(VTOL_fixedsize
    ${grampc_source_files}
    main_VTOL.c
    probfct_VTOL.c
)
target_include_directories(VTOL_fixedsize PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(VTOL_fixedsize PUBLIC FIXEDSIZE)
target_link_libraries(VTOL_fixedsize ${grampc_link_math})