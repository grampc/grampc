add_executable(PMSM
    main_PMSM.c
    probfct_PMSM.c
)
target_link_libraries(PMSM grampc ${grampc_link_math})

add_executable(PMSM_fixedsize
    ${grampc_source_files}
    main_PMSM.c
    probfct_PMSM.c
)
target_include_directories(PMSM_fixedsize PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(PMSM_fixedsize PUBLIC FIXEDSIZE)
target_link_libraries(PMSM_fixedsize ${grampc_link_math})
