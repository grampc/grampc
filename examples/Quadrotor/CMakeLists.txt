add_executable(Quadrotor
    main_QUADROTOR.c
    probfct_QUADROTOR.c
)
target_link_libraries(Quadrotor grampc ${grampc_link_math})

add_executable(Quadrotor_fixedsize
    ${grampc_source_files}
    main_QUADROTOR.c
    probfct_QUADROTOR.c
)
target_include_directories(Quadrotor_fixedsize PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(Quadrotor_fixedsize PUBLIC FIXEDSIZE)
target_link_libraries(Quadrotor_fixedsize ${grampc_link_math})