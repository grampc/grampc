add_executable(Reactor_PDE
    main_REACTOR_PDE.c
    probfct_REACTOR_PDE.c
)
target_link_libraries(Reactor_PDE grampc ${grampc_link_math})

add_executable(Reactor_PDE_fixedsize
    ${grampc_source_files}
    main_REACTOR_PDE.c
    probfct_REACTOR_PDE.c
)
target_include_directories(Reactor_PDE_fixedsize PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(Reactor_PDE_fixedsize PUBLIC FIXEDSIZE)
target_link_libraries(Reactor_PDE_fixedsize ${grampc_link_math})