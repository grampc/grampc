add_executable(Reactor_CSTR
    main_REACTOR_CSTR.c
    probfct_REACTOR_CSTR.c
)
target_link_libraries(Reactor_CSTR grampc ${grampc_link_math})

add_executable(Reactor_CSTR_fixedsize
    ${grampc_source_files}
    main_REACTOR_CSTR.c
    probfct_REACTOR_CSTR.c
)
target_include_directories(Reactor_CSTR_fixedsize PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(Reactor_CSTR_fixedsize PUBLIC FIXEDSIZE)
target_link_libraries(Reactor_CSTR_fixedsize ${grampc_link_math})