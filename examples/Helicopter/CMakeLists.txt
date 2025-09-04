add_executable(Helicopter
    main_HELICOPTER.c
    probfct_HELICOPTER.c
)
target_link_libraries(Helicopter grampc ${grampc_link_math})

add_executable(Helicopter_fixedsize
    ${grampc_source_files}
    main_HELICOPTER.c
    probfct_HELICOPTER.c
)
target_include_directories(Helicopter_fixedsize PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(Helicopter_fixedsize PUBLIC FIXEDSIZE)
target_link_libraries(Helicopter_fixedsize ${grampc_link_math})