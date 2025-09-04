add_executable(BallOnPlate
    main_BALL_ON_PLATE.c
    probfct_BALL_ON_PLATE.c
)
target_link_libraries(BallOnPlate grampc ${grampc_link_math})

add_executable(BallOnPlate_fixedsize
    ${grampc_source_files}
    main_BALL_ON_PLATE.c
    probfct_BALL_ON_PLATE.c
)
target_include_directories(BallOnPlate_fixedsize PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(BallOnPlate_fixedsize PUBLIC FIXEDSIZE)
target_link_libraries(BallOnPlate_fixedsize ${grampc_link_math})