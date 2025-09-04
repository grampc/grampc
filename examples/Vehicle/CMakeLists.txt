add_executable(Vehicle
    main_VEHICLE.c
    probfct_VEHICLE.c
)
target_link_libraries(Vehicle grampc ${grampc_link_math})

add_executable(Vehicle_fixedsize
    ${grampc_source_files}
    main_VEHICLE.c
    probfct_VEHICLE.c
)
target_include_directories(Vehicle_fixedsize PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(Vehicle_fixedsize PUBLIC FIXEDSIZE)
target_link_libraries(Vehicle_fixedsize ${grampc_link_math})