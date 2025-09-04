add_executable(Helicopter_Discrete main_HELICOPTER.c 
    probfct_HELICOPTER.c)
target_link_libraries(Helicopter_Discrete grampc ${grampc_link_math})
