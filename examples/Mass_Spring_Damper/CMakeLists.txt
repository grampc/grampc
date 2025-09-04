add_executable(Mass_Spring_Damper main_MSD.c 
    probfct_MSD.c)
target_link_libraries(Mass_Spring_Damper grampc ${grampc_link_math})

