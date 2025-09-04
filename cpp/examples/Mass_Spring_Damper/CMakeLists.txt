add_executable(Mass_Spring_Damper_cpp
    main_MSD.cpp
    problem_description_MSD.cpp
)
target_link_libraries(Mass_Spring_Damper_cpp
    grampc_cpp
)
