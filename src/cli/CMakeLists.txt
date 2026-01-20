add_executable(spacesim2 main.cpp)
target_link_libraries(
    spacesim2
    spacesim2_core
    spacesim2_physics
    spacesim2_sim
    spacesim2_model
)
