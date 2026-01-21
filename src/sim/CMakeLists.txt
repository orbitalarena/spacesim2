add_library(spacesim2_sim STATIC
    scenario.cpp
    sim.cpp
    expand_scenario.cpp
)

target_include_directories(spacesim2_sim PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(spacesim2_sim PUBLIC
    spacesim2_physics
    spacesim2_core
)
