add_library(spacesim2_model STATIC)

target_sources(spacesim2_model PRIVATE
    model.cpp
    tle_spawn.cpp
    tle_report.cpp
    rocket_model.cpp
)

target_include_directories(spacesim2_model PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(spacesim2_model PUBLIC
    spacesim2_sim
    spacesim2_physics
    spacesim2_core
)
