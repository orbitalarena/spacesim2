add_library(spacesim2_model STATIC)

target_sources(spacesim2_model PRIVATE
    model.cpp
    rocket_model.cpp
    tle_report.cpp
    tle_spawn.cpp
    lambert_demo.cpp
)

target_include_directories(spacesim2_model PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(spacesim2_model PUBLIC
    spacesim2_core
    spacesim2_physics
)
