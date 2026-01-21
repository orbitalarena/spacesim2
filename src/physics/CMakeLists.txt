add_library(spacesim2_physics STATIC)

target_sources(spacesim2_physics PRIVATE
    engine.cpp
    orbit.cpp
    sun.cpp
    rocket.cpp
)

target_include_directories(spacesim2_physics PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(spacesim2_physics PUBLIC
    spacesim2_core
)
