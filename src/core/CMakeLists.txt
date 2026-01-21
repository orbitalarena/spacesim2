add_library(spacesim2_core STATIC)

target_sources(spacesim2_core PRIVATE
    output.cpp
    geodesy.cpp
    environment.cpp
    vector.cpp
    tle.cpp
    tle_to_coe.cpp
    orbit/lambert.cpp
)

target_include_directories(spacesim2_core PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)
