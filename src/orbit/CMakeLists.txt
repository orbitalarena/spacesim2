add_library(spacesim_orbit STATIC lambert.cpp)
target_include_directories(spacesim_orbit PUBLIC ${CMAKE_SOURCE_DIR}/include)
