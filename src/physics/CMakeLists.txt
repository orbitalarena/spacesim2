add_library(spacesim2_physics STATIC engine.cpp orbit.cpp sun.cpp)
target_include_directories(spacesim2_physics PUBLIC ${PROJECT_SOURCE_DIR}/include)
