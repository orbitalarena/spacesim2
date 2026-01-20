add_library(spacesim2_core STATIC environment.cpp solarsystem.cpp vector.cpp)
target_include_directories(spacesim2_core PUBLIC ${PROJECT_SOURCE_DIR}/include)
