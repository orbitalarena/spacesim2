add_library(spacesim2_sim STATIC scenario.cpp sim.cpp)
target_include_directories(spacesim2_sim PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(spacesim2_sim spacesim2_physics spacesim2_core)
