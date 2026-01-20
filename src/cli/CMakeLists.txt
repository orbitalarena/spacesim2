add_executable(spacesim2 main.cpp)

target_include_directories(spacesim2 PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(spacesim2 PRIVATE
    spacesim2_model
)
