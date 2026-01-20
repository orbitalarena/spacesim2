add_library(${PROJECT_NAME}_${CMAKE_CURRENT_LIST_DIR##*/} STATIC)
target_include_directories(${PROJECT_NAME}_${CMAKE_CURRENT_LIST_DIR##*/} PUBLIC ${PROJECT_SOURCE_DIR}/include)
