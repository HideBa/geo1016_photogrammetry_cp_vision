cmake_minimum_required(VERSION 2.9)

project(3dv_demos)
set(CMAKE_CXX_STANDARD 11)
find_package(OpenCV REQUIRED)

set(SRC_DIR	"${CMAKE_SOURCE_DIR}/src")
set(BIN_DIR	"${CMAKE_SOURCE_DIR}/bin")

file(GLOB APP_SOURCES "${SRC_DIR}/*.cpp")
foreach(app_source ${APP_SOURCES})
    string(REPLACE ".cpp" "" app_name ${app_source})
    string(REPLACE "${SRC_DIR}/" "" app_name ${app_name})
    add_executable(${app_name} ${app_source})
    target_link_libraries(${app_name} ${OpenCV_LIBS})
    install(TARGETS ${app_name} DESTINATION ${BIN_DIR})

    if (Ceres_FOUND)
        target_link_libraries(${app_name} ${CERES_LIBRARIES})
    endif ()

    # The resources directory
    target_compile_definitions(${app_name} PRIVATE "DATA_DIR=\"${CMAKE_CURRENT_LIST_DIR}/data/\""
    )

endforeach(app_source ${APP_SOURCES})
