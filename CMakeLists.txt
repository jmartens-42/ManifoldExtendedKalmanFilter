cmake_minimum_required(VERSION 3.16)

project(ManifoldExtendedKalmanFilter)

find_package(Eigen3 REQUIRED)
# find_package(ManifoldExtendedKalmanFilter REQUIRED)

set(SRC_FILES 
    src/attitude_kalman_9dof.cpp
    src/tools.cpp
)

add_library(${PROJECT_NAME} ${SRC_FILES})

# target_compile_definitions(${PROJECT_NAME}
#     PUBLIC 
#         EIGEN_NO_MALLOC
#         EIGEN_DONT_PARALLELIZE
#         EIGEN_NO_DEBUG
#         EIGEN_NO_AUTOMATIC_RESIZING
#         EIGEN_NO_IO
#         EIGEN_DONT_VECTORIZE
#
# )

target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/>
)

target_link_libraries(${PROJECT_NAME} PUBLIC
    Eigen3::Eigen
)

install(DIRECTORY include
    DESTINATION ${CMAKE_INSTALL_PREFIX}
)

install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}Targets
    INCLUDES DESTINATION include
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)
# export(
#     TARGETS ${PROJECT_NAME}
#     # NAMESPACE ${PROJECT_NAME}::
#     FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake"
# )
include(GNUInstallDirs)
install(EXPORT
    ${PROJECT_NAME}Targets
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
    # NAMESPACE ${PROJECT_NAME}::
)

include(CMakePackageConfigHelpers)
configure_package_config_file( 
  "cmake/${PROJECT_NAME}Config.cmake.in" 
  "${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  PATH_VARS
    CMAKE_INSTALL_LIBDIR
  )

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  VERSION 1.0.0
  COMPATIBILITY SameMajorVersion
  )

### Install Config and ConfigVersion files
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

