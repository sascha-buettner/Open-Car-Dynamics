cmake_minimum_required(VERSION 3.12)
set(PROJECT_NAME ocd_point_mass_model_cpp)
project(${PROJECT_NAME})

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_package(ament_cmake REQUIRED)

# Header-only: point mass on the 3D track with exact first and second derivatives
add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME}
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)

install(
  DIRECTORY include/
  DESTINATION include
)
install(
  TARGETS ${PROJECT_NAME}
  EXPORT export_${PROJECT_NAME}
  DESTINATION lib/${PROJECT_NAME}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)
ament_export_targets(export_${PROJECT_NAME} HAS_LIBRARY_TARGET)

option(OCD_BUILD_MODEL_PROBES "Build optional synthetic model probes" OFF)
if(OCD_BUILD_MODEL_PROBES)
  enable_testing()
  add_executable(point_mass_derivatives_test test/derivatives.cpp)
  target_link_libraries(point_mass_derivatives_test ${PROJECT_NAME})
  add_test(NAME point_mass_derivatives_test COMMAND point_mass_derivatives_test)
endif()

ament_package()
