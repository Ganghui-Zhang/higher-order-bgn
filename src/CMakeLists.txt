
 add_executable(higher-order-bgn_surface higher-order-bgn_surface.cc)
 target_link_dune_default_libraries(higher-order-bgn_surface)
 dune_target_enable_all_packages(higher-order-bgn_surface)
 dune_default_include_directories(higher-order-bgn_surface PRIVATE)
 target_compile_definitions("higher-order-bgn_surface" PRIVATE DUNE_GRID_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/grids/\")

add_executable(higher-order-bgn_curve higher-order-bgn_curve.cc)
target_link_dune_default_libraries(higher-order-bgn_curve)
dune_target_enable_all_packages(higher-order-bgn_curve)
dune_default_include_directories(higher-order-bgn_curve PRIVATE)
target_compile_definitions("higher-order-bgn_curve" PRIVATE DUNE_GRID_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/grids/\")

 add_executable(eoc_curve_test eoc_curve_test.cc)
 target_link_dune_default_libraries(eoc_curve_test)
 dune_target_enable_all_packages(eoc_curve_test)
 dune_default_include_directories(eoc_curve_test PRIVATE)
 target_compile_definitions("eoc_curve_test" PRIVATE DUNE_GRID_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/grids/\")

 add_executable(eoc_surface_test eoc_surface_test.cc)
 target_link_dune_default_libraries(eoc_surface_test)
 dune_target_enable_all_packages(eoc_surface_test)
 dune_default_include_directories(eoc_surface_test PRIVATE)
 target_compile_definitions("eoc_surface_test" PRIVATE DUNE_GRID_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/grids/\")

 add_executable(eoc_curve_test_ref eoc_curve_test_ref.cc)
 target_link_dune_default_libraries(eoc_curve_test_ref)
 dune_target_enable_all_packages(eoc_curve_test_ref)
 dune_default_include_directories(eoc_curve_test_ref PRIVATE)
 target_compile_definitions("eoc_curve_test_ref" PRIVATE DUNE_GRID_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/grids/\")
