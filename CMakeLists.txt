cmake_minimum_required (VERSION 3.25)

project (calculisto_finite_difference
    VERSION 0.0.0
    DESCRIPTION "Central finite difference approximations of derivatives"
    LANGUAGES CXX
)

# The library: header-only.
add_library (calculisto_finite_difference INTERFACE)
add_library (calculisto::finite_difference ALIAS calculisto_finite_difference)
set_target_properties (calculisto_finite_difference PROPERTIES EXPORT_NAME finite_difference)
target_include_directories (calculisto_finite_difference INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_compile_features (calculisto_finite_difference INTERFACE cxx_std_20)

# The tests, only when this is the main project.
if (PROJECT_IS_TOP_LEVEL)
    include (CTest)
    if (BUILD_TESTING)
        add_subdirectory (tests)
    endif ()
endif ()
