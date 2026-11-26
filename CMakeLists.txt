cmake_minimum_required (VERSION 3.25)

project (calculisto_finite_difference
    VERSION 0.0.0
    DESCRIPTION "Central finite difference approximations of derivatives"
    LANGUAGES CXX
)

include (GNUInstallDirs)

# The library: header-only.
add_library (calculisto_finite_difference INTERFACE)
add_library (calculisto::finite_difference ALIAS calculisto_finite_difference)
set_target_properties (calculisto_finite_difference PROPERTIES EXPORT_NAME finite_difference)
target_include_directories (calculisto_finite_difference INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_features (calculisto_finite_difference INTERFACE cxx_std_20)

# Installation, for use with find_package (calculisto_finite_difference): by
# default only when this is the main project, not when it is used by another
# one (add_subdirectory, FetchContent).
option (CALCULISTO_FINITE_DIFFERENCE_INSTALL "Generate the installation rules" ${PROJECT_IS_TOP_LEVEL})
if (CALCULISTO_FINITE_DIFFERENCE_INSTALL)
    include (CMakePackageConfigHelpers)
    # Header-only, so independent of the architecture: under share/.
    set (config_directory "${CMAKE_INSTALL_DATADIR}/cmake/calculisto_finite_difference")
    install (TARGETS calculisto_finite_difference EXPORT calculisto_finite_difference_targets)
    install (DIRECTORY include/ DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
    install (EXPORT calculisto_finite_difference_targets
        NAMESPACE calculisto::
        FILE calculisto_finite_differenceTargets.cmake
        DESTINATION "${config_directory}"
    )
    configure_package_config_file (
        cmake/calculisto_finite_differenceConfig.cmake.in
        "${PROJECT_BINARY_DIR}/calculisto_finite_differenceConfig.cmake"
        INSTALL_DESTINATION "${config_directory}"
    )
    # Version 0.x: a new minor version may break the API.
    write_basic_package_version_file (
        "${PROJECT_BINARY_DIR}/calculisto_finite_differenceConfigVersion.cmake"
        COMPATIBILITY SameMinorVersion
        ARCH_INDEPENDENT
    )
    install (FILES
        "${PROJECT_BINARY_DIR}/calculisto_finite_differenceConfig.cmake"
        "${PROJECT_BINARY_DIR}/calculisto_finite_differenceConfigVersion.cmake"
        DESTINATION "${config_directory}"
    )
endif ()

# The tests, only when this is the main project.
if (PROJECT_IS_TOP_LEVEL)
    include (CTest)
    if (BUILD_TESTING)
        add_subdirectory (tests)
    endif ()
endif ()
