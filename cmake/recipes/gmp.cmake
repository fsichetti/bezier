# Try to find the GNU Multiple Precision Arithmetic Library (GMP)
# See http://gmplib.org/

if(TARGET gmp::gmp)
    return()
endif()

# Download precompiled .dll on Windows
if(WIN32)
    message(STATUS "Third-party: downloading gmp + mpfr")

    include(CPM)
    CPMAddPackage(
        NAME gmp_mpfr
        URL "https://github.com/CGAL/cgal/releases/download/v5.2.1/CGAL-5.2.1-win64-auxiliary-libraries-gmp-mpfr.zip"
        DOWNLOAD_ONLY YES
    )

    # For CGAL and Cork
    set(ENV{GMP_DIR} "${gmp_mpfr_SOURCE_DIR}/gmp")
    set(ENV{MPFR_DIR} "${gmp_mpfr_SOURCE_DIR}/gmp")
else()
    # On Linux/macOS, gmp+mpfr should be installed system-wide
endif()

message(STATUS "Third-party: creating target 'gmp::gmp'")

# Find_package will look for our downloaded lib on Windows, and system-wide on Linux/macOS
find_package(GMP REQUIRED)

if(NOT TARGET gmp::gmp)
    message(FATAL_ERROR "Creation of target 'gmp::gmp' failed")
endif()
