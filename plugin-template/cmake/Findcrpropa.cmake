# QUIMBY_INCLUDE_DIR = path to Quimby directory
# QUIMBY_LIBRARY = libQuimby.a
# QUIMBY_FOUND = true if Quimby is found

find_path(crpropa_INCLUDE_DIR CRPropa.h)
find_library(crpropa_LIBRARY crpropa)
find_path(crpropa_SHARE "share/crpropa")

set(crpropa_FOUND FALSE)
if(crpropa_INCLUDE_DIR AND crpropa_LIBRARY AND crpropa_SHARE)
    set(crpropa_FOUND TRUE)
    set(crpropa_SHARE ${crpropa_SHARE}/share/crpropa)
    MESSAGE(STATUS "CRPropa: Found!")
else()
    MESSAGE(STATUS "CRPropa: NOT Found!")    
endif()

MESSAGE(STATUS "  Include:     ${crpropa_INCLUDE_DIR}")
MESSAGE(STATUS "  Library:     ${crpropa_LIBRARY}")
MESSAGE(STATUS "  Share:       ${crpropa_SHARE}")

mark_as_advanced(crpropa_INCLUDE_DIR crpropa_LIBRARY crpropa_SHARE crpropa_FOUND)
