cmake_minimum_required(VERSION 3.14)

# set the version to whatever you want
# also change the name
project(myPlugin VERSION 1.0.0)

# should be at least 17
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# set other needed compiler flags:
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdiagnostics-color=always")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")

if(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--as-needed")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,--as-needed")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--as-needed")
  message(STATUS "Use --as-needed linker flags!")
endif(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)

option(ENABLE_PYTHON "Create python bindings per pybind11 (use pip or setup.py)" ON)
option(ENABLE_TESTING "Build tests and enable test target" ON)

option(USE_CRPROPA_BUILTIN "Use builtin CRPropa instead of in CRPROPA_PATH installed CRPropa" OFF)

# add other extra libraries by appending to these variables (see below on the example of CRPropa):
set(${CMAKE_PROJECT_NAME}_EXTRA_INCLUDES)
set(${CMAKE_PROJECT_NAME}_EXTRA_LIBRARIES)
set(${CMAKE_PROJECT_NAME}_EXTRA_SOURCES)

# build crpropa and add resulting libcrpropa.so as library
if(USE_CRPROPA_BUILTIN)
	if(NOT crpropa_POPULATED)
		# needed to fetch builtin version	
		include(FetchContent)
		
		# pulls repository specified and can even be pinned on specific commit
		FetchContent_Declare(
			crpropa
			GIT_REPOSITORY "https://github.com/CRPropa/CRPropa3.git"
			# GIT_TAG f5cf064cf60b8a4771deab6d963b01d44df13c3f
		)
		
		# contains .cmake files necessary for finding modules
		list(APPEND CMAKE_MODULE_PATH "${crpropa_SOURCE_DIR}/cmake")
		
		# set CMAKE_ARGS per set befor add_subdirectory (FORCE option overwrites existing entries):
		# set(ENABLE_GIT OFF CACHE BOOL "")
		# set(ENABLE_GALACTICMAGNETICLENS OFF CACHE BOOL "")
		# set(OMP_SCHEDULE "dynamic" CACHE STRING "")
		# set(FAST_WAVES ON CACHE BOOL "")
		# set(SIMD_EXTENSIONS "avx+fma" CACHE STRING "")
		
		FetchContent_MakeAvailable(crpropa)
	endif()

	list(APPEND ${CMAKE_PROJECT_NAME}_EXTRA_INCLUDES ${crpropa_SOURCE_DIR}/include)
	list(APPEND ${CMAKE_PROJECT_NAME}_EXTRA_INCLUDES ${CRPROPA_EXTRA_INCLUDES})
	list(APPEND ${CMAKE_PROJECT_NAME}_EXTRA_LIBRARIES crpropa)
	set(CRPROPA_SWIG_INTERFACE_FILES ${crpropa_SOURCE_DIR}/python)
else()
	# find crpropa:
	find_package(crpropa REQUIRED)

	# append header to includes, library to libraries and remember location of swig_interface files
	list(APPEND ${CMAKE_PROJECT_NAME}_EXTRA_INCLUDES ${crpropa_INCLUDE_DIR})
	list(APPEND ${CMAKE_PROJECT_NAME}_EXTRA_LIBRARIES ${crpropa_LIBRARY})
	set(CRPROPA_SWIG_INTERFACE_FILES ${crpropa_SHARE}/swig_interface)

	# googletest (provided by CRPropa3 repo, see https://code.google.com/p/googletest/wiki/FAQ
	#     Why is it not recommended use a pre-compiled copy of Google Test?)
	# is needed for better tests, however only need to be build if not using buitin option
	if(ENABLE_TESTING)
		add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/libs/gtest)
	endif(ENABLE_TESTING)
endif()

# googletest (provided, see https://code.google.com/p/googletest/wiki/FAQ
#     Why is it not recommended use a pre-compiled copy of Google Test?)
# Here are only the includes and apple fixes for gtest
if(ENABLE_TESTING)
	list(APPEND ${CMAKE_PROJECT_NAME}_EXTRA_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/libs/gtest/include)
	if(APPLE)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
	endif(APPLE)

	# temporary workaround for newer clang versions due to its handling of unwinding
	#   see: https://github.com/google/googletest/issues/3062
	if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
		if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 15)
		add_definitions(-DCRPROPA_TESTS_SKIP_EXCEPTIONS)
		endif (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 15)
	endif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
endif(ENABLE_TESTING)

# global includes, try to avoid this style of adding includes,
# if possible use target_include_directories instead
include_directories(
	"${CMAKE_CURRENT_BINARY_DIR}"
	"${CMAKE_CURRENT_SOURCE_DIR}/include"
)

# add the source files located in ${CMAKE_CURRENT_SOURCE_DIR}/src
add_library(${CMAKE_PROJECT_NAME} SHARED
	${CMAKE_CURRENT_SOURCE_DIR}/src/plugin.cpp

	${${CMAKE_PROJECT_NAME}_EXTRA_SOURCES}
)
# include any other needed include, the standard ${CMAKE_CURRENT_SOURCE_DIR}/include folder was included above
target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${${CMAKE_PROJECT_NAME}_EXTRA_INCLUDES})
# link crpropa and other needed libraries
target_link_libraries(${CMAKE_PROJECT_NAME} ${${CMAKE_PROJECT_NAME}_EXTRA_LIBRARIES})

# specify installation instructions
install(TARGETS ${CMAKE_PROJECT_NAME} DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/" DESTINATION "${CMAKE_INSTALL_PREFIX}/include" FILES_MATCHING PATTERN "*.h")

# Fix for apple devices, same as the one in CRPropa
# ----------------------------------------------------------------------------
# Fix Apple RPATH
# ----------------------------------------------------------------------------
set(CMAKE_MACOSX_RPATH 1)
option(USE_ABSOLUTE_RPATH "Add absolute rpath to all libraries and executables" ON)
if(APPLE OR USE_ABSOLUTE_RPATH)
  set(CMAKE_SKIP_BUILD_RPATH  FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  set(ABSOLUTE_RPATH "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")
  if(NOT IS_ABSOLUTE ${ABSOLUTE_RPATH})
    set(ABSOLUTE_RPATH ${CMAKE_BINARY_DIR}/${ABSOLUTE_RPATH})
  endif()

  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${ABSOLUTE_RPATH}" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    message(STATUS "Use absolute RPATH ${ABSOLUTE_RPATH}")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:${ABSOLUTE_RPATH}")
  endif("${isSystemDir}" STREQUAL "-1")
endif()

# Python interface, this part should work basically fully automated.
if(ENABLE_PYTHON)
	# find python, preferably the one in the current active virtual environment:
	set(Python_FIND_VIRTUALENV "FIRST")
	set(Python_FIND_STRATEGY "LOCATION")
	find_package(Python 3.0 COMPONENTS
		Interpreter
		Development
		NumPy
	)

	if(Python_FOUND)
		# If we found python we also search for SWIG
	 	find_package(SWIG 3.0 REQUIRED)
		include(UseSWIG)
		add_definitions(-D${CMAKE_PROJECT_NAME}_HAVE_PYTHON)
		list(APPEND ${CMAKE_PROJECT_NAME}_SWIG_DEFINES -D${CMAKE_PROJECT_NAME}_HAVE_PYTHON)

		# print Python info in detail
		if(NOT USE_CRPROPA_BUILTIN)
			message(STATUS "Python: Found!")
			message(STATUS "  version ${Python_VERSION}")
			message(STATUS "  executable: ${Python_EXECUTABLE}")
			message(STATUS "  libraries: ${Python_LIBRARIES}")
			message(STATUS "  headers: ${Python_INCLUDE_DIRS}")
			message(STATUS "  site packages: ${Python_SITELIB}")
			if(Python_Development_FOUND)
				message(STATUS "  development libraries: Found!")
			else()
				message(STATUS "  development libraries: NOT found!")
			endif()
			# numpy is required, otherwise Python_FOUND = false
			message(STATUS "NumPy: Found!")
			message(STATUS "  headers: ${Python_NumPy_INCLUDE_DIRS} (version ${Python_NumPy_VERSION})")
		endif()

		
		#  use Python_INSTALL_PACKAGE_DIR if provided; otherwise, install in Python_SITELIB
		set(Python_INSTALL_PACKAGE_DIR "${Python_SITELIB}" CACHE PATH "folder in which the python package is installed")
		message(STATUS "  package install directory: ${Python_INSTALL_PACKAGE_DIR}")

		set(CMAKE_SWIG_FLAGS -DWITHNUMPY ${CRP})
		list(APPEND ${CMAKE_PROJECT_NAME}_SWIG_DEFINES -DWITHNUMPY)
		list(APPEND ${CMAKE_PROJECT_NAME}_SWIG_DEFINES -doxygen)

		# use SWIG's builtin (maybe should be abandoned if everything will be tested and working with -builtin)
		option(ENABLE_SWIG_BUILTIN "Use SWIG builtin option" ON)

		if(ENABLE_SWIG_BUILTIN)
			set(BUILTIN "-builtin")
			list(APPEND ${CMAKE_PROJECT_NAME}_SWIG_DEFINES "-builtin")
		else()
			set(BUILTIN "")
		endif(ENABLE_SWIG_BUILTIN)

		file(GLOB_RECURSE ${CMAKE_PROJECT_NAME}_SWIG_INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/python/*.i)

		set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/python/${CMAKE_PROJECT_NAME}.i
			PROPERTIES
			CPLUSPLUS ON
			SWIG_MODULE_NAME "${CMAKE_PROJECT_NAME}"
			SWIG_FLAGS "${${CMAKE_PROJECT_NAME}_SWIG_DEFINES};-dirprot;${CMAKE_SWIG_FLAGS}"
			OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}
			INCLUDE_DIRECTORIES "${CRPROPA_SWIG_INTERFACE_FILES};${${CMAKE_PROJECT_NAME}_SWIG_INPUTS};${${CMAKE_PROJECT_NAME}_EXTRA_INCLUDES};"
		)
		
		swig_add_library(${CMAKE_PROJECT_NAME}-swig TYPE SHARED LANGUAGE python
			SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/python/${CMAKE_PROJECT_NAME}.i
		)
		
		target_include_directories(${CMAKE_PROJECT_NAME}-swig PUBLIC
			${CMAKE_CURRENT_SOURCE_DIR}/include
			${${CMAKE_PROJECT_NAME}_EXTRA_INCLUDES}
			${Python_INCLUDE_DIRS}
			${Python_NumPy_INCLUDE_DIRS}
			${CRPROPA_SWIG_INTERFACE_FILES}
		)
		
		target_link_libraries(${CMAKE_PROJECT_NAME}-swig ${CMAKE_PROJECT_NAME} ${Python_LIBRARIES} ${Python_LIBRARY})
		
		set_target_properties(${CMAKE_PROJECT_NAME}-swig PROPERTIES
			PREFIX ""
			OUTPUT_NAME "_${CMAKE_PROJECT_NAME}"
			SUFFIX ".so"
		)

		if(APPLE)
			set_target_properties(${CMAKE_PROJECT_NAME}-swig PROPERTIES
			INSTALL_RPATH "@loader_path;@loader_path/../../${PYTHON_PROJECT}/.libs"
			)
			set_property(TARGET ${CMAKE_PROJECT_NAME}-swig APPEND PROPERTY
			LINK_FLAGS "-flat_namespace -undefined suppress"
			)
		endif()

		install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/python/${CMAKE_PROJECT_NAME}" DESTINATION "${Python_INSTALL_PACKAGE_DIR}")
		install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}.py" DESTINATION "${Python_INSTALL_PACKAGE_DIR}/${CMAKE_PROJECT_NAME}")
		install(TARGETS ${CMAKE_PROJECT_NAME}-swig LIBRARY DESTINATION "${Python_INSTALL_PACKAGE_DIR}/${CMAKE_PROJECT_NAME}")
		install(FILES ${${CMAKE_PROJECT_NAME}_SWIG_INPUTS} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/${CMAKE_PROJECT_NAME}/swig_interface)

		find_program(PYBIND11_STUBGEN_EXECUTABLE NAMES "pybind11-stubgen")
		if(PYBIND11_STUBGEN_EXECUTABLE)
			message(STATUS "Found pybind11-stubgen at ${PYBIND11_STUBGEN_EXECUTABLE}")
			install(CODE "execute_process(COMMAND ${PYBIND11_STUBGEN_EXECUTABLE} --ignore-all-errors -o ${Python_SITELIB} ${CMAKE_PROJECT_NAME})")
		install(CODE "execute_process(COMMAND echo -- Generated python stubs over pybind11-stubgen for ${CMAKE_PROJECT_NAME})")
		endif()

	endif()

endif()

# the testing part. You should create proper tests for your code here:
if(ENABLE_TESTING)
	enable_testing()

	# first add the cpp tests:
	add_executable(testPluginCPP ${CMAKE_CURRENT_SOURCE_DIR}/test/testPlugin.cpp)
	target_link_libraries(testPluginCPP ${CMAKE_PROJECT_NAME} gtest gtest_main pthread)
	add_test(testPluginCPP testPluginCPP)

	if(ENABLE_PYTHON AND Python_FOUND)
		file(
			COPY
			${CMAKE_CURRENT_SOURCE_DIR}/test/
			DESTINATION 
			${CMAKE_CURRENT_BINARY_DIR}
			FILES_MATCHING PATTERN "*.py"
		)

		# add python tests here, just copy paste the add_test statements and change test name and path
		if(USE_CRPROPA_BUILTIN)
			add_test(NAME testPluginPython COMMAND bash -c "cp --target-directory ${CMAKE_CURRENT_BINARY_DIR} ${crpropa_BINARY_DIR}/crpropa.py ${crpropa_BINARY_DIR}/_crpropa.so; \
				${Python_EXECUTABLE} testPlugin.py"
			)
		else()
			add_test(NAME testPluginPython COMMAND ${Python_EXECUTABLE} testPlugin.py)
		endif()
	endif()

endif()