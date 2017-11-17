
option(TESTVISP off)
add_definitions(--std=c++11)
if(TESTVISP)
	find_package(VISP)
	include_directories(${VISP_INCLUDE_DIR})
	link_libraries(${VISP_CORE} ${VISP_VISION})
	add_definitions(-DUSEVISP)
endif()

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(testtsai testtsai.cpp tsai.cpp)

add_executable(testtsai2 testtsai2.cpp tsai.cpp)
