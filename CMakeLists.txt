cmake_minimum_required(VERSION 3.10.0)
project(cmake_testapp)
message("Building QTLSim")

set(CMAKE_C_STANDARD 11)


find_library(GSL_LIB gsl)
if(NOT GSL_LIB)
  message(FATAL_ERROR "gsl library not found")
endif()

find_library(NLOPT_LIB nlopt)
if(NOT NLOPT_LIB)
  message(FATAL_ERROR "nlopt library not found")
endif()


file(GLOB QTLSim_SRC CONFIGURE_DEPENDS "src/*.h" "src/*.c")

add_executable(QTLSim ${QTLSim_SRC})
target_link_libraries(QTLSim gsl)
target_compile_options(QTLSim PRIVATE -O3)
