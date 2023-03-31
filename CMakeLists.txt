cmake_minimum_required(VERSION 3.24)
project(Bonnie_JacTest)

set(CMAKE_CXX_STANDARD 14)

find_package(Eigen3 REQUIRED NO_MODULE)
list(APPEND CMAKE_PREFIX_PATH "/opt/openrobots")
find_package(pinocchio REQUIRED)

add_executable(Bonnie_JacTest main.cpp Leg.h Leg.cpp calculate.cpp calculate.h cTypes.h cppTypes.h Pinocchio_Utilities.cpp Pinocchio_Utilities.h FileOperator.h FileOperator.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen pinocchio::pinocchio)
