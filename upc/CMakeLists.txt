cmake_minimum_required(VERSION 3.14)
project(upc LANGUAGES CXX)

find_package(UPCXX REQUIRED)

# Build the upc executable
add_executable(upc main.cpp serial.cpp common.h)
target_link_libraries(upc PRIVATE UPCXX::upcxx)
target_compile_definitions(upc PRIVATE)

# Copy the job scripts
# configure_file(job-upc job-upc COPYONLY)
