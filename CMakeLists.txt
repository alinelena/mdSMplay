project(mdsmplay)
cmake_minimum_required(VERSION 2.8)
if (CMAKE_VERSION VERSION_LESS 3.1.0)
	set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
else()
	set(CMAKE_CXX_STANDARD 11)
endif()

find_package(Boost  COMPONENTS random program_options)
if(Boost_FOUND)
 include_directories(${Boost_INCLUDE_DIRS})

add_executable(mdsmplay source/utils.cpp 
source/solution.cpp 
source/reading.cpp 
source/gutenberg.cpp 
source/generalinputs.cpp 
source/collectivevariable.cpp 
source/parseCLI.cpp
source/atom.cpp source/main.cpp)
link_directories (${LIBRARY_OUTPUT_PATH})
target_link_libraries(mdsmplay ${Boost_LIBRARIES})
endif()


#add_subdirectory(source)
