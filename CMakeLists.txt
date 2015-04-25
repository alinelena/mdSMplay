project(mdsmplay)
cmake_minimum_required(VERSION 2.8)

find_package(Boost  COMPONENTS random )
if(Boost_FOUND)
   include_directories(${Boost_INCLUDE_DIRS})

add_executable(mdsmplay source/utils.cpp source/solution.cpp source/reading.cpp source/gutenberg.cpp source/generalinputs.cpp source/collectivevariable.cpp source/atom.cpp source/main.cpp)
endif()


#add_subdirectory(source)
