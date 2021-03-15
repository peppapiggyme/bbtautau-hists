cmake_minimum_required(VERSION 3.2)

project(bbtautau-hists LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -Wextra -g")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
SET(CMAKE_C_ARCHIVE_FINISH   "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
SET(CMAKE_CXX_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")

find_package(ROOT REQUIRED COMPONENTS Core Hist Tree RIO Net 
    Matrix Physics MathCore HistFactory RooFit RooFitCore RooFitMore RooStats)

include_directories(SYSTEM ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/examples)

add_library(UtilsLib SHARED
    src/Utils.cpp 
    src/Utils_WS.cpp)

add_library(ConfigLib SHARED
    src/Config.cpp
    src/Processes.cpp
    src/Regions.cpp
    src/Variables.cpp
    src/Systematics.cpp)

add_library(HistToolLib SHARED 
    src/HistTool.cpp
    src/DrawStackTool.cpp
    src/CompTool.cpp
    src/SystCompTool.cpp)

add_library(WorkspaceToolLib SHARED
    src/WorkspaceTool.cpp
    src/Pulls.cpp)

target_link_libraries(UtilsLib PUBLIC ${ROOT_LIBRARIES})
target_link_libraries(ConfigLib PUBLIC stdc++fs ${ROOT_LIBRARIES})
target_link_libraries(HistToolLib PUBLIC stdc++fs ${ROOT_LIBRARIES})
target_link_libraries(WorkspaceToolLib PUBLIC stdc++fs ${ROOT_LIBRARIES})

add_executable(bbtautau-hists exec/main.cpp 
    examples/ExamplesInclude.h)

target_link_libraries(bbtautau-hists 
    ConfigLib HistToolLib UtilsLib WorkspaceToolLib)
