cmake_minimum_required(VERSION 3.2)

project(bbtautau-hists LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -Wextra -g")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_C_ARCHIVE_FINISH   "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
set(CMAKE_CXX_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")

find_package(ROOT REQUIRED COMPONENTS Core Hist Tree RIO Net 
    Matrix Physics MathCore HistFactory RooFit RooFitCore RooFitMore RooStats)

find_package(Threads REQUIRED)

include_directories(SYSTEM ${ROOT_INCLUDE_DIR} 
    ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/ext ${CMAKE_SOURCE_DIR}/examples)

root_generate_dictionary(${PROJECT_NAME}Dict RooTwoSidedCBShape.h LINKDEF LinkDef.h)

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
    src/SystCompTool.cpp
    src/LephadTTBarCompTool.cpp
    src/AutoBinningTool.cpp)

add_library(StatsLib SHARED
    src/WorkSpace.cpp
    src/PullsTool.cpp
    src/RankingTool.cpp
    ext/roofitUtils.cpp
    ext/RooTwoSidedCBShape.cpp
    ext/LinkDef.h
    ${PROJECT_NAME}Dict.cxx)

add_library(ExamplesLib SHARED
    examples/Example.cpp
    examples/Example_klambda.cpp
    examples/Example_NonRes.cpp
    examples/Example_Pulls.cpp
    examples/Example_Ranking.cpp
    examples/Example_ttbarRew.cpp
    examples/Example_ttbarSyst.cpp
    examples/Example_WS.cpp
    examples/Example_WtDS.cpp
    examples/Example_WtGen.cpp
    examples/Example_WtOTF.cpp
    examples/Example_yield.cpp
    examples/Example_ZtautauMG.cpp
    examples/HadHadSR_v15_binning.cpp
    examples/HadHadSR_v15.cpp
    examples/HadHadSR_v15_yield.cpp)

target_link_libraries(UtilsLib PUBLIC ${ROOT_LIBRARIES})
target_link_libraries(ConfigLib PUBLIC stdc++fs ${ROOT_LIBRARIES})
target_link_libraries(HistToolLib PUBLIC stdc++fs ${ROOT_LIBRARIES})
target_link_libraries(StatsLib PUBLIC stdc++fs Threads::Threads ${ROOT_LIBRARIES})
target_link_libraries(ExamplesLib PUBLIC UtilsLib ConfigLib HistToolLib StatsLib ${ROOT_LIBRARIES})

add_executable(run-test exec/test.cpp)
add_executable(run-main exec/main.cpp)
add_executable(run-pulls exec/run_pulls.cpp)
add_executable(run-ranking exec/run_ranking.cpp)

target_link_libraries(run-test ConfigLib HistToolLib UtilsLib StatsLib ExamplesLib)
target_link_libraries(run-main ConfigLib HistToolLib UtilsLib StatsLib ExamplesLib)
target_link_libraries(run-pulls UtilsLib StatsLib ExamplesLib)
target_link_libraries(run-ranking UtilsLib StatsLib ExamplesLib)

add_custom_target(link_target ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink run-main bbtautau-hists)