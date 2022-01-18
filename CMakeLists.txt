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

root_generate_dictionary(${PROJECT_NAME}Dict RooTwoSidedCBShape.h RooExpandedFitResult.h LINKDEF LinkDef.h)

add_library(UtilsLib SHARED
    src/Colors.cpp
    src/Utils.cpp 
    src/Utils_WS.cpp
    ext/catch.cpp)

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
    src/DistributionTool.cpp
    src/AutoBinningTool.cpp
    ext/WSMakerBinning.cpp)

add_library(StatsLib SHARED
    src/WorkSpace.cpp
    src/PullsTool.cpp
    src/RankingTool.cpp
    src/WorkspaceToRoot.cpp
    ext/roofitUtils.cpp
    ext/RooExpandedFitResult.cpp
    ext/RooTwoSidedCBShape.cpp
    ext/LinkDef.h
    ${PROJECT_NAME}Dict.cxx)

file(GLOB EXAMPLES "examples/Example_*.cpp")
file(GLOB HADHADEXAMPLES "examples/HadHadSR_*.cpp")
file(GLOB THESISEXAMPLES "examples/thesis/Example_*.cpp")

add_library(ExamplesLib SHARED
    examples/Example.cpp
    ${EXAMPLES} ${HADHADEXAMPLES} ${THESISEXAMPLES})

target_link_libraries(UtilsLib PUBLIC ${ROOT_LIBRARIES})
target_link_libraries(ConfigLib PUBLIC stdc++fs ${ROOT_LIBRARIES})
target_link_libraries(HistToolLib PUBLIC stdc++fs ${ROOT_LIBRARIES})
target_link_libraries(StatsLib PUBLIC stdc++fs Threads::Threads ${ROOT_LIBRARIES})
target_link_libraries(ExamplesLib PUBLIC UtilsLib ConfigLib HistToolLib StatsLib ${ROOT_LIBRARIES})

add_executable(run-test exec/test.cpp)
add_executable(run-main exec/main.cpp)
add_executable(run-pulls exec/run_pulls.cpp)
add_executable(run-ranking exec/run_ranking.cpp)
add_executable(run-pd exec/run_pseudodata.cpp)
add_executable(gen-pd-zcr exec/gen_pseudodata_zcr.cpp)
add_executable(ul exec/upper_limit_calc.cpp)
add_executable(ul-rootstat exec/upper_limit_calc_roostat.cpp)

add_executable(catch-test-binning exec/unit_tests/utest_binning.cpp)

target_link_libraries(run-test ConfigLib HistToolLib UtilsLib StatsLib ExamplesLib)
target_link_libraries(run-main ConfigLib HistToolLib UtilsLib StatsLib ExamplesLib)
target_link_libraries(run-pulls UtilsLib StatsLib ExamplesLib)
target_link_libraries(run-ranking UtilsLib StatsLib ExamplesLib)
target_link_libraries(run-pd UtilsLib StatsLib ExamplesLib)
target_link_libraries(gen-pd-zcr UtilsLib)
target_link_libraries(ul UtilsLib StatsLib)
target_link_libraries(ul-rootstat UtilsLib StatsLib)

target_link_libraries(catch-test-binning UtilsLib ConfigLib HistToolLib StatsLib ExamplesLib)

add_custom_target(link_target ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink run-main bbtautau-hists)