cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 14)


#project name
project(TauLeapMethods)

#executables
add_executable(${PROJECT_NAME} main.cpp)

#download submodules
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()


add_subdirectory(src/cpp)

target_include_directories(${PROJECT_NAME} 
    PUBLIC src/hpp
    PUBLIC external/eigen
    PUBLIC external/getpot
)

target_link_directories(${PROJECT_NAME} PRIVATE src/cpp)

target_link_libraries(${PROJECT_NAME} 
    ChebyshevMethods ChemicalSystem ImplicitTauLeap ChemicalSystem MonteCarlo Parameters ProblemsList
    R-tau-ROCK SK-tau-ROCK Solver SplitStepImplicitTau SSA StabilizedTauLeapMethod tau-ROCK TauLeap
    TauLeapMethod TrapezoidalTauLeap
)

