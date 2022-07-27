include_directories(../hpp ../../external/eigen ../../external/getpot)
set(CMAKE_CXX_STANDARD 14)


add_library(ChebyshevMethods ChebyshevMethods.cpp)
add_library(ChemicalSystem ChemicalSystem.cpp)
add_library(ImplicitTauLeap ImplicitTauLeap.cpp)
add_library(MonteCarlo MonteCarlo.cpp)
add_library(Parameters Parameters.cpp)
add_library(ProblemsList ProblemsList.cpp)
add_library(R-tau-ROCK R-tau-ROCK.cpp)
add_library(SK-tau-ROCK SK-tau-ROCK.cpp)
add_library(Solver Solver.cpp)
add_library(SplitStepImplicitTau SplitStepImplicitTau.cpp)
add_library(SSA SSA.cpp)
add_library(StabilizedTauLeapMethod StabilizedTauLeapMethod.cpp)
add_library(tau-ROCK tau-ROCK.cpp)
add_library(TauLeap TauLeap.cpp)
add_library(TauLeapMethod TauLeapMethod.cpp)
add_library(TrapezoidalTauLeap TrapezoidalTauLeap.cpp)

# I don't know why but gcc doesn't find these dependencies automatically'
target_link_libraries(Parameters PUBLIC ChemicalSystem ImplicitTauLeap)
target_link_libraries(R-tau-ROCK PUBLIC ChebyshevMethods)
