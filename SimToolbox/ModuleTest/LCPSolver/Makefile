include ../MakefileInc.mk

SRC =   testLCPSolver.cpp  ../../Collision/CPSolver.cpp  ../../Trilinos/TpetraUtil.cpp
INC =   ../../Collision/CPSolver.hpp  ../../Trilinos/TpetraUtil.hpp
         
# Definitions
EXE =   TestLCPSolver.X
OBJ :=   $(SRC:.cpp=.o)


# Link rule
$(EXE):	$(OBJ) 
	$(LINK) $(OBJ) -o $(EXE) $(LINKFLAGS) $(SYSLIB) $(LIBRARY_DIRS) $(LIBRARIES)

# Compilation rules
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -c $*.cpp -o $*.o


all: $(EXE)

clean: 
	rm ./$(OBJ)
	rm ./$(EXE)
