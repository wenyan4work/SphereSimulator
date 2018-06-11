include ../MakefileInc.mk

SRC = testTCMAT.cpp ../../Trilinos/TpetraUtil.cpp
INC = ../../Trilinos/TpetraUtil.hpp
         
# Definitions

EXE =   TestTCMAT.X
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
