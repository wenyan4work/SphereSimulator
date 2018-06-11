include ../MakefileInc.mk

SRC =   testParSort.cpp
INC =  ../../sctl/sctl.hpp
         
# Definitions

EXE =   TestParSort.X
OBJ :=   $(SRC:.cpp=.o)


# Link rule
$(EXE):	$(OBJ) 
	$(LINK) $(OBJ) -o $(EXE) $(LINKFLAGS) $(SYSLIB) $(LIBRARY_DIRS) $(LIBRARIES)

# Compilation rules
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -c $*.cpp -o $*.o -DDEBUGVTU


all: $(EXE)

clean: 
	rm -f ./$(OBJ)
	rm -f ./$(EXE)
