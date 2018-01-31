include ../../MakefileInc.mk

INCLUDE_DIRS = $(Trilinos_INCLUDE_DIRS) $(Trilinos_TPL_INCLUDE_DIRS) $(USERINCLUDE)
LIBRARY_DIRS = $(Trilinos_LIBRARY_DIRS) $(Trilinos_TPL_LIBRARY_DIRS) 
LIBRARIES = $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) $(USERLIB)

CXX =mpicxx
LD=$(CXX)

CXXFLAGS= -std=c++11 -fopenmp -O3 -DNDEBUG
LDFLAGS= $(CXXFLAGS)

SRC =   TestSyncManager.cpp 
INC =   ../../SyncManager.hpp
         
# Definitions

EXE =   TestSyncManager.X
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