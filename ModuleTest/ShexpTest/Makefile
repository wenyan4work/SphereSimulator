include ../../MakefileInc.mk

INCLUDE_DIRS = $(Trilinos_INCLUDE_DIRS) $(Trilinos_TPL_INCLUDE_DIRS) $(USERINCLUDE) -I../../
LIBRARY_DIRS = $(Trilinos_LIBRARY_DIRS) $(Trilinos_TPL_LIBRARY_DIRS) 
LIBRARIES = $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) $(USERLIB)

SRC =   testShexp.cpp ../../Sphere/Shexp.cpp ./PointDistribution.cpp ../../Util/Base64.cpp
INC =   ../../Sphere/Shexp.hpp ./PointDistribution.h
         
# Definitions
EXE =   TestShexp.X
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