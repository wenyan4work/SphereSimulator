include MakefileInc.mk

SRC =   testShexp.cpp ../../Sphere/Shexp.cpp PointDistribution.cpp ../../SimToolbox/Util/Base64.cpp
INC =   ../../Sphere/Shexp.hpp PointDistribution.h
         
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
