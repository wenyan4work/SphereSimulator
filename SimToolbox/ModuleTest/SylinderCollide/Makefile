include ../MakefileInc.mk

SRC =   testCollisionSylinder.cpp ../../Sylinder/Sylinder.cpp ../../Util/Base64.cpp ../../Trilinos/TpetraUtil.cpp
INC =   ../../Sylinder/Sylinder.hpp ../../Util/Base64.hpp.hpp ../../Util/Buffer.hpp ../../Trilinos/TpetraUtil.hpp
         
# Definitions

EXE =   TestCollisionSylinder.X
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