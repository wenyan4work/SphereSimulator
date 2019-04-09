include ../MakefileInc.mk

SRC =   testDump.cpp ../../Sphere/Sphere.cpp ../../Sphere/Shexp.cpp ../../Util/Base64.cpp
INC =   ../../Sphere/Sphere.hpp ../../Sphere/Shexp.hpp ../../Util/Buffer.hpp
         
# Definitions
EXE =   TestDump.X
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
