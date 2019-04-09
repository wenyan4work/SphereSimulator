include ../MakefileInc.mk

SRC =   testInteractionManager.cpp
INC =  ../../sctl/sctl.hpp ../../MPI/InteractionManager.hpp ../../MPI/NearInteraction.hpp ../../MPI/NearInteraction.txx
         
# Definitions

EXE =   TestInteractionManager.X
OBJ :=   $(SRC:.cpp=.o)


# Link rule
$(EXE):	$(OBJ) 
	$(LINK) $(OBJ) -o $(EXE) $(LINKFLAGS) $(SYSLIB) $(LIBRARY_DIRS) $(LIBRARIES)

# Compilation rules
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -c $*.cpp -o $*.o # -DDEBUGVTU -DDEBUGINTERACT


all: $(EXE)

clean: 
	rm -f ./$(OBJ)
	rm -f ./$(EXE)
