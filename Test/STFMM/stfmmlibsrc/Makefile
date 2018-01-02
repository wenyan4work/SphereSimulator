HOST=linux

include $(PVFMM_DIR)/MakeVariables

FC=mpifc  
FCFLAGS= $(CXXFLAGS_PVFMM) 

# Files
SRC :=   $(wildcard *.f)
OBJ :=   $(SRC:.f=.o)

all: $(OBJ)

%.o: %.f
	$(FC) $(FCFLAGS) -c $*.f -o $*.o

#
clean: 
	rm -f $(OBJS)
	
