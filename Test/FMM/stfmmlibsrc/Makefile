HOST=linux

FC=mpifc  
FCFLAGS= -ipo

# Files
SRC :=   $(wildcard *.f)
OBJ :=   $(SRC:.f=.o)

all: $(OBJ)

%.o: %.f
	$(FC) $(FCFLAGS) -c $*.f -o $*.o

#
clean: 
	rm -f $(OBJS)
	