FC=gfortran-mp-7
FCFLAGS= -O3 -march=native -Wa,-q  #-static # -fno-underscoring

# DO NOT add -fopenmp flag

# Files
SRC :=  sphtrans_dr.f prinm.f sphtrans.f rotproj_cmpl.f rotviajrec.f rotviarecur3.f \
    yrecursion.f legeexps.f chebexps.f dfft.f fftnext.f prini.f hkrand.f dlaran.f
OBJ :=   $(SRC:.f=.o)
EXE := TestF.X

all: $(OBJ)
	$(FC) $(FCFLAGS) -o $(EXE) ./*.o

%.o: %.f
	$(FC) $(FCFLAGS) -c $*.f -o $*.o -lgfortran


# remove compilation products
clean: 
	rm ./$(OBJ)
	rm ./$(EXE)
