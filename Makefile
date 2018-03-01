include ./MakefileInc.mk

# set of libraries to correctly link to the target

INCLUDE_DIRS = $(Trilinos_INCLUDE_DIRS) $(Trilinos_TPL_INCLUDE_DIRS) $(USERINCLUDE) -I$(CURDIR) 
LIBRARY_DIRS = $(Trilinos_LIBRARY_DIRS) $(Trilinos_TPL_LIBRARY_DIRS) $(USERLIB_DIRS)
LIBRARIES = $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES) $(USERLIBS)

# System-specific settings
SHELL = /bin/sh
SYSLIB =	
SIZE =	size


# Files
SRC = Config.cpp main.cpp SphereSystem.cpp \
	./Collision/CollisionSolver.cpp \
	./Collision/CPSolver.cpp \
	./Sphere/Sphere.cpp \
	./Sphere/SPHExp.cpp \
	./STKFMM/STKFMM.cpp \
	./Trilinos/TpetraUtil.cpp \
	./Util/Base64.cpp \
	./sctl/legendre_rule.cpp

# Definitions
EXE :=   SphereSimulator.X
OBJ :=   $(SRC:.cpp=.o)

all: $(EXE)

# pull in dependency info for *existing* .o files
-include $(OBJ:.o=.d)

LFLAG = $(LINKFLAGS) $(SYSLIB) $(LIBRARY_DIRS) $(LIBRARIES)

# Link rule
$(EXE):	$(OBJ)
#	$(LINK) $(OBJ)  -o $(EXE) $(subst -fopenmp, ,$(LFLAG))
	$(LINK) $(OBJ)  -o $(EXE) $(LINKFLAGS) $(SYSLIB) $(LIBRARY_DIRS) $(LIBRARIES)
	$(SIZE) $(EXE)


# use the trick from
# http://scottmcpeak.com/autodepend/autodepend.html
# handle header dependency automatically
# compile and generate dependency info;
# more complicated dependency computation, so all prereqs listed
# will also become command-less, prereq-less targets
#   sed:    strip the target (everything before colon)
#   sed:    remove any continuation backslashes
#   fmt -1: list words one per line
#   sed:    strip leading spaces
#   sed:    add trailing colons
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -c $*.cpp -o $*.o
	$(CXX) -MM $(CXXFLAGS) $(INCLUDE_DIRS) -c $*.cpp > $*.d
	@cp -f $*.d $*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp


# remove compilation products
clean: 
	rm -f ./$(OBJ)
	rm -f ./$(EXE)
	rm -f ./*.d
	rm -f ./*~
