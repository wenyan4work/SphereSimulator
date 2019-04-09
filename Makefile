include MakefileInc.mk

# Files
SRC =  \
SRC/Config.cpp \
SRC/main.cpp \
SRC/SphereSystem.cpp \
Sphere/Sphere.cpp \
Sphere/Shexp.cpp \
BIE/SphereSTKMobMat.cpp \
BIE/SphereSTKSHOperator.cpp \
$(SIMTOOLBOX)/Collision/CollisionSolver.cpp \
$(SIMTOOLBOX)/Collision/CPSolver.cpp \
$(SIMTOOLBOX)/Trilinos/TpetraUtil.cpp \
$(SIMTOOLBOX)/Util/Base64.cpp \
$(STKFMM)/STKFMM.cpp 

# Definitions
EXE :=   SphereSimulator.X
OBJ :=   $(SRC:.cpp=.o)

all: $(EXE)

# pull in dependency info for *existing* .o files
-include $(OBJ:.o=.d)

LFLAG = $(LINKFLAGS) $(SYSLIB) $(LIBRARY_DIRS) $(LIBRARIES)

# Link rule
$(EXE):	$(OBJ)
	$(LINK) $(OBJ) -o $(EXE) $(LFLAG)
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
	rm -vf ./$(OBJ)
	rm -vf ./$(EXE)
	rm -vf ./*.d
	rm -vf ./*~
