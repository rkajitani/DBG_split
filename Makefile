CXX = g++
CXXFLAGS = -g -std=c++0x -O3 -funroll-loops -Wall -Wno-sign-compare -fopenmp -finline-limit-50000 -lm -Dnullptr=0


PRG = platanus_allee
OBJ = main.o assemble.o scaffold.o scaffoldGraph.o gapClose.o common.o baseCommand.o seqlib.o mapper.o gapCloseOLC.o merge.o iterate.o polish.o solveDBG.o pairedDBG.o phase.o consensus.o divide.o calc_cov.o


all: $(PRG)

$(PRG): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^
.cpp.o:
	$(CXX) -o $@ -c $< $(CXXFLAGS)

clean:
	rm -f $(PRG) $(OBJ)

