INCLUDE 	= -I . -I ./Eigen -I ./Spectra
CC      	= g++
FLAGS   	= -fpermissive -std=c++17 -O3 -w -static -m64 -fopenmp -march=x86-64 -DNDEBUG -DENABLE_RELATEDNESS -msse -msse2 -msse3 -msse4.1 -msse4.2 -mlzcnt -mpopcnt
HEADER  	= vcfpop.h global.h misc.h file.h hash.h math2.h math512.h mathAVX.h mathSSE.h string2.h statistics.h matrix.h parameters.h function.h menu.h ml.h ml2.h mlbin.h mom.h mom2.h dre.h  Eigen/Eigen Spectra/SymEigsSolver.h
OBJ     	= vcfpop.o global.o misc.o file.o hash.o math2.o math512.o mathAVX.o mathSSE.o string2.o statistics.o matrix.o parameters.o function.o menu.o ml.o ml2.o mlbin.o mom.o mom2.o dre.o
BIN_PATH	= ../bin/
TARGET		= vcfpop

vcfpop: $(OBJ)
	$(CC) $(INCLUDE) $(OBJ) -m64 -fopenmp -march=x86-64 -DNDEBUG -DENABLE_RELATEDNESS -pthread -L/usr/local/lib -lm -lpthread -ldl -lz -o vcfpop

vcfpop.o: vcfpop.cpp $(HEADER)
	$(CC) $(INCLUDE) -c vcfpop.cpp $(FLAGS)

global.o: global.cpp $(HEADER)
	$(CC) $(INCLUDE) -c global.cpp $(FLAGS)

misc.o: misc.cpp $(HEADER)
	$(CC) $(INCLUDE) -c misc.cpp $(FLAGS)

file.o: file.cpp $(HEADER)
	$(CC) $(INCLUDE) -c file.cpp $(FLAGS)

hash.o: hash.cpp $(HEADER)
	$(CC) $(INCLUDE) -c hash.cpp $(FLAGS)

math2.o: math2.cpp $(HEADER)
	$(CC) $(INCLUDE) -c math2.cpp $(FLAGS)

math512.o: math512.cpp $(HEADER)
	$(CC) $(INCLUDE) -c math512.cpp $(FLAGS) -mavx -mavx2 -mfma -mavx512f -mavx512bw

mathAVX.o: mathAVX.cpp $(HEADER)
	$(CC) $(INCLUDE) -c mathAVX.cpp $(FLAGS) -mavx -mavx2 -mfma

mathSSE.o: mathSSE.cpp $(HEADER)
	$(CC) $(INCLUDE) -c mathSSE.cpp $(FLAGS)

string2.o: string2.cpp $(HEADER)
	$(CC) $(INCLUDE) -c string2.cpp $(FLAGS)

statistics.o: statistics.cpp $(HEADER)
	$(CC) $(INCLUDE) -c statistics.cpp $(FLAGS)

matrix.o: matrix.cpp $(HEADER)
	$(CC) $(INCLUDE) -c matrix.cpp $(FLAGS)

parameters.o: parameters.cpp $(HEADER)
	$(CC) $(INCLUDE) -c parameters.cpp $(FLAGS)

function.o: function.cpp $(HEADER)
	$(CC) $(INCLUDE) -c function.cpp $(FLAGS)

menu.o: menu.cpp $(HEADER)
	$(CC) $(INCLUDE) -c menu.cpp $(FLAGS)

ml.o: ml.cpp $(HEADER)
	$(CC) $(INCLUDE) -c ml.cpp $(FLAGS)

ml2.o: ml2.cpp $(HEADER)
	$(CC) $(INCLUDE) -c ml2.cpp $(FLAGS)

mlbin.o: mlbin.cpp $(HEADER)
	$(CC) $(INCLUDE) -c mlbin.cpp $(FLAGS)

mom.o: mom.cpp $(HEADER)
	$(CC) $(INCLUDE) -c mom.cpp $(FLAGS)

mom2.o: mom2.cpp $(HEADER)
	$(CC) $(INCLUDE) -c mom2.cpp $(FLAGS)

dre.o: dre.cpp $(HEADER)
	$(CC) $(INCLUDE) -c dre.cpp $(FLAGS)
	
install: $(TARGET)
	cp $(TARGET) $(BIN_PATH)
	
clean:
	rm -rf $(TARGET) $(OBJ)