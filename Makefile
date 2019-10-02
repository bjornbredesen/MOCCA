CXX ?= g++
CFLAGS = -Wall -Wconversion -O3 -fPIC
SHVER = 2

all: svm.o mocca

mocca: ./src/main.cpp ./src/motifs.hpp ./src/motifs.cpp ./src/models/baseclassifier.hpp ./src/models/baseclassifier.cpp ./src/models/sequenceclassifier.hpp ./src/models/sequenceclassifier.cpp ./src/models/features.hpp ./src/models/features.cpp ./src/models/dummypredictor.hpp ./src/models/dummypredictor.cpp ./src/models/cpredictor.hpp ./src/models/cpredictor.cpp ./src/models/seqsvm.hpp ./src/models/seqsvm.cpp ./src/models/seqlo.hpp ./src/models/seqlo.cpp ./src/models/seqdummy.hpp ./src/models/seqdummy.cpp ./src/models/svmmocca.hpp ./src/models/svmmocca.cpp ./src/sequencelist.hpp ./src/sequencelist.cpp ./src/sequences.hpp ./src/sequences.cpp ./src/validation.hpp ./src/validation.cpp ./src/config.hpp ./src/config.cpp ./src/aux.hpp ./src/aux.cpp
	$(CXX) $(CFLAGS) ./src/main.cpp ./src/aux.cpp ./src/config.cpp ./src/validation.cpp ./src/motifs.cpp ./src/sequences.cpp ./src/sequencelist.cpp ./src/models/baseclassifier.cpp ./src/models/sequenceclassifier.cpp ./src/models/features.cpp ./src/models/cpredictor.cpp ./src/models/seqsvm.cpp ./src/models/seqlo.cpp ./src/models/seqdummy.cpp ./src/models/dummypredictor.cpp ./src/models/svmmocca.cpp -std=c++0x -o ./bin/mocca ./bin/svm.o -lm -O2


svm.o: ./src/lib/libsvm-3.17/svm.cpp ./src/lib/libsvm-3.17/svm.h
	$(CXX) $(CFLAGS) -c ./src/lib/libsvm-3.17/svm.cpp -o ./bin/svm.o

clean:
	rm -f *~ ./bin/svm.o ./bin/mocca

