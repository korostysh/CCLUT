SHELL := /bin/bash


#default compiler settings
CC = g++
OPT = -O3 -g -std=c++0x
LDFLAGS = -lm

# set pattern conversion name
Gen_EXE = pat_gen
Gen_SRC = main.cc
Gen_OBJ = Test.o

# compilation for runs
all: clean track generate link

 track : CCLUT.cc CCLUT.h
	$(CC) $(OPT) -c CCLUT.cc -o CCLUT.o

 generate : main.cc CCLUT.cc CCLUT.h
	$(CC) $(OPT) -c main.cc -o main.o
#	$(CC)  $(OPT) $(Gen_SRC) -o $(Gen_EXE) $(LDFLAGS)

 link : CCLUT.o main.o
	$(CC) main.o CCLUT.o -o $(Gen_EXE) -lm
	/bin/rm -rf *.o

clean:
	rm -rf *.o $(Gen_EXE)
