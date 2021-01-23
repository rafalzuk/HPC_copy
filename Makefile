CXXflags = -Wall
LDLIBS = -lblas
OBJECTS = model.o burgers.o main.o

default: compile
	./myprog

all: diff advy advx burg


model.o: model.cpp MODEL.h
	mpicxx $(CXXflags) -o model.o -c model.cpp

burgers.o: burgers.cpp BURGERS.h
	mpicxx $(CXXflags) -o burgers.o -c burgers.cpp

main.o: main.cpp
	mpicxx $(CXXflags) -o main.o -c main.cpp

compile: main.o burgers.o model.o
	mpicxx -o myprog main.o burgers.o model.o $(LDLIBS)

diff: compile
	./myprog 0 0 0 1 1 500 10 50 10 50

advx: compile
	./myprog 1 0 0 0 1 500 10 50 10 50

advy: compile 
	./myprog 0 1 0 0 1 500 10 50 10 50

burg: compile 
	./myprog 1.0 0.5 1.0 0.02 1 500 10 50 10 50




diffp: compile
	mpiexec -np 2 ./myprog 0 0 0 1 1 500 10 50 10 50 2 1

advxp: compile
	mpiexec -np 2 ./myprog 1 0 0 0 1 500 10 50 10 50 2 1

advyp: compile
	mpiexec -np 2 ./myprog 0 1 0 0 1 500 10 50 10 50 2 1

burgp: compile
	mpiexec -np 2 ./myprog 1.0 0.5 1.0 0.02 1 500 10 50 10 50 2 1

custom: compile
	mpiexec -np 16 ./myprog 1.0 0.5 1.0 0.02 1 4001 10 2001 10 2001 4 4

.PHONY: clean 
clean: 
	-rm -f myprog $(OBJECTS)