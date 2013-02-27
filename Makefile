FLAGS = -lrt -std=c++0x
OBJS = partition.o circuit.o main.o

all: $(OBJS)
	g++ $(FLAGS) -o tpar $(OBJS)

partition.o: src/partition.cpp
	g++ -c $(FLAGS) src/partition.cpp

circuit.o: src/circuit.cpp src/topt.cpp
	g++ -c $(FLAGS) src/circuit.cpp

main.o: src/main.cpp
	g++ -c $(FLAGS) src/main.cpp

clean: 
	rm *.o
