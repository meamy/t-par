FLAGS = -O3 -lrt -std=c++0x
OBJS = partition.o util.o circuit.o main.o
CXX = g++

all: $(OBJS)
	$(CXX) $(FLAGS) -o tpar $(OBJS)

partition.o: src/partition.cpp
	$(CXX) -c $(FLAGS) src/partition.cpp

util.o: src/util.cpp
	$(CXX) -c $(FLAGS) src/util.cpp

circuit.o: src/circuit.cpp
	$(CXX) -c $(FLAGS) src/circuit.cpp

main.o: src/main.cpp
	$(CXX) -c $(FLAGS) src/main.cpp

clean: 
	rm *.o
