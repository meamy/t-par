FLAGS = -O3 -std=c++0x
OBJS = partition.o util.o circuit.o equiv_check.o main.o
CXX = g++

all: $(OBJS)
	$(CXX) $(FLAGS) -o tpar $(OBJS) -lrt -lz3

partition.o: src/partition.cpp
	$(CXX) -c $(FLAGS) src/partition.cpp

util.o: src/util.cpp
	$(CXX) -c $(FLAGS) src/util.cpp

circuit.o: src/circuit.cpp
	$(CXX) -c $(FLAGS) src/circuit.cpp

equiv_check.o: src/equiv_check.cpp
	$(CXX) -c $(FLAGS) src/equiv_check.cpp

main.o: src/main.cpp
	$(CXX) -c $(FLAGS) src/main.cpp

clean: 
	rm *.o
