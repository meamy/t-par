FLAGS = -I/opt/local/include -Wall -O3 -std=c++1y
OBJS = partition.o util.o circuit.o main.o
CXX = clang++

all: $(OBJS)
	$(CXX) $(FLAGS) -o t-par $(OBJS)

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
