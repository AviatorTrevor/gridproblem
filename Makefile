CFLAGS = -g -Wall
TARGET = GridProblem

all: $(TARGET)

$(TARGET): main.o Point.o GridSolver.o
	g++ main.o Point.o GridSolver.o -o $(TARGET)

Point.o: Point.cpp
	g++ $(CFLAGS) -c Point.cpp

GridSolver.o: GridSolver.cpp
	g++ $(CFLAGS) -c GridSolver.cpp

main.o: main.cpp
	g++ $(CFLAGS) -c main.cpp

clean:
	rm *.o GridProblem