CXX = g++
CXXFLAGS = -Wall -std=c++17

TARGET = a.out
OBJS = main.o matrix_io.o gauss_solver.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

main.o: main.cpp matrix_io.h gauss_solver.h
matrix_io.o: matrix_io.cpp matrix_io.h
gauss_solver.o: gauss_solver.cpp gauss_solver.h

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
