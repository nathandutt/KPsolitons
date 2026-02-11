CXX = g++
CXXFLAGS = -Wall -std=c++20 -g

TARGET = main

SRC = main.cpp matrix.cpp polynomial.cpp cchain.cpp cnumber.cpp

OBJ = $(SRC:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
