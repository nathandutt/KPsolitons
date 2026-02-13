CXX = g++
CXXFLAGS = -Wall -std=c++20 -g

TARGET = twosoliton

SRC = twosoliton.cpp polynomials.cpp monomials.cpp

OBJ = $(SRC:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
