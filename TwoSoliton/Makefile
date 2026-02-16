CXX = g++
CXXFLAGS = -Wall -std=c++20 -g

TARGET = twosoliton

SRC = twosoliton.cpp pointevolution.cpp basicpolynomial.cpp polynomials.cpp monomials.cpp complexchain.cpp complexnumber.cpp

OBJ = $(SRC:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
