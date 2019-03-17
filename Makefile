COMPILER=g++
FLAGS=-Wall -O3
SRC=Matrix.cpp main.cpp
TARGET=decompositions

all:
	$(COMPILER) $(FLAGS) $(SRC) -o $(TARGET)

clean:
	rm -rf $(TARGET)