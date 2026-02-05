COBJ = polynomes.cpp
OBJ = polynomes
CSV = pole_evolution.csv
PTHN = plot_poles.py

all: build

build: $(COBJ)
	g++ $(COBJ) -o $(OBJ)

run: $(OBJ) csvclean
	./$(OBJ)

plot:
	python3 $(PTHN)

csvclean:
	rm -rf $(CSV)

clean:
	rm -f $(OBJ) $(CSV)
