COBJ = Zcloud.cpp
OBJ = Zcloud
CSV = pole_evolution.csv
PTHN = plot_poles.py

all: build run plot

build: $(COBJ)
	g++ $(COBJ) -o $(OBJ)

run: $(OBJ) csvclean
	./Zcloud

plot:
	python3 $(PTHN)

csvclean:
	rm -rf $(CSV)

clean:
	rm -f $(OBJ) $(CSV)
