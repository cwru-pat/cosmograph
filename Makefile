
all: main.cc
	g++ -march=native -ffast-math -O3  main.cc -lm -lrt

clean:
	rm -f a.out
