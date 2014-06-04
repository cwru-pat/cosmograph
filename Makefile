
all: main.cc
	g++ main.cc -lm -lrt -ffast-math -O3

clean:
	rm -f a.out
