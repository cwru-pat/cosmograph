
all: main.cc
	g++ main.cc -lm -lrt

clean:
	rm -f a.out
