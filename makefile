all: bwtsearch.o
	g++ bwtsearch.o -o bwtsearch
o_bwtsearch2: bwtsearch.cpp
	g++ -c bwtsearch.cpp
clean:
	rm -f *.o bwtsearch

