# Simple GNU Makefile for the polygon-triangulation C++ program.
# Don E. Willcox, 2025

all:
	g++ -std=c++20 -g -o polygon-triangulation.exe polygon-triangulation.cpp

test: all
	./polygon-triangulation.exe test_square_10-10.csv

clean:
	rm *.exe
