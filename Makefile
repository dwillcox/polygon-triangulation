# Simple GNU Makefile for the polygon-triangulation C++ program.
# Don E. Willcox, 2025

all:
	g++ -o polygon-triangulation.exe polygon-triangulation.cpp

test: all
	./polygon-triangulation.exe test_inputs_1
	./polygon-triangulation.exe test_inputs_1
	./polygon-triangulation.exe test_inputs_1

clean:
	rm *.exe
