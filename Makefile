# Simple GNU Makefile for the PolygonTriangulation C++ program.
# Don E. Willcox, 2025

all:
	g++ -std=c++20 -g -o PolygonTriangulation.exe *.cpp -I.

test: all
	./PolygonTriangulation.exe test_square_10-10.csv
	./PolygonTriangulation.exe simple_concave_poly.csv
	./PolygonTriangulation.exe concave_poly.csv

clean:
	rm *.exe
