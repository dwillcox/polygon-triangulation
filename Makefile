# Simple GNU Makefile for the PolygonTriangulation C++ program.
# Don E. Willcox, 2025

SOURCES += PolygonTriangulation_Algorithm.cpp
SOURCES += PolygonTriangulation_Geometry_Math.cpp
SOURCES += PolygonTriangulation_Geometry_Shapes.cpp
SOURCES += PolygonTriangulation_Testing.cpp

all:
	g++ -std=c++20 -O3 -o PolygonTriangulation.exe PolygonTriangulation.cpp $(SOURCES) -I.
	g++ -std=c++20 -g -o RunUnitTests.exe RunUnitTests.cpp $(SOURCES) -I.

test: all
	./PolygonTriangulation.exe test_square_10-10.csv
	./PolygonTriangulation.exe simple_concave_poly.csv
	./PolygonTriangulation.exe concave_poly.csv

clean:
	rm *.exe
