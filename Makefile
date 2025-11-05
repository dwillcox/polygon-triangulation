# Simple GNU Makefile for the PolygonTriangulation C++ program.
# Don E. Willcox, 2025

SOURCES += PolygonTriangulation_Algorithm.cpp
SOURCES += PolygonTriangulation_Geometry_Math.cpp
SOURCES += PolygonTriangulation_Geometry_Shapes.cpp
SOURCES += PolygonTriangulation_Testing.cpp

all:
	g++ -std=c++20 -O3 -o PolygonTriangulation.exe PolygonTriangulation.cpp $(SOURCES) -I.

tests:
	g++ -std=c++20 -g -o RunUnitTests.exe RunUnitTests.cpp $(SOURCES) -I.
	g++ -std=c++20 -g -o RunRegressionTests.exe RunRegressionTests.cpp $(SOURCES) -I.

clean:
	rm *.exe
