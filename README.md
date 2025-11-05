# Polygon Triangulation

Don E. Willcox, 2025.

This program triangulates a simple polygon using the Ear Clipping Method.

## Specifying the Polygon as a CSV file

The simple polygon is specified as a single command line argument to this
function. It should contain comma-separated x,y pairs. One pair appearing
on each line of the CSV file.

## Building this program

PolygonTriangulation requires only C++20 and the C++ standard library.

To build this program, do:

```
$ make
```

Build tested on Linux with the following software versions:

- GCC 15.2.1
- GNU Make 4.4.1
- Linux kernel 6.17.2

## Running this program

Please run this program as './PolygonTriangulation.exe [CSV file]'

The program will write out the resulting triangulation to the console.

In addition, the program will print the total triangulated area of the polygon.

## Running Unit and Regression Tests

GitHub Actions is set up to check the build as well as run the included unit
and regressions tests for any pull request.

Here are instructions to build and run the tests manually.

Build the unit and regression tests by doing:

```
$ make tests
```

Then run unit tests:

```
$ ./RunUnitTests.exe
```

And run regression tests:

```
$ ./RunRegressionTests.exe
```

