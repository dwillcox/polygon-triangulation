#ifndef POLY_TRI_GEOMETRY_SHAPES_H_
#define POLY_TRI_GEOMETRY_SHAPES_H_

#include <array>
#include <vector>
#include <string>

struct Coordinates {
    double x, y;

    Coordinates()
        : x(0), y(0)
    {}

    // Set the coordinate components explicitly
    Coordinates(double u, double v)
        : x(u), y(v)
    {}

    // Set the coordinate components for the vector
    // oriented from coordinate 'a' to coordinate 'b'
    Coordinates(Coordinates a, Coordinates b)
        : x(b.x-a.x),
          y(b.y-a.y)
    {}
};


struct Triangle {
    std::array<Coordinates, 3> vertices;

    Triangle()
    {
        vertices[0] = Coordinates();
        vertices[1] = Coordinates();
        vertices[2] = Coordinates();
    }

    // Store the given vertices, passed as arguments
    Triangle(const Coordinates& a, const Coordinates& b, const Coordinates& c)
    {
        vertices[0] = a;
        vertices[1] = b;
        vertices[2] = c;
    }
};


struct Polygon {
    std::vector<Coordinates> vertices;
    int orientation;

    Polygon()
        : vertices({}),
          orientation(0)
    {}

    Polygon(std::vector<Coordinates>& polygon_points)
        : vertices(polygon_points),
          orientation(0)
    {}

    int get_num_vertices() const { return vertices.size(); }

    bool read_from_csv(std::string csv_filename);
};

#endif // POLY_TRI_GEOMETRY_SHAPES_H_
