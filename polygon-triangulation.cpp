#include <iostream>
#include <array>
#include <vector>
#include <string>

/*
 * We get a polygon specified in CSV format as x, y coordinate pairs. We will
 * need to read these coordinates into a list of 2D points in a Polygon class.
 *
 * We are told the polygon is simple, but not whether the vertices are given in
 * clockwise order or in counter-clockwise order, so we will need to determine
 * the orientation and store that in the polygon class.
 *
 * We will need a list of vertex labels. Integer labels corresponding to their
 * Point index will do. This list will be accessed frequently for iteration, so
 * for optimization it could be a doubly linked list. Let's go with a std
 * vector first though.
 *
 * We will need a list of triangles we identified during triangulation. Each
 * Triangle object needs to store the integer vertex labels for its vertices.
 *
 * We will need geometry computation utility functions for the following tasks:
 * - calculate the area of a triangle from its vertex coordinates
 * - calculate the angle in radians of a 2D vector, where angle is in [0, 2*pi)
 * - calculate whether or not a point lies within a given triangle
 * - sum the area of triangles to calculate total polygon area (consider esum)
 */

struct Point {
    double x, y; 
};


struct Triangle {
    std::array<Point, 3> vertices;

    double compute_area() const
    {
        double area;
        // compute area of this triangle
        return area;
    }

    bool contains_point(Point const& point)
    {
        // returns true or false depending on whether
        // this triangle contains the given point.
        return false;
    }
};


class Polygon {
    private:
        std::vector<Point> m_points;
    public:
        Polygon(std::vector<Point>& polygon_points)
            : m_points(polygon_points)
        {}

        Point& get_point(int i) { return m_points[i]; }
        int get_num_points() { return m_points.size(); }
};


void triangulate(Polygon& polygon, std::vector<Triangle>& triangulation)
{
    // Clear the triangulation buffer
    triangulation.resize(0);

    // Create a list of vertex indices for the polygon points
    const int num_points = polygon.get_num_points();

    std::vector<int> vertex_ids;
    vertex_ids.resize(num_points);

    for (int i = 0; i < num_points; i++)
    {
        vertex_ids[i] = i;
    }

    // Iterate over the vertex indices until we have clipped away all but two
    // vertices. At that point the polygon is fully triangulated.
    while (vertex_ids.size() > 2)
    {
        // Find an Ear

        // Clip the Ear, adding to our list of Triangles

        // Remove the Ear vertex from our list of vertices
    }

    // We are finished with triangulation. The list of triangles is now given
    // by the triangulation Triangle buffer.
}


double compute_triangulated_area(std::vector<Triangle>& triangulation)
{
    // Loop over the Triangles to calculate area and accumulate the total area.
    double area = 0.0;

    for (const auto& triangle : triangulation)
    {
        area += triangle.compute_area();
    }

    return area;
}


int main()
{
    return 0;
}
