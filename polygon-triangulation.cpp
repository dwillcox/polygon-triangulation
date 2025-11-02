#include <iostream>
#include <array>
#include <vector>
#include <numbers>
#include <cmath>
#include <string>

/*
 * We get a polygon specified in CSV format as x, y coordinate pairs, one per line.
 * We will read these coordinates into a list of 2D points in a Polygon class.
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

struct Coordinates {
    double x, y; 

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

namespace VectorMath2D {
    double dot_product(Coordinates v1, Coordinates v2)
    {
        // Compute the dot product of v1 [dot] v2.
        return v1.x * v2.x + v1.y * v2.y;
    }

    double cross_product(Coordinates v1, Coordinates v2)
    {
        // Compute the cross product of v1 [cross] v2.
        //
        // Since v1 and v2 are 2D vectors, the cross product will lie
        // perpendicular to their plane and we return its
        // component along the 3rd dimension.
        return v1.x * v2.y - v1.y * v2.x;
    }

    double magnitude(Coordinates v)
    {
        // Compute the magnitude of the vector v
        return std::sqrt(dot_product(v,v));
    }

    double angle_between_vectors(Coordinates v1, Coordinates v2)
    {
        // We wish to know the angle theta subtended by vector v1 if rotated
        // to orient along vector v2.
        //
        // We require theta lie in (-pi, pi).
        //
        // We get theta as follows:
        //
        // compute the cross product v1 [cross] v2.
        double v1_cross_v2 = cross_product(v1, v2);

        // compute the angle phi between v1 and v2 using asin -> (-pi/2, pi/2).
        double sin_phi = v1_cross_v2/(magnitude(v1) * magnitude(v2));
        double phi = std::asin(sin_phi);

        // compute the dot product v1 [dot] v2.
        double v1_dot_v2 = dot_product(v1, v2);

        // project the angle phi into (-pi, pi) using the dot product sign.
        double theta = phi;

        if (v1_dot_v2 < 0.0) {
            if (phi > 0.0) {
                theta = std::numbers::pi - phi;
            } else if (phi < 0.0) {
                theta = -phi - std::numbers::pi;
            } else {
                theta = 0.0;
            }
        }

        return theta;
    }
};

struct Triangle {
    std::array<Coordinates, 3> vertices;

    Triangle(const Coordinates& a, const Coordinates& b, const Coordinates& c)
    {
        // Store the given vertices, passed as arguments
        vertices[0] = a;
        vertices[1] = b;
        vertices[2] = c;
    }
};

namespace TriangleMath
{
    double compute_triangle_area(const Triangle& triangle)
    {
        // Compute area of this triangle using the determinant
        // and the coordinates of its vertices.
        double area = 0.5 * (vertices[0].x * (vertices[1].y - vertices[2].y) +
                             vertices[1].x * (vertices[2].y - vertices[0].y) +
                             vertices[2].x * (vertices[0].y - vertices[1].y));
        return area;
    }

    bool triangle_contains_coordinates(const Triangle& triangle, const Coordinates& point)
    {
        // Returns true or false depending on whether
        // this triangle contains the given point's coordinates.
        //
        // This triangle contains the given point if the 
        // following holds.
        //
        // Consider extending each side of the triangle into a
        // line that bisects the x-y plane.
        //
        // Then, for each side of the triangle, the point must
        // lie either on that line or on the same side of the line
        // as the corner of the triangle facing the given side.
        //
        // In algebraic terms, consider triangle A-B-C and point P.
        //
        // Let AB denote the vector from point A to point B and so on.
        //
        // Then each of the following is true if A-B-C contains P:
        //
        // sign(AC cross AB) == sign(AP cross AB) OR (AP cross AB == 0)
        //
        // AND CYCLIC PERMUTATIONS IN THE SEQUENCE A->B->C->A->B->C ...
        //
        // So:
        bool triangle_has_point = true;

        // Containment condition for side AB
        const Vec2D vec_ac(vertices[0], vertices[2]);
        const Vec2D vec_ab(vertices[0], vertices[1]);
        const Vec2D vec_ap(vertices[0], point);

        const double ap_cross_ab = VectorMath2D::cross_product(vec_ap, vec_ab);
        const double ac_cross_ab = VectorMath2D::cross_product(vec_ac, vec_ab);

        const bool ab_check_passed = (ap_cross_ab == 0.0) ||
                                     (std::signbit(ac_cross_ab) == std::signbit(ap_cross_ab));

        triangle_has_point = triangle_has_point && ab_check_passed;

        // first cyclic permutation: Containment condition for side BC
        if (triangle_has_point)
        {
            const Vec2D vec_ba(vertices[1], vertices[0]);
            const Vec2D vec_bc(vertices[1], vertices[2]);
            const Vec2D vec_bp(vertices[1], point);

            const double bp_cross_bc = VectorMath2D::cross_product(vec_bp, vec_bc);
            const double ba_cross_bc = VectorMath2D::cross_product(vec_ba, vec_bc);

            const bool bc_check_passed = (bp_cross_bc == 0.0) ||
                                         (std::signbit(ba_cross_bc) == std::signbit(bp_cross_bc));

            triangle_has_point = triangle_has_point && bc_check_passed;
        }

        // second cyclic permutation: Containment condition for side CA
        if (triangle_has_point)
        {
            const Vec2D vec_cb(vertices[2], vertices[1]);
            const Vec2D vec_ca(vertices[2], vertices[0]);
            const Vec2D vec_cp(vertices[2], point);

            const double cp_cross_ca = VectorMath2D::cross_product(vec_cp, vec_ca);
            const double cb_cross_ca = VectorMath2D::cross_product(vec_cb, vec_ca);

            const bool ca_check_passed = (cp_cross_ca == 0.0) ||
                                         (std::signbit(cb_cross_ca) == std::signbit(cp_cross_ca));

            triangle_has_point = triangle_has_point && ca_check_passed;
        }

        return triangle_has_point;
    }
};


struct Polygon {
    std::vector<Coordinates> vertices;
    int orientation;

    Polygon(std::vector<Coordinates>& polygon_points)
        : vertices(polygon_points),
          orientation(0)
    {}

    int get_num_vertices() { return vertices.size(); }
};

namespace PolygonMath
{
    int compute_polygon_orientation(const Polygon& polygon)
    {
        // Computes whether the (simple) polygon is left-handed or right-handed.
        // Left-handed: points are ordered clockwise, return -1
        // Right-handed: points are ordered counter-clockwise, return +1
        //
        // To do this, we iterate through the polygon vertices, adding up the
        // subtended angle relative to each edge. Then we check its sign.
        //
        // Mathematically, this is a stencil operation on points
        // A=p[i], B=p[i+1], and C=p[i+2]. 
        // Consider the vectors v1=AB and v2=BC.
        // For a valid polygon, we know A, B, and C are all unique points.
        // We wish to know the angle subtended by vector AB if it is rotated
        // to orient along vector BC, where we require the angle in (-pi, pi).
        //
        // This is computed by the VectorMath2D::angle_between_vectors function.
        //
        // Next all we need do is add up the total angle subtended around the
        // entire polygon. Its sign is the orientation of the polygon.
        int num_vertices = polygon.get_num_vertices();
        double total_polygon_angle = 0.0;
        for (int i = 0; i < num_vertices; i++)
        {
            // We need the following two polygon vertex indices
            // so shift indices to the start of the polygon if we
            // are near the starting point.
            int ip1 = i+1 % num_vertices;
            int ip2 = i+2 % num_vertices;

            Coordinates v1(polygon.vertices[ip], polygon.vertices[ip1]);
            Coordinates v2(polygon.vertices[ip1], polygon.vertices[ip2]);
            double theta = VectorMath2D::angle_between_vectors(v1, v2);
            total_polygon_angle += theta;
        } 

        // Return the signed orientation of the polygon
        if (total_polygon_angle > 0.0) {
            return 1;
        } else if (total_polygon_angle < 0.0) {
            return -1;
        } else {
            return 0;
        }
    }
};


void triangulate(Polygon& polygon, std::vector<Triangle>& triangulation)
{
    // Clear the triangulation buffer
    triangulation.resize(0);

    // Create a list of vertex indices for the polygon points
    const int num_vertices = polygon.get_num_vertices();

    std::vector<int> vertex_ids;
    vertex_ids.resize(num_vertices);

    for (int i = 0; i < num_vertices; i++)
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
        area += TriangleMath::compute_triangle_area(triangle);
    }

    return area;
}


int main()
{
    // Read CSV file into a Polygon
    //
    // Compute the Polygon orientation
    //
    // Compute the Polygon triangulation
    //
    // Compute the triangulation total area
    //
    // Output total area, triangulation
    return 0;
}
