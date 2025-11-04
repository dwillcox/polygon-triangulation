#ifndef POLY_TRI_GEOMETRY_MATH_H_
#define POLY_TRI_GEOMETRY_MATH_H_

#include <PolygonTriangulation_Geometry_Shapes.h>

namespace VectorMath {
    double dot_product(Coordinates v1, Coordinates v2);

    double cross_product(Coordinates v1, Coordinates v2);

    double magnitude(Coordinates v);

    double angle_between_vectors(Coordinates v1, Coordinates v2);
};


namespace TriangleMath
{
    double compute_triangle_area(const Triangle& triangle);

    bool triangle_contains_coordinates(const Triangle& triangle,
                                       const Coordinates& point);
};


namespace PolygonMath
{
    int compute_polygon_orientation(const Polygon& polygon);
};


namespace TriangulationMath
{
    double compute_triangulated_area(std::vector<Triangle>& triangulation);
};


#endif // POLY_TRI_GEOMETRY_MATH_H_
