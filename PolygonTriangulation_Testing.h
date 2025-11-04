#ifndef POLY_TRI_TESTING_H_
#define POLY_TRI_TESTING_H_

namespace Testing
{
    double relative_error(double a, double b);

    namespace Shapes
    {
        int test_coordinate_initialization();

        int test_triangle_initialization();
    };

    namespace Vectors
    {
        int test_dot_product();
        int test_cross_product();
        int test_magnitude();
        int test_angle_between_vectors();
    };

    namespace Triangles
    {
        int test_triangle_contains_coordinates();
        int test_compute_triangle_area();
    };

    namespace Polygons
    {
        int test_compute_polygon_orientation();
    };
};

#endif // POLY_TRI_TESTING_H_
