#ifndef POLY_TRI_TESTING_H_
#define POLY_TRI_TESTING_H_

namespace PolygonTriangulation
{

    namespace Testing
    {
        double relative_error(double a, double b);

        namespace Shapes
        {
            bool test_coordinate_initialization();

            bool test_triangle_initialization();
        };

        namespace Vectors
        {
            bool test_dot_product();
            bool test_cross_product();
            bool test_magnitude();
            bool test_angle_between_vectors();
        };

        namespace Triangles
        {
            bool test_triangle_contains_coordinates();
            bool test_compute_triangle_area();
        };

        namespace Polygons
        {
            bool test_compute_polygon_orientation();
        };
    };

};

#endif // POLY_TRI_TESTING_H_
