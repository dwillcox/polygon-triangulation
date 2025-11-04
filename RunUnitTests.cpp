#include <PolygonTriangulation_Testing.h>

#include <iostream>

int main()
{
    Testing::Shapes::test_coordinate_initialization();

    Testing::Shapes::test_triangle_initialization();

    Testing::Vectors::test_dot_product();

    Testing::Vectors::test_cross_product();

    Testing::Vectors::test_magnitude();

    Testing::Vectors::test_angle_between_vectors();

    Testing::Triangles::test_triangle_contains_coordinates();

    Testing::Triangles::test_compute_triangle_area();

    Testing::Polygons::test_compute_polygon_orientation();

    std::cout << "All Tests Pass!\n";

    return 0;
}
