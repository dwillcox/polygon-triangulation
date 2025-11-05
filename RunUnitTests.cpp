#include <PolygonTriangulation_Testing.h>

#include <iostream>
#include <functional>
#include <cassert>

int main()
{
    using namespace PolygonTriangulation;

    auto test_and_report = [] (std::string label, std::function<bool()> test) {
        std::cout << "Test: " << label << " ... ";

        bool success = test();

        if (success == true) {
            std::cout << "SUCCESS!\n";
        } else {
            std::cout << "FAILED!\n";
        }

        assert(success == true);
    };

    test_and_report("Initialize Coordinates",
                    Testing::Shapes::test_coordinate_initialization);

    test_and_report("Initialize Triangle",
                    Testing::Shapes::test_triangle_initialization);

    test_and_report("Calculating Vector Dot Product",
                    Testing::Vectors::test_dot_product);

    test_and_report("Calculating Vector Cross Product",
                    Testing::Vectors::test_cross_product);

    test_and_report("Calculating Vector Magnitude",
                    Testing::Vectors::test_magnitude);

    test_and_report("Calculating Angle Between Vectors",
                    Testing::Vectors::test_angle_between_vectors);

    test_and_report("Calculating Whether Triangle Contains Coordinates",
                    Testing::Triangles::test_triangle_contains_coordinates);

    test_and_report("Calculating Triangle Area",
                    Testing::Triangles::test_compute_triangle_area);

    test_and_report("Calculating Polygon Orientation",
                    Testing::Polygons::test_compute_polygon_orientation);

    std::cout << "All Tests Pass!\n";

    return 0;
}
