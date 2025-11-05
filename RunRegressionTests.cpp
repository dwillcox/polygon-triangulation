#include <PolygonTriangulation_Geometry_Shapes.h>
#include <PolygonTriangulation_Geometry_Math.h>
#include <PolygonTriangulation_Algorithm.h>
#include <PolygonTriangulation_Testing.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cassert>

int main(int argc, char** argv)
{
    using namespace PolygonTriangulation;

    if (argc != 1) {
        std::cout << "Please run this program as './RunRegressionTests.exe\n";
        return -1;
    }

    std::string csv_filename = "test_square_10-10.csv";
    std::cout << "Reading Polygon from CSV file: " << csv_filename << std::endl;

    // Read CSV file into a Polygon
    Polygon polygon;
    polygon.read_from_csv(csv_filename);

    // Compute the Polygon orientation
    polygon.orientation = PolygonMath::compute_polygon_orientation(polygon);

    // Compute the Polygon triangulation
    std::vector<Triangle> triangulation;
    triangulate_polygon(polygon, triangulation);

    // Compute the triangulation total area
    double total_area = TriangulationMath::compute_triangulated_area(triangulation);

    // Check triangulation
    std::cout << "Checking Triangulation ... ";

    std::vector<Triangle> reference_triangulation;
    reference_triangulation.emplace_back(Coordinates(-5, 5), Coordinates(-5, -5), Coordinates(5, -5));
    reference_triangulation.emplace_back(Coordinates(-5, 5), Coordinates(5, -5), Coordinates(5, 5));
    assert(reference_triangulation.size() == triangulation.size());
    for (int i = 0; i < triangulation.size(); i++) {
        const auto& vertices = triangulation[i].vertices;
        const auto& ref_vert = reference_triangulation[i].vertices;
        assert(vertices[0].x == ref_vert[0].x && vertices[0].y == ref_vert[0].y);
        assert(vertices[1].x == ref_vert[1].x && vertices[1].y == ref_vert[1].y);
        assert(vertices[2].x == ref_vert[2].x && vertices[2].y == ref_vert[2].y);
    }
    std::cout << "SUCCESS!\n";

    // Check total area
    std::cout << "Checking Total Area ... ";
    assert(Testing::relative_error(total_area, 100.0) <= 1.0e-14);
    std::cout << "SUCCESS!\n";

    return 0;
}
