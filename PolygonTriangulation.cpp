#include <PolygonTriangulation_Geometry_Shapes.h>
#include <PolygonTriangulation_Geometry_Math.h>
#include <PolygonTriangulation_Algorithm.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

/* Polygon Triangulation Overview
 *
 * We get a polygon specified in CSV format as x, y coordinate pairs, one per
 * line. We read these coordinates into a list of 2D points in a Polygon class.
 *
 * We are told the polygon is simple, but not whether the vertices are given in
 * clockwise order or in counter-clockwise order, so we next determine the
 * orientation and store that in the Polygon object.
 *
 * For the core triangulation algorithm, we use the Ear Clipping Method.
 *
 * The Ear Clipping Method identifies polygon vertices which qualify as "ears"
 * and removes them from the polygon by "clipping" them.
 *
 * A vertex is an "ear" if it appears as the "B" vertex in a triangle formed
 * by consecutive polygon vertices A-B-C such that the following hold:
 *
 * Ear Condition 1: Triangle (Ear) A-B-C contains no other polygon vertices.
 *
 * Ear Condition 2: Triangle (Ear) A-B-C lies within the polygon.
 *
 *                  Since A-B and B-C edges are polygon edges as well,
 *                  this is true iff line segment A-C is within the polygon.
 *
 * When clipping an "ear" we record the triangle formed by that ear in a list.
 *
 * We then identify and clip the next ear until the polygon is fully
 * triangulated.
 *
 * We obtain a list of triangles we identified during triangulation. Each
 * Triangle object stores its vertex coordinates for easy calculation.
 *
 * Finally, we report the triangulation and total area to the user.
 */

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cout << "Please run this program as './PolygonTriangulation.exe [CSV file]'\n";
        return -1;
    }

    std::string csv_filename = argv[1];
    std::cout << "Reading Polygon from CSV file: " << csv_filename << std::endl;

    // Read CSV file into a Polygon
    Polygon polygon;
    std::cout << "hi1\n";
    polygon.read_from_csv(csv_filename);
    std::cout << "hi2\n";

    // Compute the Polygon orientation
    polygon.orientation = PolygonMath::compute_polygon_orientation(polygon);
    std::cout << "hi3\n";

    // Compute the Polygon triangulation
    std::vector<Triangle> triangulation;
    std::cout << "hi4\n";
    triangulate_polygon(polygon, triangulation);
    std::cout << "hi5\n";

    // Compute the triangulation total area
    double total_area = TriangulationMath::compute_triangulated_area(triangulation);
    std::cout << "hi6\n";

    // Output triangulation
    std::cout << "Triangulation: list of triangle vertices as (x,y), (x,y), (x,y)\n";
    for (const auto& triangle : triangulation) {
        std::cout << "(" << triangle.vertices[0].x << "," << triangle.vertices[0].y << ")";
        std::cout << ", ";
        std::cout << "(" << triangle.vertices[1].x << "," << triangle.vertices[1].y << ")";
        std::cout << ", ";
        std::cout << "(" << triangle.vertices[2].x << "," << triangle.vertices[2].y << ")";
        std::cout << std::endl;
    }

    // Output total area
    std::cout << "Total Area: " << std::setprecision(20) << total_area << std::endl;

    return 0;
}
