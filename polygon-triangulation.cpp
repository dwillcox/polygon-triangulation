#include <iostream>
#include <iomanip>
#include <fstream>
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


std::vector<Coordinates> read_polygon_from_csv(std::string csv_filename)
{
    std::cout << "Got csv filename: '" << csv_filename << "'\n";

    std::vector<Coordinates> xy_pairs;

    // read x,y pairs from the csv file
    // assuming they are stored one pair per line
    std::ifstream csv_fstream(csv_filename, std::ios::in);

    if (!csv_fstream.is_open()) {
        std::cout << "Could not open csv file for reading.\n";
    } else {
        std::cout << "Opened csv file for reading.\n";
        std::string line = "";
        const std::string delimiter = ",";
        while (std::getline(csv_fstream, line)) {
            std::cout << "read line '" << line << "' from file.\n";
            size_t delimiter_index = line.find(delimiter);
            // if this line contains a comma, interpret as x,y pair
            // otherwise do nothing
            if (delimiter_index != line.npos) {
                std::string xstr = line.substr(0, delimiter_index);
                std::string ystr = line.substr(delimiter_index+1, line.npos);
                double x = std::stod(xstr);
                double y = std::stod(ystr);
                Coordinates xy(x,y);
                xy_pairs.push_back(xy);
            }
        }
        csv_fstream.close();

        std::cout << "Read " << xy_pairs.size() << " x,y pairs:\n";
        for (const auto& parr : xy_pairs) {
            std::cout << parr.x << ", " << parr.y << std::endl;
        }
    }

    // If the first and last points are identical, then delete
    // the last point, as we store only the unique points that
    // define the polygon.
    if (xy_pairs.front().x == xy_pairs.back().x &&
        xy_pairs.front().y == xy_pairs.back().y) {
        xy_pairs.pop_back();
    }

    return xy_pairs;
}


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

        // project the angle phi into [-pi, pi) using the dot product sign.
        double theta = phi;

        if (v1_dot_v2 < 0.0) {
            if (phi > 0.0) {
                theta = std::numbers::pi - phi;
            } else if (phi < 0.0) {
                theta = -phi - std::numbers::pi;
            } else {
                theta = -std::numbers::pi;
            }
        }

        return theta;
    }
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


namespace TriangleMath
{
    double compute_triangle_area(const Triangle& triangle)
    {
        // Compute area of this triangle using the determinant
        // and the coordinates of its vertices.
        const auto& vertices = triangle.vertices;
        double area = 0.5 * (vertices[0].x * (vertices[1].y - vertices[2].y) +
                             vertices[1].x * (vertices[2].y - vertices[0].y) +
                             vertices[2].x * (vertices[0].y - vertices[1].y));

        // Depending on the coordinates, the determinant may be negative,
        // so we specify the absolute value of the determinant as the area.
        area = std::abs(area);

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
        // In the triangle vertices, A-B-C are indexed 0-1-2.
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
        const auto& vertices = triangle.vertices;

        // Containment condition for side AB
        const Coordinates vec_ac(vertices[0], vertices[2]);
        const Coordinates vec_ab(vertices[0], vertices[1]);
        const Coordinates vec_ap(vertices[0], point);

        const double ap_cross_ab = VectorMath2D::cross_product(vec_ap, vec_ab);
        const double ac_cross_ab = VectorMath2D::cross_product(vec_ac, vec_ab);

        const bool ab_check_passed = (ap_cross_ab == 0.0) ||
                                     (std::signbit(ac_cross_ab) == std::signbit(ap_cross_ab));

        triangle_has_point = triangle_has_point && ab_check_passed;

        // Check a special case in which the triangle ABC is one-dimensionally flat.
        const bool triangle_is_flat = (ac_cross_ab == 0.0);
        if (triangle_is_flat) {
            // In this case ac_cross_ab == 0, so there are two cases to check:
            if (ap_cross_ab != 0.0) {
                // * if ap_cross_ab is nonzero the point P is outside the triangle
                triangle_has_point = false;
            } else {
                // * if ap_cross_ab == 0, then P may still be outside the triangle if it does
                //   not lie on any of the line segments AB, BC, or CA.
                //   - first compute the vectors PA, PB, and PC
                const Coordinates vec_pa(point, vertices[0]);
                const Coordinates vec_pb(point, vertices[1]);
                const Coordinates vec_pc(point, vertices[2]);
                const double mag_pa = VectorMath2D::magnitude(vec_pa);
                const double mag_pb = VectorMath2D::magnitude(vec_pb);
                const double mag_pc = VectorMath2D::magnitude(vec_pc);
                //   - then, compute the dot products PA [dot] PB, PB [dot] PC, PC [dot] PA.
                const double pa_dot_pb = VectorMath2D::dot_product(vec_pa, vec_pb);
                const double pb_dot_pc = VectorMath2D::dot_product(vec_pb, vec_pc);
                const double pc_dot_pa = VectorMath2D::dot_product(vec_pc, vec_pa);
                //   - P is outside the triangle iff none of |PA|, |PB|, or |PC| are zero
                //     and all three dot products in the previous step have the same sign.
                const bool nonzero_magnitudes = (mag_pa != 0.0 && mag_pb != 0.0 && mag_pc != 0.0);
                const bool same_sign_dot_prod = (std::signbit(pa_dot_pb) == std::signbit(pb_dot_pc) &&
                                                 std::signbit(pa_dot_pb) == std::signbit(pc_dot_pa));
                if (nonzero_magnitudes && same_sign_dot_prod) {
                    triangle_has_point = false;
                }
            }
        }

        // first cyclic permutation: Containment condition for side BC
        // (no need to check the special case of a flat triangle again)
        if (triangle_has_point)
        {
            const Coordinates vec_ba(vertices[1], vertices[0]);
            const Coordinates vec_bc(vertices[1], vertices[2]);
            const Coordinates vec_bp(vertices[1], point);

            const double bp_cross_bc = VectorMath2D::cross_product(vec_bp, vec_bc);
            const double ba_cross_bc = VectorMath2D::cross_product(vec_ba, vec_bc);

            const bool bc_check_passed = (bp_cross_bc == 0.0) ||
                                         (std::signbit(ba_cross_bc) == std::signbit(bp_cross_bc));

            triangle_has_point = triangle_has_point && bc_check_passed;
        }

        // second cyclic permutation: Containment condition for side CA
        // (no need to check the special case of a flat triangle again)
        if (triangle_has_point)
        {
            const Coordinates vec_cb(vertices[2], vertices[1]);
            const Coordinates vec_ca(vertices[2], vertices[0]);
            const Coordinates vec_cp(vertices[2], point);

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

    Polygon()
        : vertices({}),
          orientation(0)
    {}

    Polygon(std::vector<Coordinates>& polygon_points)
        : vertices(polygon_points),
          orientation(0)
    {}

    int get_num_vertices() const { return vertices.size(); }
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
        std::cout << "got num_vertices: " << num_vertices << "\n";
        double total_polygon_angle = 0.0;
        for (int i = 0; i < num_vertices; i++)
        {
            // We need the following two polygon vertex indices
            // so shift indices to the start of the polygon if we
            // are near the starting point.
            int ip1 = (i+1) % num_vertices;
            int ip2 = (i+2) % num_vertices;

            std::cout << "calculating theta from indices: " << i << " " << ip1 << " " << ip2 << "\n";
            Coordinates v1(polygon.vertices[i], polygon.vertices[ip1]);
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


void triangulate_polygon(Polygon& polygon, std::vector<Triangle>& triangulation)
{
    // Clear the triangulation buffer
    triangulation.resize(0);

    // Create a list of vertex indices for the polygon points
    int num_vertices = polygon.get_num_vertices();

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
        // Utilities
        int num_vertices_remaining = vertex_ids.size();
        std::cout << "iterating vertex_ids with N vertices: " << num_vertices_remaining << "\n";

        // lambda for converting indices in the working list of vertices
        // to ring-based indices of vertices in the Polygon
        auto ring_poly_index = [&, num_vertices_remaining] (int i) -> int {
            int iRing = (i + num_vertices_remaining) % num_vertices_remaining;
            return vertex_ids[iRing];
        };

        // Candidate Triangle for the Ear at this index
        Triangle candidate_triangle;

        // lambda for checking if this vertex index is an Ear
        auto vertex_is_ear = [&] (int i) -> bool {
            // Given index i in the working list of unclipped vertex indices,
            // check if vertex i is an Ear using the surrounding
            // unclipped vertices at i-1 (A), i (B), and i+1 (C).
            //
            // It must satisfy the following:
            // Ear Condition 1: AC must lie within the polygon.
            // Ear Condition 2: Triangle ABC must not contain any vertices.
            //
            //// Check Ear Condition 1: AC_within_polygon
            // Let D be the vertex at i-2.
            // We check the angle of AC relative to DA (angleAC)
            // and compare with the angle of AB relative to DA (angleAB)
            //
            int iA = ring_poly_index(i-1);
            int iB = ring_poly_index(i);
            int iC = ring_poly_index(i+1);
            int iD = ring_poly_index(i-2);
            Coordinates vecDA(polygon.vertices[iD], polygon.vertices[iA]);
            Coordinates vecAC(polygon.vertices[iA], polygon.vertices[iC]);
            Coordinates vecAB(polygon.vertices[iA], polygon.vertices[iB]);
            double angleAC = VectorMath2D::angle_between_vectors(vecDA, vecAC);
            double angleAB = VectorMath2D::angle_between_vectors(vecDA, vecAB);

            // Now, both angleAC and angleAB are given in (-pi,pi)
            // with respect to the vector from D to A.
            //
            // We know the polygon edge proceeds from D to A to B.
            //
            // So which side of the path DAB does AC lie on?
            //
            // If angleAC is in (-pi, angleAB) then it is within the polygon
            // only if the polygon is left-handed, or clockwise-oriented.
            //
            // If angleAC is in (angleAB, pi) then it is within the polygon
            // only if the polygon is right-handed, or counter-clockwise-oriented.
            //
            // We know for a simple polygon, it is an error if angleAC == angleAB,
            // as that would require C and A be collocated.
            bool AC_within_polygon = false;

            std::cout << "A: " << polygon.vertices[iA].x << "," << polygon.vertices[iA].y << "\n";
            std::cout << "B: " << polygon.vertices[iB].x << "," << polygon.vertices[iB].y << "\n";
            std::cout << "C: " << polygon.vertices[iC].x << "," << polygon.vertices[iC].y << "\n";
            std::cout << "D: " << polygon.vertices[iD].x << "," << polygon.vertices[iD].y << "\n";
            std::cout << "angleAB: " << angleAB << "\n";
            std::cout << "angleAC: " << angleAC << "\n";
            std::cout << "ori: " << polygon.orientation << "\n";

            if (iD == iC && num_vertices_remaining == 3) {
                AC_within_polygon = true;
            } else {
                if (polygon.orientation == 1) {
                    // condition for right-handed polygon
                    //
                    // note that if angleAB == angleAC then BC is parallel to AB
                    // and B may be clipped, though it contributes zero total area.
                    if ((angleAB <= angleAC) && (angleAC < std::numbers::pi)) {
                        AC_within_polygon = true;
                    }
                } else if (polygon.orientation == -1) {
                    // condition for left-handed polygon
                    //
                    // note that if angleAB == angleAC then BC is parallel to AB
                    // and B may be clipped, though it contributes zero total area.
                    if ((-std::numbers::pi < angleAC) && (angleAC <= angleAB)) {
                        AC_within_polygon = true;
                    }
                } else {
                    // ERROR: polygon.orientation must be 1 or -1
                }
            }

            //// Check Ear Condition 2: ABC_is_empty
            bool ABC_is_empty = true;
            if (AC_within_polygon)
            {
                Triangle triangleABC(polygon.vertices[iA],
                                     polygon.vertices[iB],
                                     polygon.vertices[iC]);
                candidate_triangle = triangleABC;
                for (int ipoly = 0; ipoly < polygon.vertices.size(); ipoly++)
                {
                    // skip checking the triangle vertices A,B,C themselves
                    if (ipoly == iA || ipoly == iB || ipoly == iC) continue;

                    // check that no polygon vertex lies within triangle ABC
                    if (TriangleMath::triangle_contains_coordinates(triangleABC,
                                                                    polygon.vertices[ipoly]))
                    {
                        std::cout << "ABC contains point at: " << polygon.vertices[ipoly].x
                                  << "," << polygon.vertices[ipoly].y << "\n";
                        ABC_is_empty = false;
                        break;
                    }
                }
            }

            // check if the vertex indicated by the working index i is an Ear
            std::cout << "AC_within_polygon = " << AC_within_polygon << "\n";
            std::cout << "ABC_is_empty = " << ABC_is_empty << "\n";
            return AC_within_polygon && ABC_is_empty;
        };
 
        // Find an Ear
        int iEar = -1;
        for (int i = 0; i < num_vertices_remaining; i++) {
            if (vertex_is_ear(i)) {
                std::cout << "Vertex " << i << " is an ear!\n";
                // Add to our list of Triangles
                triangulation.push_back(candidate_triangle);

                iEar = i;
                break;
            } else {
                std::cout << "Vertex " << i << " is NOT an ear!\n";
            }
        }

        // ERROR if iEar == -1 and Ear not found
        if (iEar == -1) {
            std::cout << "could not find an ear with num_vertices_remaining = " << num_vertices_remaining << "\n";
            break;
        }

        // Last before we loop:
        // Remove the Ear vertex from our working list of vertices
        std::cout << "found iEar = " << iEar << "\n";
        vertex_ids.erase(vertex_ids.begin() + iEar);
    }

    if (vertex_ids.size() != 2) {
        std::cout << "ERROR: incomplete triangulation. " << vertex_ids.size() << " vertices remaining.\n";
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


int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cout << "Please run this program as './polygon-triangulation.exe [CSV file]'\n";
        return -1;
    }

    std::string csv_filename = argv[1];
    std::cout << "Reading Polygon from CSV file: " << csv_filename << std::endl;

    // Read CSV file into a Polygon
    std::vector<Coordinates> polygon_coordinates = read_polygon_from_csv(csv_filename);
    std::cout << "hi1\n";
    Polygon polygon(polygon_coordinates);
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
    double total_area = compute_triangulated_area(triangulation);
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
