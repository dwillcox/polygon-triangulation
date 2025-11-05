#include <PolygonTriangulation_Geometry_Math.h>

#include <PolygonTriangulation_Geometry_Shapes.h>

#include <iostream>
#include <numbers>
#include <cmath>

namespace PolygonTriangulation
{

    namespace VectorMath {
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

            const double ap_cross_ab = VectorMath::cross_product(vec_ap, vec_ab);
            const double ac_cross_ab = VectorMath::cross_product(vec_ac, vec_ab);

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
                    const double mag_pa = VectorMath::magnitude(vec_pa);
                    const double mag_pb = VectorMath::magnitude(vec_pb);
                    const double mag_pc = VectorMath::magnitude(vec_pc);
                    //   - then, compute the dot products PA [dot] PB, PB [dot] PC, PC [dot] PA.
                    const double pa_dot_pb = VectorMath::dot_product(vec_pa, vec_pb);
                    const double pb_dot_pc = VectorMath::dot_product(vec_pb, vec_pc);
                    const double pc_dot_pa = VectorMath::dot_product(vec_pc, vec_pa);
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

                const double bp_cross_bc = VectorMath::cross_product(vec_bp, vec_bc);
                const double ba_cross_bc = VectorMath::cross_product(vec_ba, vec_bc);

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

                const double cp_cross_ca = VectorMath::cross_product(vec_cp, vec_ca);
                const double cb_cross_ca = VectorMath::cross_product(vec_cb, vec_ca);

                const bool ca_check_passed = (cp_cross_ca == 0.0) ||
                                             (std::signbit(cb_cross_ca) == std::signbit(cp_cross_ca));

                triangle_has_point = triangle_has_point && ca_check_passed;
            }

            return triangle_has_point;
        }
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
            // This is computed by the VectorMath::angle_between_vectors function.
            //
            // Next all we need do is add up the total angle subtended around the
            // entire polygon. Its sign is the orientation of the polygon.
            int num_vertices = polygon.get_num_vertices();

#if defined(DEBUGPRINT)
            std::cout << "got num_vertices: " << num_vertices << "\n";
#endif

            double total_polygon_angle = 0.0;
            for (int i = 0; i < num_vertices; i++)
            {
                // We need the following two polygon vertex indices
                // so shift indices to the start of the polygon if we
                // are near the starting point.
                int ip1 = (i+1) % num_vertices;
                int ip2 = (i+2) % num_vertices;

#if defined(DEBUGPRINT)
                std::cout << "calculating theta from indices: " << i << " " << ip1 << " " << ip2 << "\n";
#endif

                Coordinates v1(polygon.vertices[i], polygon.vertices[ip1]);
                Coordinates v2(polygon.vertices[ip1], polygon.vertices[ip2]);
                double theta = VectorMath::angle_between_vectors(v1, v2);
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


    namespace TriangulationMath
    {
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
    };

};

