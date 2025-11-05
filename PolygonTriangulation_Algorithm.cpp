#include <PolygonTriangulation_Algorithm.h>

#include <PolygonTriangulation_Geometry_Shapes.h>
#include <PolygonTriangulation_Geometry_Math.h>

#include <iostream>
#include <vector>
#include <cassert>

namespace PolygonTriangulation
{

    void triangulate_polygon(Polygon& polygon, std::vector<Triangle>& triangulation)
    {
        // Clear the triangulation buffer
        triangulation.resize(0);

        // Create a list of vertex indices for the polygon points
        int num_vertices = polygon.get_num_vertices();

        if (num_vertices <= 2) {
            std::cout << "In order to triangulate the polygon, the polygon requres more than two unique vertices.\n";
            assert(num_vertices > 2);
        }

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

#if defined(DEBUGPRINT)
            std::cout << "iterating vertex_ids with N vertices: " << num_vertices_remaining << "\n";
#endif

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
                double angleAC = VectorMath::angle_between_vectors(vecDA, vecAC);
                double angleAB = VectorMath::angle_between_vectors(vecDA, vecAB);

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

#if defined(DEBUGPRINT)
                std::cout << "A: " << polygon.vertices[iA].x << "," << polygon.vertices[iA].y << "\n";
                std::cout << "B: " << polygon.vertices[iB].x << "," << polygon.vertices[iB].y << "\n";
                std::cout << "C: " << polygon.vertices[iC].x << "," << polygon.vertices[iC].y << "\n";
                std::cout << "D: " << polygon.vertices[iD].x << "," << polygon.vertices[iD].y << "\n";
                std::cout << "angleAB: " << angleAB << "\n";
                std::cout << "angleAC: " << angleAC << "\n";
                std::cout << "ori: " << polygon.orientation << "\n";
#endif

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
#if defined(DEBUGPRINT)
                            std::cout << "ABC contains point at: " << polygon.vertices[ipoly].x
                                      << "," << polygon.vertices[ipoly].y << "\n";
#endif
                            ABC_is_empty = false;
                            break;
                        }
                    }
                }

                // check if the vertex indicated by the working index i is an Ear

#if defined(DEBUGPRINT)
                std::cout << "AC_within_polygon = " << AC_within_polygon << "\n";
                std::cout << "ABC_is_empty = " << ABC_is_empty << "\n";
#endif

                return AC_within_polygon && ABC_is_empty;
            };
     
            // Find an Ear
            int iEar = -1;
            for (int i = 0; i < num_vertices_remaining; i++) {
                if (vertex_is_ear(i)) {
#if defined(DEBUGPRINT)
                    std::cout << "Vertex " << i << " is an ear!\n";
#endif

                    // Add to our list of Triangles
                    triangulation.push_back(candidate_triangle);

                    iEar = i;
                    break;
                } else {
#if defined(DEBUGPRINT)
                    std::cout << "Vertex " << i << " is NOT an ear!\n";
#endif
                }
            }

            // ERROR if iEar == -1 and Ear not found
            if (iEar == -1) {
#if defined(DEBUGPRINT)
                std::cout << "could not find an ear with num_vertices_remaining = " << num_vertices_remaining << "\n";
#endif
                assert(iEar != -1);
                break;
            }

            // Last before we loop:
            // Remove the Ear vertex from our working list of vertices
#if defined(DEBUGPRINT)
            std::cout << "found iEar = " << iEar << "\n";
#endif
            vertex_ids.erase(vertex_ids.begin() + iEar);
        }

        if (vertex_ids.size() != 2) {
            std::cout << "ERROR: incomplete triangulation. " << vertex_ids.size() << " vertices remaining.\n";
            assert(vertex_ids.size() == 2);
        }

        // We are finished with triangulation. The list of triangles is now given
        // by the triangulation Triangle buffer.
    }

};

