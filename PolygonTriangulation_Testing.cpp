#include <PolygonTriangulation_Testing.h>

#include <PolygonTriangulation_Geometry_Shapes.h>
#include <PolygonTriangulation_Geometry_Math.h>

#include <iostream>
#include <numbers>
#include <cassert>
#include <cmath>

namespace Testing
{
    double relative_error(double a, double b)
    {
        if (a == b) {
            return 0.0;
        } else {
            return std::abs(a-b)/std::max(std::abs(a), std::abs(b));
        }
    }

    namespace Shapes
    {
        int test_coordinate_initialization()
        {
            // Test Coordinate initialization

            {
                // test 1: default
                Coordinates c;
                assert(c.x == 0.0 && c.y == 0.0 && "SUCCESS: default initialization");
            }

            {
                // test 2: x and y coordinates provided
                Coordinates c(1.0, -7.0);
                assert(c.x == 1.0 && c.y == -7.0 && "SUCCESS: coordinate pair specified");
            }

            {
                // test 3: initialize relative coordinates from the first
                //         to the second pair of provided coordinates
                Coordinates c_from(3,3);
                Coordinates c_to(6,0);
                Coordinates c(c_from, c_to);
                assert(c.x == 3.0 && c.y == -3.0 && "SUCCESS: relative coordinate initialization");
            }
                
            return 0;
        }


        int test_triangle_initialization()
        {
            // test Triangle initialization

            {
                // default init
                Triangle t;
                const auto& vertices = t.vertices;
                assert(vertices[0].x == 0.0 && vertices[0].y == 0.0 &&
                       vertices[1].x == 0.0 && vertices[1].y == 0.0 &&
                       vertices[2].x == 0.0 && vertices[2].y == 0.0 &&
                       "SUCCESS: triangle default init");
            }

            {
                // init vertices with specified coordinates
                Coordinates a(1,1);
                Coordinates b(-5.0, 4.0);
                Coordinates c(33.0, 6);
                Triangle t(a,b,c);
                const auto& vertices = t.vertices;
                assert(vertices[0].x == 1.0 && vertices[0].y == 1.0 &&
                       vertices[1].x == -5.0 && vertices[1].y == 4.0 &&
                       vertices[2].x == 33.0 && vertices[2].y == 6.0 &&
                       "SUCCESS: triangle vertex init");
            }

            return 0;
        }
    };

    namespace Vectors
    {
        int test_dot_product()
        {
            Coordinates u(1.0, 0.0);
            Coordinates v(1.0, 1.0);
            Coordinates w(-1.0, -1.0);

            double u_dot_v = VectorMath::dot_product(u, v);
            assert(u_dot_v == 1.0 && "SUCCESS: acute dot product");

            double u_dot_w = VectorMath::dot_product(u, w);
            assert(u_dot_w == -1.0 && "SUCCESS: obtuse dot product");

            return 0;
        }


        int test_cross_product()
        {
            Coordinates u(1.0, 0.0);
            Coordinates v(1.0, 1.0);
            Coordinates w(-1.0, -1.0);
            Coordinates q(0.0, 1.0);

            double u_cross_v = VectorMath::cross_product(u, v);
            assert(u_cross_v == 1.0 && "SUCCESS: positive cross product");

            double u_cross_w = VectorMath::cross_product(u, w);
            assert(u_cross_w == -1.0 && "SUCCESS: negative cross product");

            double u_cross_q = VectorMath::cross_product(u, q);
            assert(u_cross_q == 1.0 && "SUCCESS: unit cross product");

            double v_cross_w = VectorMath::cross_product(v, w);
            assert(v_cross_w == 0.0 && "SUCCESS: zero cross product");

            return 0;
        }


        int test_magnitude()
        {
            Coordinates u(1.0, 0.0);
            Coordinates v(1.0, 1.0);
            Coordinates w(-1.0, -1.0);
            Coordinates q(0.0, 1.0);

            double mag_u = VectorMath::magnitude(u);
            assert(mag_u == 1.0);

            double mag_v = VectorMath::magnitude(v);
            assert(mag_v == std::sqrt(2.0));

            double mag_w = VectorMath::magnitude(w);
            assert(mag_w == std::sqrt(2.0));

            double mag_q = VectorMath::magnitude(q);
            assert(mag_q == 1.0);

            return 0;
        }


        int test_angle_between_vectors()
        {
            Coordinates u(1.0, 0.0);
            Coordinates v(1.0, 1.0);
            Coordinates w(-1.0, -1.0);
            Coordinates q(0.0, 1.0);

            double angle_uv = VectorMath::angle_between_vectors(u, v);
            assert(relative_error(angle_uv, std::numbers::pi/4.0) <= 1.0e-14);

            double angle_uw = VectorMath::angle_between_vectors(u, w);
            assert(relative_error(angle_uw, -3.0*std::numbers::pi/4.0) <= 1.0e-14);

            double angle_uq = VectorMath::angle_between_vectors(u, q);
            assert(relative_error(angle_uq, std::numbers::pi/2.0) <= 1.0e-14);

            return 0;
        }
    };


    namespace Triangles
    {
        int test_compute_triangle_area()
        {
            {
                // test an easy triangle
                Coordinates a(0,0);
                Coordinates b(2,0);
                Coordinates c(2,1);
                Triangle triangle(a,b,c);

                double area = TriangleMath::compute_triangle_area(triangle);
                assert(relative_error(area, 1.0) <= 1.0e-14);
            }

            {
                // test a flat triangle
                Coordinates a(0,0);
                Coordinates b(1,0);
                Coordinates c(2,0);
                Triangle triangle(a,b,c);

                double area = TriangleMath::compute_triangle_area(triangle);
                assert(relative_error(area, 0.0) <= 1.0e-14);
            }

            return 0;
        }

        int test_triangle_contains_coordinates()
        {
             {
                // test an easy triangle
                Coordinates a(0,0);
                Coordinates b(2,0);
                Coordinates c(2,1);
                Triangle triangle(a,b,c);

                // try a point outside
                {
                    Coordinates point(3,3);
                    bool point_inside = TriangleMath::triangle_contains_coordinates(triangle, point);
                    assert(!point_inside);
                }

                // try a point inside
                {
                    Coordinates point(1.5,0.1);
                    bool point_inside = TriangleMath::triangle_contains_coordinates(triangle, point);
                    assert(point_inside);
                }

                // try a point on an edge
                {
                    Coordinates point(1.5,0.0);
                    bool point_inside = TriangleMath::triangle_contains_coordinates(triangle, point);
                    assert(point_inside);
                }
            }

            {
                // test a flat triangle
                Coordinates a(0,0);
                Coordinates b(1,0);
                Coordinates c(2,0);
                Triangle triangle(a,b,c);

                // try a point outside
                {
                    Coordinates point(3,3);
                    bool point_inside = TriangleMath::triangle_contains_coordinates(triangle, point);
                    assert(!point_inside);
                }

                // try a point inside
                {
                    Coordinates point(1.5,0.0);
                    bool point_inside = TriangleMath::triangle_contains_coordinates(triangle, point);
                    assert(point_inside);
                }
            }

            return 0;
        }
    };

    namespace Polygons
    {
        int test_compute_polygon_orientation()
        {
            {
                // test a counter-clockwise polygon
                std::vector<Coordinates> vertices;

                // vertices of a counter-clockwise square
                vertices.emplace_back(0,0);
                vertices.emplace_back(1,0);
                vertices.emplace_back(1,1);
                vertices.emplace_back(0,1);

                Polygon polygon(vertices);

                int orientation = PolygonMath::compute_polygon_orientation(polygon);
                assert(orientation == 1);
            }

            {
                // test a clockwise polygon
                std::vector<Coordinates> vertices;

                // vertices of a clockwise square
                vertices.emplace_back(0,0);
                vertices.emplace_back(0,1);
                vertices.emplace_back(1,1);
                vertices.emplace_back(1,0);

                Polygon polygon(vertices);

                int orientation = PolygonMath::compute_polygon_orientation(polygon);
                assert(orientation == -1);
            }

            return 0;
        }
    };
};


