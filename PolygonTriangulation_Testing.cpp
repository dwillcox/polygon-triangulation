#include <PolygonTriangulation_Testing.h>

#include <PolygonTriangulation_Geometry_Shapes.h>

#include <cassert>

namespace Testing
{
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
};


