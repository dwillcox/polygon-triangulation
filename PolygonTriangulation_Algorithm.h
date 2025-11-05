#ifndef POLY_TRI_ALGORITHM_H_
#define POLY_TRI_ALGORITHM_H_

#include <PolygonTriangulation_Geometry_Shapes.h>

#include <vector>

namespace PolygonTriangulation
{

    void triangulate_polygon(Polygon& polygon, std::vector<Triangle>& triangulation);

};


#endif // POLY_TRI_ALGORITHM_H_
