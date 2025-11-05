#include <PolygonTriangulation_Geometry_Shapes.h>

#include <iostream>
#include <fstream>
#include <string>

namespace PolygonTriangulation
{

    bool Polygon::read_from_csv(std::string csv_filename)
    {
#if defined(DEBUGPRINT)
        std::cout << "Got csv filename: '" << csv_filename << "'\n";
#endif

        // clear any existing vertices and reset existing orientation
        vertices.resize(0);
        orientation = 0;

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

#if defined(DEBUGPRINT)
                std::cout << "read line '" << line << "' from file.\n";
#endif

                size_t delimiter_index = line.find(delimiter);
                // if this line contains a comma, interpret as x,y pair
                // otherwise do nothing
                if (delimiter_index != line.npos) {
                    std::string xstr = line.substr(0, delimiter_index);
                    std::string ystr = line.substr(delimiter_index+1, line.npos);
                    double x = std::stod(xstr);
                    double y = std::stod(ystr);
                    Coordinates xy(x,y);
                    vertices.push_back(xy);
                }
            }
            csv_fstream.close();

#if defined(DEBUGPRINT)
            std::cout << "Read " << vertices.size() << " x,y pairs:\n";

            for (const auto& parr : vertices) {
                std::cout << parr.x << ", " << parr.y << std::endl;
            }
#endif
        }

        // If the first and last points are identical, then delete
        // the last point, as we store only the unique points that
        // define the polygon.
        if (vertices.front().x == vertices.back().x &&
            vertices.front().y == vertices.back().y) {
            vertices.pop_back();
        }

        if (vertices.size() > 0)
            return true;
        else
            return false;
    }

};

