#ifndef MAPGENERATOR_H
#define MAPGENERATOR_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "extents2d.h"
#include "dcel.h"
#include "poissondiscsampler.h"
#include "voronoi.h"

#include "jsoncons/json.hpp"

namespace gen {

class MapGenerator {

public:
    MapGenerator();
    MapGenerator(Extents2d extents, double resolution);

    void initialize();
    void outputVoronoiDiagram(std::string filename);

private:

    Extents2d _extents;
    double _resolution;

    dcel::DCEL _voronoi;
    bool _isInitialized = false;

};

}

#endif