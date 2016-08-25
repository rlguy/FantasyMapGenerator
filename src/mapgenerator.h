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
#include "vertexmap.h"

#include "jsoncons/json.hpp"

namespace gen {

class MapGenerator {

public:
    MapGenerator();
    MapGenerator(Extents2d extents, double resolution);

    void initialize();
    void outputVoronoiDiagram(std::string filename);
    void outputVertices(std::string filename);
    void outputEdgeVertices(std::string filename);
    void outputInteriorVertices(std::string filename);

private:
    jsoncons::json _getExtentsJSON();
    void _outputVertices(std::vector<dcel::Vertex> &verts, 
                         std::string filename);

    Extents2d _extents;
    double _resolution;

    dcel::DCEL _voronoi;
    VertexMap _vertexMap;
    bool _isInitialized = false;

    double _samplePadFactor = 1.5;

};

}

#endif