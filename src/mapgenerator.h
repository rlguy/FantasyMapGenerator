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
#include "nodemap.h"

#include "jsoncons/json.hpp"

namespace gen {

class MapGenerator {

public:
    MapGenerator();
    MapGenerator(Extents2d extents, double resolution);

    void initialize();
    void normalize();
    void round();
    void relax();
    void setSeaLevel(double level);
    void setSeaLevelToMedian();
    void addHill(double px, double py, double radius, double height);
    void addCone(double px, double py, double radius, double height);
    void addSlope(double px, double py, double dirx, double diry, 
                  double radius, double height);
    void erode(double amount);
    void erode();

    void outputVoronoiDiagram(std::string filename);
    void outputVertices(std::string filename);
    void outputEdgeVertices(std::string filename);
    void outputInteriorVertices(std::string filename);
    void outputHeightMap(std::string filename);
    void outputContour(std::string filename, double isolevel);

private:
    jsoncons::json _getExtentsJSON();
    void _outputVertices(std::vector<dcel::Vertex> &verts, 
                         std::string filename);
    std::vector<double> _computeFaceHeights(NodeMap<double> &heightMap);
    std::vector<double> _computeContour(double isolevel);
    bool _isEdgeInMap(dcel::HalfEdge &h);
    bool _isContourEdge(dcel::HalfEdge &h, 
                        std::vector<double> &faceheights, 
                        double isolevel);
    void _calculateErosionMap(NodeMap<double> &erosionMap);
    void _fillDepressions();
    void _calculateFlowMap(NodeMap<int> &flowMap);
    void _calculateFluxMap(NodeMap<double> &fluxMap);
    double _calculateFluxCap(NodeMap<double> &fluxMap);
    void _calculateSlopeMap(NodeMap<double> &slopeMap);
    double _calculateSlope(int i);

    Extents2d _extents;
    double _resolution;

    dcel::DCEL _voronoi;
    VertexMap _vertexMap;
    NodeMap<double> _heightMap;
    NodeMap<double> _fluxMap;
    bool _isInitialized = false;

    double _samplePadFactor = 2.5;
    double _fluxCapPercentile = 0.995;
    double _maxErosionRate = 50.0;
    double _erosionRiverFactor = 500.0;
    double _erosionCreepFactor = 100.0;
    double _defaultErodeAmount = 0.1;

};

}

#endif