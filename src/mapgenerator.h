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

    void addCity();

    void outputVoronoiDiagram(std::string filename);
    void outputHeightMap(std::string filename);
    void outputDrawData(std::string filename);

private:
    typedef std::vector<dcel::Vertex> VertexList;

    struct Segment {
        dcel::Point p1;
        dcel::Point p2;
    };

    struct City {
        dcel::Point position;
    };

    jsoncons::json _getExtentsJSON();
    void _outputVertices(std::vector<dcel::Vertex> &verts, 
                         std::string filename);
    std::vector<double> _computeFaceValues(NodeMap<double> &heightMap);
    std::vector<dcel::Point> _computeFacePositions();
    dcel::Point _computeFacePosition(int fidx);
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

    void _getContourDrawData(std::vector<std::vector<double> > &data);
    void _getContourPaths(std::vector<VertexList> &paths);
    void _getLandFaces(std::vector<bool> &isLandFace);
    void _getFaceHeights(std::vector<double> &faceHeights);
    bool _isLand(double isolevel);
    void _cleanupLandFaces(std::vector<bool> &isLandFace);
    void _getConnectedFaces(int seed, std::vector<bool> &isLandFace,
                            std::vector<bool> &isFaceProcessed,
                            std::vector<int> &faces);
    bool _isContourEdge(dcel::HalfEdge &h, std::vector<bool> &isLandFace);
    bool _isContourEdge(dcel::Vertex &v1, dcel::Vertex &v2, 
                        std::vector<bool> &isLandFace);
    void _getContourPath(int seed, std::vector<bool> &isContourVertex, 
                                   std::vector<bool> &isEndVertex, 
                                   std::vector<bool> &isLandFace,
                                   std::vector<bool> &isVertexInContour, 
                                   VertexList &path);

    void _getRiverDrawData(std::vector<std::vector<double> > &data);
    void _getRiverPaths(std::vector<VertexList> &paths);
    void _getRiverVertices(VertexList &vertices);
    bool _isLandVertex(int vidx, std::vector<bool> &isLandFace);
    bool _isCoastVertex(int vidx, std::vector<bool> &isLandFace);
    void _getFixedRiverVertices(VertexList &riverVertices, 
                                VertexList &fixedVertices);
    VertexList _smoothPath(VertexList &path,
                           double factor);

    void _getSlopeDrawData(std::vector<double> &data);
    void _getSlopeSegments(std::vector<Segment> &segments);
    void _calculateHorizontalSlopeMap(NodeMap<double> &slopeMap);
    double _calculateHorizontalSlope(int i);
    void _calculateVerticalSlopeMap(NodeMap<double> &slopeMap);
    double _calculateVerticalSlope(int i);
    void _calculateVertexNormal(int vidx, double *nx, double *ny, double *nz);

    void _getCityDrawData(std::vector<double> &data);
    dcel::Point _getCityLocation();
    void _getCityScores(NodeMap<double> &cityScores);
    double _getPointDistance(dcel::Point &p1, dcel::Point &p2);
    double _pointToEdgeDistance(dcel::Point p);

    Extents2d _extents;
    double _resolution;

    dcel::DCEL _voronoi;
    VertexMap _vertexMap;
    NodeMap<double> _heightMap;
    NodeMap<double> _fluxMap;
    NodeMap<int> _flowMap;
    bool _isInitialized = false;

    double _samplePadFactor = 2.5;
    double _fluxCapPercentile = 0.995;
    double _maxErosionRate = 50.0;
    double _erosionRiverFactor = 500.0;
    double _erosionCreepFactor = 500.0;
    double _defaultErodeAmount = 0.1;
    double _riverFluxThreshold = 0.05;
    double _riverSmoothingFactor = 0.5;
    double _isolevel = 0.0;
    double _minIslandFaceThreshold = 35;

    double _minSlopeThreshold = 0.09;
    double _minSlope = 0.0;
    double _maxSlope = 0.7;
    double _minSlopeAngle = 0.2;
    double _maxSlopeAngle = 1.5;
    double _minSlopeLength = 0.75;
    double _maxSlopeLength = 1.0;
    double _minVerticalSlope = -0.25;
    double _maxVerticalSlope = 0.05;

    double _fluxScoreBonus = 2.0;
    double _nearEdgeScorePenalty = 0.1;
    double _nearCityScorePenalty = 2.0;
    double _maxPenaltyDistance = 4.0;

    std::vector<City> _cities;
};

}

#endif