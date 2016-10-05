#ifndef MAPGENERATOR_H
#define MAPGENERATOR_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
#include <string>

#include "extents2d.h"
#include "dcel.h"
#include "poissondiscsampler.h"
#include "voronoi.h"
#include "vertexmap.h"
#include "nodemap.h"
#include "fontface.h"
#include "spatialpointgrid.h"

#include "jsoncons/json.hpp"

namespace gen {

class MapGenerator {

public:
    MapGenerator();
    MapGenerator(Extents2d extents, double resolution, 
                 int imgwidth, int imgheight);
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
    void addTown();

    void outputVoronoiDiagram(std::string filename);
    void outputHeightMap(std::string filename);
    std::vector<char> getDrawData();

private:
    typedef std::vector<dcel::Vertex> VertexList;

    struct Segment {
        dcel::Point p1;
        dcel::Point p2;
    };

    struct CityLocation {
        dcel::Point position;
        int faceid;
    };

    struct City {
        dcel::Point position;
        int faceid;
        std::vector<double> movementCosts;
    };

    struct Town {
        dcel::Point position;
        int faceid;
    };

    struct CollisionData {
        int id;
    };

    struct LabelCandidate {
        std::string text;
        std::string fontface;
        int fontsize;
        dcel::Point position;
        Extents2d extents;
        std::vector<Extents2d> charextents;
        int cityid;

        double orientationScore;
        double edgeScore;
        double markerScore;
        double contourScore;
        double riverScore;
        double borderScore;
        double baseScore;

        int parentIdx;
        int collisionIdx;
        std::vector<CollisionData> collisionData;
    };

    struct Label {
        std::string text;
        std::string fontface;
        int fontsize;
        dcel::Point position;

        std::vector<LabelCandidate> candidates;
        int candidateIdx = -1;
        double score = 0.0;
    };

    struct LabelOffset {
        dcel::Point offset;
        double score = 0.0;

        LabelOffset() {}
        LabelOffset(dcel::Point p, double s) : offset(p), score(s) {}
    };

    void _initializeVoronoiData();
    void _initializeMapData();
    void _initializeNeighbourMap();
    void _initializeFaceNeighbours();
    void _initializeFaceVertices();
    void _initializeFaceEdges();
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
    bool _isLandFace(int fidx);
    void _getLandFaces(std::vector<bool> &isLandFace);
    void _getFaceHeights(std::vector<double> &faceHeights);
    bool _isLand(double isolevel);
    void _initializeLandFaceTable();
    void _cleanupLandFaces(std::vector<bool> &isLandFace);
    void _getConnectedFaces(int seed, std::vector<bool> &isLandFace,
                            std::vector<bool> &isFaceProcessed,
                            std::vector<int> &faces);
    bool _isContourEdge(dcel::HalfEdge &h);
    bool _isContourEdge(dcel::Vertex &v1, dcel::Vertex &v2);
    void _getContourPath(int seed, std::vector<bool> &isContourVertex, 
                                   std::vector<bool> &isEndVertex,
                                   std::vector<bool> &isVertexInContour, 
                                   VertexList &path);

    void _getRiverDrawData(std::vector<std::vector<double> > &data);
    void _getRiverPaths(std::vector<VertexList> &paths);
    void _getRiverVertices(VertexList &vertices);
    bool _isLandVertex(int vidx);
    bool _isCoastVertex(int vidx);
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
    void _getTownDrawData(std::vector<double> &data);
    CityLocation _getCityLocation();
    void _getCityScores(NodeMap<double> &cityScores);
    double _getPointDistance(dcel::Point &p1, dcel::Point &p2);
    double _pointToEdgeDistance(dcel::Point p);
    void _updateCityMovementCost(City &city);

    void _getTerritoryDrawData(std::vector<std::vector<double> > &data);
    void _getTerritoryBorders(std::vector<VertexList> &borders);
    void _getFaceTerritories(std::vector<int> &faceTerritories);
    void _cleanupFaceTerritories(std::vector<int> &faceTerritories);
    void _smoothTerritoryBoundaries(std::vector<int> &faceTerritories);
    void _getConnectedTerritories(std::vector<int> &faceTerritories,
                                  std::vector<std::vector<int> > &connected);
    void _getConnectedTerritory(int fidx, 
                                std::vector<int> &faceTerritories,
                                std::vector<bool> &isFaceProcessed,
                                std::vector<int> &connectedFaces);
    void _getDisjointTerritories(std::vector<int> &faceTerritories,
                                 std::vector<std::vector<int> > &connected, 
                                 std::vector<std::vector<int> > &disjoint);
    void _claimDisjointTerritories(std::vector<std::vector<int> > &disjoint,
                                   std::vector<int> &faceTerritories);
    int _getTerritoryOwner(std::vector<int> &territory,
                            std::vector<int> &faceTerritories);
    void _getBorderPaths(std::vector<int> &faceTerritories, 
                         std::vector<VertexList> &borders);
    void _getBorderEdges(std::vector<int> &faceTerritories, 
                         std::vector<dcel::HalfEdge> &borderEdges);
    bool _isBorderEdge(dcel::HalfEdge &h, std::vector<int> &faceTerritories);
    bool _isBorderEdge(dcel::Vertex &v1, dcel::Vertex &v2, 
                       std::vector<int> &faceTerritories);
    void _getBorderPath(int vidx, std::vector<int> &faceTerritories, 
                                  std::vector<bool> &isEndVertex, 
                                  std::vector<bool> &isVertexProcessed,
                                  VertexList &path);

    void _getLabelDrawData(std::vector<jsoncons::json> &data);
    void _initializeLabels(std::vector<Label> &labels);
    void _initializeMarkerLabels(std::vector<std::string> names, 
                                 std::vector<Label> &labels);
    void _initializeAreaLabels(std::vector<std::string> names, 
                               std::vector<Label> &labels);
    void _initializeCityLabel(City &city, std::string &name, Label &label);
    void _initializeTownLabel(Town &town, std::string &name, Label &label);
    void _initializeAreaLabel(City &city, std::string &name, Label &label);
    std::vector<std::string> _getLabelNames(int num);
    std::vector<LabelCandidate> _getMarkerLabelCandidates(Label label, 
                                                          double markerRadius);
    std::vector<LabelCandidate> _getAreaLabelCandidates(Label label,
                                                        City &city);
    void _getAreaLabelSamples(City &city, std::vector<dcel::Point> &samples);
    void _shuffleVector(std::vector<int> &vector);
    dcel::Point _getPixelCoordinates(dcel::Point &p);
    dcel::Point _getMapCoordinates(dcel::Point &p);
    Extents2d _getTextExtents(std::string text, dcel::Point pos);
    std::vector<Extents2d> _getCharacterExtents(std::string text, dcel::Point pos);
    jsoncons::json _getLabelJSON(LabelCandidate &label);
    dcel::Point _normalizeMapCoordinate(dcel::Point &p);
    dcel::Point _normalizeMapCoordinate(double x, double y);
    std::vector<LabelOffset> _getLabelOffsets(Label label, double radius);
    void _initializeMarkerLabelScores(std::vector<Label> &labels);
    void _initializeAreaLabelScores(std::vector<Label> &labels);
    void _initializeLabelEdgeScores(std::vector<Label> &labels);
    static bool _sortAreaLabelsByScore(LabelCandidate label1, LabelCandidate label2);
    double _getEdgeScore(Extents2d extents);
    void _initializeLabelMarkerScores(std::vector<Label> &labels);
    double _computeLabelMarkerScore(Extents2d extents);
    void _initializeAreaLabelMarkerScores(std::vector<Label> &labels);
    double _computeAreaLabelMarkerScore(Extents2d extents);
    void _initializeLabelContourScores(std::vector<Label> &labels);
    void _getDataPoints(std::vector<std::vector<double> > &data,
                        std::vector<dcel::Point> &points);
    void _computeContourScores(Label &label, SpatialPointGrid &grid);
    int _getLabelPointCount(LabelCandidate &c, SpatialPointGrid &grid);
    void _initializeLabelRiverScores(std::vector<Label> &labels);
    void _computeRiverScores(Label &label, SpatialPointGrid &grid);
    void _initializeLabelBorderScores(std::vector<Label> &labels);
    void _computeBorderScores(Label &label, SpatialPointGrid &grid);
    void _initializeAreaLabelOrientationScores(std::vector<Label> &labels);
    void _initializeAreaLabelOrientationScore(Label &label);
    double _calculationAreaLabelOrientationScore(LabelCandidate &label,
                                               SpatialPointGrid &territoryGrid,
                                               SpatialPointGrid &enemyGrid,
                                               SpatialPointGrid &waterGrid);
    void _initializeLabelBaseScores(std::vector<Label> &labels);
    double _computeLabelBaseScore(LabelCandidate &label);
    void _generateLabelPlacements(std::vector<Label> &labels);
    void _randomizeLabelPlacements(std::vector<Label> &labels);
    int _randomRangeInt(int minval, int maxval);
    double _randomRangeDouble(double minval, double maxval);
    void _initializeLabelCollisionData(std::vector<Label> &labels);
    void _initializeLabelCollisionData(std::vector<Label> &labels, 
                                       LabelCandidate &label);
    double _calculateLabelPlacementScore(std::vector<Label> &labels);
    double _calculateLabelPlacementScore(Label &label,
                                         std::vector<bool> &isCandidateActive);
    bool _isLabelOverlapping(LabelCandidate &label1, LabelCandidate &label2);
    bool _isExtentsOverlapping(Extents2d &e1, Extents2d &e2);

    Extents2d _extents;
    double _resolution;
    int _imgwidth;
    int _imgheight;
    int _defaultImageHeight = 1080;
    double _defaultResolution = 0.1;
    double _defaultExtentsWidth = 20.0*1.7777;
    double _defaultExtentsHeight = 20.0;

    dcel::DCEL _voronoi;
    VertexMap _vertexMap;
    NodeMap<std::vector<int> > _neighbourMap;
    std::vector<std::vector<int> > _faceNeighbours;
    std::vector<std::vector<int> > _faceVertices;
    std::vector<std::vector<int> > _faceEdges;
    NodeMap<double> _heightMap;
    NodeMap<double> _fluxMap;
    NodeMap<int> _flowMap;
    bool _isInitialized = false;

    std::vector<bool> _isLandFaceTable;
    bool _isLandFaceTableInitialized = false;

    double _samplePadFactor = 3.5;
    int _poissonSamplerKValue = 25;
    double _fluxCapPercentile = 0.995;
    double _maxErosionRate = 50.0;
    double _erosionRiverFactor = 500.0;
    double _erosionCreepFactor = 500.0;
    double _defaultErodeAmount = 0.1;
    double _riverFluxThreshold = 0.06;
    double _riverSmoothingFactor = 0.5;
    double _isolevel = 0.0;
    double _minIslandFaceThreshold = 35;

    double _minSlopeThreshold = 0.07;
    double _minSlope = 0.0;
    double _maxSlope = 0.7;
    double _minSlopeAngle = 0.2;
    double _maxSlopeAngle = 1.5;
    double _minSlopeLength = 0.75;
    double _maxSlopeLength = 1.3;
    double _minVerticalSlope = -0.25;
    double _maxVerticalSlope = 0.05;

    double _fluxScoreBonus = 2.0;
    double _nearEdgeScorePenalty = 0.5;
    double _nearCityScorePenalty = 2.0;
    double _nearTownScorePenalty = 1.5;
    double _maxPenaltyDistance = 4.0;

    double _landDistanceCost = 0.2;
    double _seaDistanceCost = 0.4;
    double _uphillCost = 0.1;
    double _downhillCost = 1.0;
    double _fluxCost = 0.8;
    double _landTransitionCost = 0.0;

    int _numTerritoryBorderSmoothingInterations = 3;
    double _territoryBorderSmoothingFactor = 0.5;

    std::vector<std::vector<double> > _contourData;
    std::vector<std::vector<double> > _riverData;
    std::vector<std::vector<double> > _borderData;
    std::vector<int> _territoryData;

    std::vector<City> _cities;
    std::vector<Town> _towns;

    FontFace _fontData;
    double _cityMarkerRadius = 10.0;    // in pixels
    double _townMarkerRadius = 5.0;
    std::string _cityLabelFontFace;
    std::string _townLabelFontFace;
    std::string _areaLabelFontFace;
    int _cityLabelFontSize = 35;
    int _townLabelFontSize = 25;
    int _areaLabelFontSize = 35;
    int _numAreaLabelSamples = 500;
    int _numAreaLabelCandidates = 120;
    double _spatialGridResolutionFactor = 5.0;
    double _labelMarkerRadiusFactor = 1.0;
    double _areaLabelMarkerRadiusFactor = 7.5;
    double _edgeScorePenalty = 4.0;
    double _markerScorePenalty = 6.0;
    double _minContourScorePenalty = 0.5;
    double _maxContourScorePenalty = 1.5;
    double _minRiverScorePenalty = 0.7;
    double _maxRiverScorePenalty = 2.0;
    double _minBorderScorePenalty = 0.8;
    double _maxBorderScorePenalty = 2.0;
    double _overlapScorePenalty = 4.0;
    double _territoryScore = 0.0;
    double _enemyScore = 6.0;
    double _waterScore = 0.2;

    double _initialTemperature = 0.91023922;     // 1.0 / log(3)
    double _annealingFactor = 0.9;
    double _maxTemperatureChanges = 100;
    int _successfulRepositioningFactor = 5;
    int _totalRepositioningFactor = 20;
};

}

#endif