#ifndef DELAUNAY_H
#define DELAUNAY_H

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <stdlib.h>

#include "dcel.h"
#include "geometry.h"

namespace Delaunay {

using namespace dcel;

DCEL triangulate(std::vector<Point> &points);

void _getSuperTriangle(std::vector<Point> &points,
                       Point *p1, Point *p2, Point *p3);
DCEL _initTriangulation(std::vector<Point> &points);
Face _locateTriangleAtPoint(Point &p, DCEL &T);
Point _computeTriangleCentroid(Face &f, DCEL &T);
bool _isSegmentIntersectingEdge(Point &p0, Point &p1, HalfEdge &h, DCEL &T);
bool _isPointInsideTriangle(Point &p, Face &f, DCEL &T);
double _pointToEdgeDistance(Point &p, HalfEdge &h, DCEL &T);
void _insertPointIntoTriangulation(Point p, Face f, DCEL &T);
void _insertPointIntoTriangle(Point p, Face f, DCEL &T);
void _insertPointIntoTriangleEdge(Point p, Face f, HalfEdge h, DCEL &T);
bool _isEdgeLegal(Vertex p0, HalfEdge e, DCEL &T);
void _legalizeEdge(Vertex p, HalfEdge h, DCEL &T);
void _cleanup(DCEL &T);
void _getCleanupInvalidFaces(DCEL &T, std::vector<Face> &invalidFaces);
void _getCleanupInvalidEdges(DCEL &T, std::vector<HalfEdge> &invalidEdges,
                                      std::vector<HalfEdge> &invalidTwins);
void _getCleanupUpdateVertices(DCEL &T, 
                               std::vector<HalfEdge> &invalidEdges,
                               std::vector<Vertex> &vertices);
void _updateCleanupVertices(DCEL &T, 
                            std::vector<Vertex> &updateVertices, 
                            std::vector<bool> &invalidEdgeTable, 
                            std::vector<bool> &invalidFaceTable);
void _removeInvalidCleanupComponents(DCEL &T,
                                     std::vector<HalfEdge> &invalidEdges,
                                     std::vector<Face> &invalidFaces);

}

#endif