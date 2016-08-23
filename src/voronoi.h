#ifndef VORONOI_H
#define VORONOI_H

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <stdlib.h>

#include "dcel.h"
#include "delaunay.h"
#include "geometry.h"

namespace Voronoi {

using namespace dcel;

DCEL voronoi(std::vector<Point> &points);
DCEL delaunayToVoronoi(DCEL &T);

typedef std::vector<std::vector<std::pair<int, int> > > VertexEdgeTable;

void _createVoronoiVertices(DCEL &T, DCEL &V, std::vector<int> &vertexToFaceTable);
Point _computeVoronoiVertex(DCEL &T, Face &f);
bool _isBoundaryVertex(DCEL &T, Vertex &v);
VertexEdgeTable _initVertexEdgeTable(
    DCEL &T, DCEL &V, std::vector<int> &delaunayFaceToVertexTable);
void _initVertexIncidentEdges(DCEL &V, VertexEdgeTable &vertexEdges);
void _getVoronoiCellEdgeLoop(Vertex &delaunayVertex,
                             DCEL &T, DCEL &V,
                             std::vector<int> &delaunayFaceToVertexTable,
                             VertexEdgeTable &vertexEdgeTable,
                             std::vector<HalfEdge> &edgeLoop);
void _initVoronoiFaceFromEdgeLoop(std::vector<HalfEdge> &edgeLoop,
                                  DCEL &V,
                                  VertexEdgeTable &vertexEdges);

}

#endif