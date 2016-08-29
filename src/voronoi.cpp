#include "voronoi.h"

dcel::DCEL Voronoi::voronoi(std::vector<Point> &points) {
    if (points.size() == 0) {
        return DCEL();
    }

    DCEL V = Delaunay::triangulate(points);
    return delaunayToVoronoi(V);
}

dcel::DCEL Voronoi::delaunayToVoronoi(DCEL &T) {
    DCEL V;
    std::vector<int> voronoiVertexToFaceTable;
    _createVoronoiVertices(T, V, voronoiVertexToFaceTable);

    std::vector<int> delaunayFaceToVertexTable(T.faces.size(), -1);
    for (unsigned int i = 0; i < voronoiVertexToFaceTable.size(); i++) {
        delaunayFaceToVertexTable[voronoiVertexToFaceTable[i]] = i;
    }

    VertexEdgeTable vertexEdges;
    vertexEdges = _initVertexEdgeTable(T, V, delaunayFaceToVertexTable);
    _initVertexIncidentEdges(V, vertexEdges);

    std::vector<Face> incidentFaces;
    std::vector<HalfEdge> edgeLoop;
    for (unsigned int vidx = 0; vidx < T.vertices.size(); vidx++) {
        Vertex v = T.vertices[vidx];
        if (_isBoundaryVertex(T, v)) {
            continue;
        }

        edgeLoop.clear();
        _getVoronoiCellEdgeLoop(v, T, V, 
                                delaunayFaceToVertexTable, 
                                vertexEdges,
                                edgeLoop);
        _initVoronoiFaceFromEdgeLoop(edgeLoop, V, vertexEdges);
    }

    return V;
}

void Voronoi::_createVoronoiVertices(DCEL &T, 
                                    DCEL &V, std::vector<int> &vertexToFaceTable) {

    for (unsigned int i = 0; i < T.faces.size(); i++) {
        if (T.faces[i].outerComponent.ref == -1) {
            continue;
        }

        Point p = _computeVoronoiVertex(T, T.faces[i]);
        V.createVertex(p);
        vertexToFaceTable.push_back(i);
    }
}

dcel::Point Voronoi::_computeVoronoiVertex(DCEL &T, Face &f) {
    // A voronoi vertex is a center of a circle that passes through
    // points p0, pi, pj
    HalfEdge h = T.outerComponent(f);
    Point p0 = T.origin(h).position;
    Point pi = T.origin(T.next(h)).position;
    Point pj = T.origin(T.prev(h)).position;

    Point p(0.5*(pi.x + pj.x), 0.5*(pi.y + pj.y));
    Point r(-(pj.y - pi.y), pj.x - pi.x);

    Point q(0.5*(pi.x + p0.x), 0.5*(pi.y + p0.y));
    Point s(-(p0.y - pi.y), p0.x - pi.x);

    Point center;
    bool isIntersection = Geometry::lineIntersection(p, r, q, s, &center);
    if (!isIntersection) {
        return p0;
    }

    return center;
}

bool Voronoi::_isBoundaryVertex(DCEL &T, Vertex &v) {
    if (v.incidentEdge.ref == -1) {
        return true;
    }

    HalfEdge h = T.incidentEdge(v);
    Ref startid = h.id;

    do {
        if (h.incidentFace.ref == -1) {
            return true;
        }

        h = T.twin(h);
        if (h.next.ref == -1) {
            return true;
        }

        h = T.next(h);
    } while (h.id != startid);

    return false;
}

Voronoi::VertexEdgeTable Voronoi::_initVertexEdgeTable(
                                  DCEL &T, DCEL &V, 
                                  std::vector<int> &delaunayFaceToVertexTable) {
    VertexEdgeTable vertexEdges;
    vertexEdges.reserve(V.vertices.size());
    for (unsigned int i = 0; i < V.vertices.size(); i++) {
        vertexEdges.push_back(std::vector<std::pair<int, int> >());
        vertexEdges[i].reserve(3);
    }

    std::vector<Face> incidentFaces;
    for (unsigned int vidx = 0; vidx < T.vertices.size(); vidx++) {
        Vertex v = T.vertices[vidx];
        if (_isBoundaryVertex(T, v)) {
            continue;
        }

        incidentFaces.clear();
        T.getIncidentFaces(v, incidentFaces);

        for (unsigned int fidx = 0; fidx < incidentFaces.size(); fidx++) {
            Face fi = incidentFaces[fidx];
            Face fj = fidx == 0 ? incidentFaces.back() : incidentFaces[fidx - 1];

            int refi = delaunayFaceToVertexTable[fi.id.ref];
            int refj = delaunayFaceToVertexTable[fj.id.ref];
            Vertex vi = V.getVertex(refi);
            Vertex vj = V.getVertex(refj);

            HalfEdge eij = V.createHalfEdge();
            eij.origin = vi.id;
            V.updateHalfEdge(eij);

            std::pair<int, int> vertexEdge(vj.id.ref, eij.id.ref);
            vertexEdges[vi.id.ref].push_back(vertexEdge);
        }
    }

    return vertexEdges;
}

void Voronoi::_initVertexIncidentEdges(DCEL &V, VertexEdgeTable &vertexEdges) {
    for (unsigned int i = 0; i < vertexEdges.size(); i++) {
        Ref refeij(vertexEdges[i][0].second);
        HalfEdge eij = V.getHalfEdge(refeij);

        Vertex vi = V.getVertex(i);
        vi.incidentEdge = eij.id;
        V.updateVertex(vi);
    }
}

void Voronoi::_getVoronoiCellEdgeLoop(Vertex &delaunayVertex,
                                       DCEL &T, DCEL &V,
                                       std::vector<int> &delaunayFaceToVertexTable,
                                       VertexEdgeTable &vertexEdges,
                                       std::vector<HalfEdge> &edgeLoop) {

    std::vector<Face> incidentFaces;
    incidentFaces.reserve(6);
    T.getIncidentFaces(delaunayVertex, incidentFaces);

    for (unsigned int fidx = 0; fidx < incidentFaces.size(); fidx++) {
        Face fi = incidentFaces[fidx];
        Face fj = fidx == 0 ? incidentFaces.back() : incidentFaces[fidx - 1];

        int refi = delaunayFaceToVertexTable[fi.id.ref];
        int refj = delaunayFaceToVertexTable[fj.id.ref];
        Vertex vi = V.getVertex(refi);
        Vertex vj = V.getVertex(refj);

        HalfEdge eij;
        for (unsigned int eidx = 0; eidx < vertexEdges[vi.id.ref].size(); eidx++) {
            std::pair<int, int> vertexEdge = vertexEdges[vi.id.ref][eidx];
            if (vertexEdge.first == vj.id.ref) {
                Ref refeij(vertexEdge.second);
                eij = V.getHalfEdge(refeij);
            }
        }

        edgeLoop.push_back(eij);
    }
}

void Voronoi::_initVoronoiFaceFromEdgeLoop(std::vector<HalfEdge> &edgeLoop,
                                            DCEL &V,
                                            VertexEdgeTable &vertexEdges) {

    Face cellface = V.createFace();
    cellface.outerComponent = edgeLoop[0].id;
    V.updateFace(cellface);

    HalfEdge eij, ejk, eri, eji;
    Vertex vi, vj;
    for (unsigned hidx = 0; hidx < edgeLoop.size(); hidx++) {
        eij = edgeLoop[hidx];
        ejk = hidx == 0 ? edgeLoop.back() : edgeLoop[hidx - 1];
        eri = hidx == edgeLoop.size() - 1 ? edgeLoop[0] : edgeLoop[hidx + 1];
        vi = V.origin(eij);
        vj = V.origin(ejk);

        bool isHalfEdgeFound = false;
        for (unsigned int eidx = 0; eidx < vertexEdges[vj.id.ref].size(); eidx++) {
            std::pair<int, int> vertexEdge = vertexEdges[vj.id.ref][eidx];
            if (vertexEdge.first == vi.id.ref) {
                eji = V.getHalfEdge(vertexEdge.second);
                isHalfEdgeFound = true;
                break;
            }
        }

        if (!isHalfEdgeFound) {
            eji = V.createHalfEdge();
            eji.origin = vj.id;
            eji.twin = eij.id;
            V.updateHalfEdge(eji);

            std::pair<int, int> vertexEdge(vi.id.ref, eji.id.ref);
            vertexEdges[vj.id.ref].push_back(vertexEdge);
        }

        eij.origin = vi.id;
        eij.twin = eji.id;
        eij.incidentFace = cellface.id;
        eij.next = ejk.id;
        eij.prev = eri.id;
        V.updateHalfEdge(eij);
    }
}
