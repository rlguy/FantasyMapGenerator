#include "delaunay.h"

dcel::DCEL Delaunay::triangulate(std::vector<dcel::Point> &points) {
    if (points.size() == 0) {
        return DCEL();
    }

    DCEL T = _initTriangulation(points);
    while (!points.empty()) {
        Point p = points.back();
        points.pop_back();

        Face f = _locateTriangleAtPoint(p, T);
        if (f.id.ref != -1) {
            _insertPointIntoTriangulation(p, f, T);
        }
    }

    _cleanup(T);

    return T;
}

void Delaunay::_getSuperTriangle(std::vector<dcel::Point> &points,
                                 dcel::Point *p1, dcel::Point *p2, dcel::Point *p3) {
    double eps = 1e-3;
    double minx = points[0].x;
    double miny = points[0].y;
    double maxx = minx + eps;
    double maxy = miny + eps;

    for (unsigned int i = 0; i < points.size(); i++) {
        Point p = points[i];
        if (p.x < minx) { minx = p.x; }
        if (p.y < miny) { miny = p.y; }
        if (p.x > maxx) { maxx = p.x; }
        if (p.y > maxy) { maxy = p.y; }
    }

    double expand = fmax(0.1 * (maxx - minx), 0.1 * (maxy - miny));
    minx -= expand;
    miny -= 5*expand;
    maxx += expand;
    maxy += expand;

    p1->x = 0.5 * (minx + maxx);
    p1->y = maxy + 0.5*(maxy - miny);

    double m = (maxy - p1->y) / (maxx - p1->x);
    p2->x = (1.0 / m) * (miny - p1->y + m * p1->x);
    p2->y = miny;

    m = (maxy - p1->y) / (minx - p1->x);
    p3->x = (1.0 / m) * (miny - p1->y + m * p1->x);
    p3->y = miny;
}

dcel::DCEL Delaunay::_initTriangulation(std::vector<dcel::Point> &points) {
    Point s1, s2, s3;
    _getSuperTriangle(points, &s1, &s2, &s3);

    DCEL T;
    Vertex p1 = T.createVertex(s1);
    Vertex p2 = T.createVertex(s2);
    Vertex p3 = T.createVertex(s3);
    HalfEdge e12 = T.createHalfEdge();
    HalfEdge e23 = T.createHalfEdge();
    HalfEdge e31 = T.createHalfEdge();
    HalfEdge e13 = T.createHalfEdge();
    HalfEdge e32 = T.createHalfEdge();
    HalfEdge e21 = T.createHalfEdge();
    Face f0 = T.createFace();

    p1.incidentEdge = e12.id;
    p2.incidentEdge = e23.id;
    p3.incidentEdge = e31.id;
    T.updateVertex(p1);
    T.updateVertex(p2);
    T.updateVertex(p3);

    e12.origin = p1.id;
    e12.twin = e21.id;
    e12.incidentFace = f0.id;
    e12.next = e23.id;
    e12.prev = e31.id;
    T.updateHalfEdge(e12);

    e23.origin = p2.id;
    e23.twin = e32.id;
    e23.incidentFace = f0.id;
    e23.next = e31.id;
    e23.prev = e12.id;
    T.updateHalfEdge(e23);

    e31.origin = p3.id;
    e31.twin = e13.id;
    e31.incidentFace = f0.id;
    e31.next = e12.id;
    e31.prev = e23.id;
    T.updateHalfEdge(e31);

    e13.origin = p1.id;
    e13.twin = e31.id;
    e13.next = e32.id;
    e13.prev = e21.id;
    T.updateHalfEdge(e13);

    e32.origin = p3.id;
    e32.twin = e23.id;
    e32.next = e21.id;
    e32.prev = e13.id;
    T.updateHalfEdge(e32);

    e21.origin = p2.id;
    e21.twin = e12.id;
    e21.next = e13.id;
    e21.prev = e32.id;
    T.updateHalfEdge(e21);

    f0.outerComponent = e12.id;
    T.updateFace(f0);

    return T;
}

// Choose a random triangle, walk toward p until the containing 
// triangle is found.
dcel::Face Delaunay::_locateTriangleAtPoint(dcel::Point &p, dcel::DCEL &T) {
    Ref randfidx(rand() % T.faces.size());
    Face f = T.getFace(randfidx);

    int count = 0;
    int maxcount = (int)(2.0*sqrt(T.faces.size()));

    int faceHistory[3] = {-1, -1, -1};
    HalfEdge h;
    Point p0;
    for (;;) {
        if (_isPointInsideTriangle(p, f, T)) {
            return f;
        }

        p0 = _computeTriangleCentroid(f, T);
        bool isNeighbourFound = false;
        h = T.outerComponent(f);
        for (int i = 0; i < 3; i++) {
            if (_isSegmentIntersectingEdge(p0, p, h, T)) {
                f = T.incidentFace(T.twin(h));
                isNeighbourFound = true;
            }
            h = T.next(h);
        }

        if (!isNeighbourFound) {
            break;
        }
    
        // Due to numerical error, a point may not report to be contained
        // in any triangle. In this case, the search will travel back and
        // forth between two faces. Tracking the last two face id's will 
        // allow cycles to be detected
        faceHistory[2] = faceHistory[1];
        faceHistory[1] = faceHistory[0];
        faceHistory[0] = f.id.ref;
        if (faceHistory[0] == faceHistory[2]) {
            break;
        }

        count++;
        if (count > maxcount) {
            break;
        }
    }

    return Face();
}

bool Delaunay::_isPointInsideTriangle(dcel::Point &p, dcel::Face &f, dcel::DCEL &T) {
    HalfEdge h = T.outerComponent(f);
    Point p0 = T.origin(h).position;
    h = T.next(h);
    Point p1 = T.origin(h).position;
    Point p2 = T.origin(T.next(h)).position;

    double area = 0.5*(-p1.y*p2.x + p0.y*(-p1.x + p2.x) + p0.x*(p1.y - p2.y) + p1.x*p2.y);
    double s = 1.0/(2.0*area)*(p0.y*p2.x - p0.x*p2.y + (p2.y - p0.y)*p.x + (p0.x - p2.x)*p.y);
    double t = 1.0/(2.0*area)*(p0.x*p1.y - p0.y*p1.x + (p0.y - p1.y)*p.x + (p1.x - p0.x)*p.y);

    return s >= 0 && t >= 0 && 1 - s - t >= 0;
}

dcel::Point Delaunay::_computeTriangleCentroid(dcel::Face &f, dcel::DCEL &T) {
    HalfEdge h = T.outerComponent(f);
    Point p0 = T.origin(h).position;
    Point p1 = T.origin(T.next(h)).position;
    Point p2 = T.origin(T.prev(h)).position;

    double frac = 1.0 / 3.0;
    return Point(frac*(p0.x + p1.x + p2.x), frac*(p0.y + p1.y + p2.y));
}

bool Delaunay::_isSegmentIntersectingEdge(dcel::Point &A, dcel::Point &B, 
                                          dcel::HalfEdge &h, dcel::DCEL &T) {
    Point C = T.origin(h).position;
    Point D = T.origin(T.twin(h)).position;
    return Geometry::lineSegmentIntersection(A, B, C, D);
}

double Delaunay::_pointToEdgeDistance(dcel::Point &p0, dcel::HalfEdge &h, dcel::DCEL &T) {
    Point p1 = T.origin(h).position;
    Point p2 = T.origin(T.twin(h)).position;

    double vx = p2.x - p1.x;
    double vy = p2.y - p1.y;
    double len = sqrt(vx*vx + vy*vy);

    return fabs((vx*(p1.y - p0.y) - (p1.x - p0.x)*vy) / len);
}

void Delaunay::_insertPointIntoTriangulation(dcel::Point p, dcel::Face f, dcel::DCEL &T) {
    double eps = 1e-9;
    int closeEdgeCount = 0;
    HalfEdge closeEdge;

    HalfEdge h = T.outerComponent(f);
    for (int i = 0; i < 3; i++) {
        double dist = _pointToEdgeDistance(p, h, T);
        if (dist < eps) {
            closeEdge = h;
            closeEdgeCount++;

            if (closeEdgeCount == 2) {
                // Point is too close to an existing vertex
                return;
            }
        }

        if (i < 2) { 
            h = T.next(h);
        }
    }

    if (closeEdgeCount == 0) {
        _insertPointIntoTriangle(p, f, T);
    } else {
        _insertPointIntoTriangleEdge(p, f, closeEdge, T);
    }
}

void Delaunay::_insertPointIntoTriangle(dcel::Point p, dcel::Face f, dcel::DCEL &T) {
    // existing components
    HalfEdge eij = T.outerComponent(f);
    HalfEdge ejk = T.next(eij);
    HalfEdge eki = T.next(ejk);

    Face f1 = f;

    Vertex pi = T.origin(eij);
    Vertex pj = T.origin(ejk);
    Vertex pk = T.origin(eki);

    // new components
    HalfEdge eri = T.createHalfEdge();
    HalfEdge eir = T.createHalfEdge();
    HalfEdge erj = T.createHalfEdge();
    HalfEdge ejr = T.createHalfEdge();
    HalfEdge erk = T.createHalfEdge();
    HalfEdge ekr = T.createHalfEdge();

    Face f2 = T.createFace();
    Face f3 = T.createFace();

    Vertex pr = T.createVertex(p);

    // update existing components
    eij.next = ejr.id;
    eij.prev = eri.id;
    T.updateHalfEdge(eij);

    ejk.incidentFace = f2.id;
    ejk.next = ekr.id;
    ejk.prev = erj.id;
    T.updateHalfEdge(ejk);

    eki.incidentFace = f3.id;
    eki.next = eir.id;
    eki.prev = erk.id;
    T.updateHalfEdge(eki);

    f1.outerComponent = eij.id;
    T.updateFace(f1);

    // initialize new components
    eri.origin = pr.id;
    eri.twin = eir.id;
    eri.incidentFace = f1.id;
    eri.next = eij.id;
    eri.prev = ejr.id;
    T.updateHalfEdge(eri);

    eir.origin = pi.id;
    eir.twin = eri.id;
    eir.incidentFace = f3.id;
    eir.next = erk.id;
    eir.prev = eki.id;
    T.updateHalfEdge(eir);

    erj.origin = pr.id;
    erj.twin = ejr.id;
    erj.incidentFace = f2.id;
    erj.next = ejk.id;
    erj.prev = ekr.id;
    T.updateHalfEdge(erj);

    ejr.origin = pj.id;
    ejr.twin = erj.id;
    ejr.incidentFace = f1.id;
    ejr.next = eri.id;
    ejr.prev = eij.id;
    T.updateHalfEdge(ejr);

    erk.origin =  pr.id;
    erk.twin = ekr.id;
    erk.incidentFace = f3.id;
    erk.next = eki.id;
    erk.prev = eir.id;
    T.updateHalfEdge(erk);

    ekr.origin = pk.id;
    ekr.twin = erk.id;
    ekr.incidentFace = f2.id;
    ekr.next = erj.id;
    ekr.prev = ejk.id;
    T.updateHalfEdge(ekr);

    f2.outerComponent = ejk.id;
    T.updateFace(f2);

    f3.outerComponent = eki.id;
    T.updateFace(f3);

    pr.incidentEdge = eri.id;
    T.updateVertex(pr);

    _legalizeEdge(pr, eij, T);
    _legalizeEdge(pr, ejk, T);
    _legalizeEdge(pr, eki, T);
}

void Delaunay::_insertPointIntoTriangleEdge(dcel::Point p, dcel::Face f, 
                                            dcel::HalfEdge h, dcel::DCEL &T) {
    // existing components
    HalfEdge eij = h;
    HalfEdge ejk = T.next(eij);
    HalfEdge eki = T.next(ejk);

    HalfEdge eji = T.twin(eij);
    HalfEdge eil = T.next(eji);
    HalfEdge elj = T.next(eil);

    Face f1 = T.incidentFace(eij);
    Face f2 = T.incidentFace(eji);

    Vertex pj = T.origin(ejk);
    Vertex pk = T.origin(eki);
    Vertex pl = T.origin(elj);

    // component replacements
    HalfEdge eir = eij;
    HalfEdge eri = eji;

    // new components
    HalfEdge erj = T.createHalfEdge();
    HalfEdge ejr = T.createHalfEdge();
    HalfEdge erk = T.createHalfEdge();
    HalfEdge ekr = T.createHalfEdge();
    HalfEdge erl = T.createHalfEdge();
    HalfEdge elr = T.createHalfEdge();

    Face f3 = T.createFace();
    Face f4 = T.createFace();

    Vertex pr = T.createVertex(p);

    // update existing components
    ejk.incidentFace = f4.id;
    ejk.next = ekr.id;
    ejk.prev = erj.id;
    T.updateHalfEdge(ejk);

    eki.next = eir.id;
    eki.prev = erk.id;
    T.updateHalfEdge(eki);

    eil.next = elr.id;
    eil.prev = eri.id;
    T.updateHalfEdge(eil);

    elj.incidentFace = f3.id;
    elj.next = ejr.id;
    elj.prev = erl.id;
    T.updateHalfEdge(elj);

    f1.outerComponent = eki.id;
    T.updateFace(f1);

    f2.outerComponent = eil.id;
    T.updateFace(f2);

    pj.incidentEdge = ejk.id;
    T.updateVertex(pj);

    // update replacement components
    eir.next = erk.id;
    eir.prev = eki.id;
    T.updateHalfEdge(eir);

    eri.origin = pr.id;
    eri.next = eil.id;
    eri.prev = elr.id;
    T.updateHalfEdge(eri);

    // initialize new components
    erj.origin = pr.id;
    erj.twin = ejr.id;
    erj.incidentFace = f4.id;
    erj.next = ejk.id;
    erj.prev = ekr.id;
    T.updateHalfEdge(erj);

    ejr.origin = pj.id;
    ejr.twin = erj.id;
    ejr.incidentFace = f3.id;
    ejr.next = erl.id;
    ejr.prev = elj.id;
    T.updateHalfEdge(ejr);

    erk.origin = pr.id;
    erk.twin = ekr.id;
    erk.incidentFace = f1.id;
    erk.next = eki.id;
    erk.prev = eir.id;
    T.updateHalfEdge(erk);

    ekr.origin = pk.id;
    ekr.twin = erk.id;
    ekr.incidentFace = f4.id;
    ekr.next = erj.id;
    ekr.prev = ejk.id;
    T.updateHalfEdge(ekr);

    erl.origin = pr.id;
    erl.twin = elr.id;
    erl.incidentFace = f3.id;
    erl.next = elj.id;
    erl.prev = ejr.id;
    T.updateHalfEdge(erl);

    elr.origin = pl.id;
    elr.twin = erl.id;
    elr.incidentFace = f2.id;
    elr.next = eri.id;
    elr.prev = eil.id;
    T.updateHalfEdge(elr);

    f3.outerComponent = elj.id;
    T.updateFace(f3);

    f4.outerComponent = ejk.id;
    T.updateFace(f4);

    pr.incidentEdge = eri.id;
    T.updateVertex(pr);

    _legalizeEdge(pr, eil, T);
    _legalizeEdge(pr, elj, T);
    _legalizeEdge(pr, ejk, T);
    _legalizeEdge(pr, eki, T);
}

bool Delaunay::_isEdgeLegal(dcel::Vertex v, dcel::HalfEdge e, dcel::DCEL &T) {
    if (T.isBoundary(T.twin(e))) {
        return true;
    }

    Point p0 = v.position;
    Point pi = T.origin(e).position;
    Point pj = T.origin(T.twin(e)).position;
    Point pk = T.origin(T.prev(T.twin(e))).position;

    /*
        Computational Geometry Algorithms and Applications
        by Mark de Berg, Otfried Cheong, Marc van Kreveld, Mark Overmars

        Theorem 9.4:

        Let edge pipj be incident to triangles p0pipj and pipjpk, and let
        C be the circle through p0, pi, and pj. The edge pipj is illegal if
        and only if the point pk lies in the interior of C.
    */

    Point p(0.5*(pi.x + pj.x), 0.5*(pi.y + pj.y));
    Point r(-(pj.y - pi.y), pj.x - pi.x);

    Point q(0.5*(pi.x + p0.x), 0.5*(pi.y + p0.y));
    Point s(-(p0.y - pi.y), p0.x - pi.x);

    Point center;
    bool isIntersection = Geometry::lineIntersection(p, r, q, s, &center);
    if (!isIntersection) {
        return false;
    }

    double dvx = p0.x - center.x;
    double dvy = p0.y - center.y;
    double crsq = dvx*dvx + dvy*dvy;

    double dkx = pk.x - center.x;
    double dky = pk.y - center.y;
    double distsq = dkx*dkx + dky*dky;
    
    return distsq >= crsq;
}

void Delaunay::_legalizeEdge(dcel::Vertex pr, dcel::HalfEdge eij, dcel::DCEL &T) {
    if (_isEdgeLegal(pr, eij, T)) {
        return;
    }

    HalfEdge ejr = T.next(eij);
    HalfEdge eri = T.next(ejr);
    HalfEdge eji = T.twin(eij);
    HalfEdge eik = T.next(eji);
    HalfEdge ekj = T.next(eik);

    Face f1 = T.incidentFace(eij);
    Face f2 = T.incidentFace(eji);

    Vertex pi = T.origin(eij);
    Vertex pj = T.origin(eji);
    Vertex pk = T.origin(ekj);

    // replacement components
    HalfEdge erk = eij;
    HalfEdge ekr = eji;

    // update existing components
    ejr.incidentFace = f2.id;
    ejr.next = erk.id;
    ejr.prev = ekj.id;
    T.updateHalfEdge(ejr);

    eri.next = eik.id;
    eri.prev = ekr.id;
    T.updateHalfEdge(eri);

    eik.incidentFace = f1.id;
    eik.next = ekr.id;
    eik.prev = eri.id;
    T.updateHalfEdge(eik);

    ekj.next = ejr.id;
    ekj.prev = erk.id;
    T.updateHalfEdge(ekj);

    f1.outerComponent = ekr.id;
    T.updateFace(f1);

    f2.outerComponent = erk.id;
    T.updateFace(f2);

    pi.incidentEdge = eik.id;
    T.updateVertex(pi);

    pj.incidentEdge = ejr.id;
    T.updateVertex(pj);

    pk.incidentEdge = ekr.id;
    T.updateVertex(pk);

    pr.incidentEdge = erk.id;
    T.updateVertex(pr);

    // update replacement components
    erk.origin = pr.id;
    erk.twin = ekr.id;
    erk.incidentFace = f2.id;
    erk.next = ekj.id;
    erk.prev = ejr.id;
    T.updateHalfEdge(erk);

    ekr.origin = pk.id;
    ekr.twin = erk.id;
    ekr.incidentFace = f1.id;
    ekr.next = eri.id;
    ekr.prev = eik.id;
    T.updateHalfEdge(ekr);

    _legalizeEdge(pr, eik, T);
    _legalizeEdge(pr, ekj, T);
}

void Delaunay::_getCleanupInvalidFaces(dcel::DCEL &T, std::vector<dcel::Face> &invalidFaces) {
    Vertex p1 = T.vertices[0];
    Vertex p2 = T.vertices[1];
    Vertex p3 = T.vertices[2];
    T.getIncidentFaces(p1, invalidFaces);
    T.getIncidentFaces(p2, invalidFaces);
    T.getIncidentFaces(p3, invalidFaces);
}

void Delaunay::_getCleanupInvalidEdges(dcel::DCEL &T, 
                                       std::vector<dcel::HalfEdge> &invalidEdges,
                                       std::vector<dcel::HalfEdge> &invalidTwins) {
    Vertex p1 = T.vertices[0];
    Vertex p2 = T.vertices[1];
    Vertex p3 = T.vertices[2];
    T.getIncidentEdges(p1, invalidEdges);
    T.getIncidentEdges(p2, invalidEdges);
    T.getIncidentEdges(p3, invalidEdges);

    invalidTwins.reserve(invalidEdges.size());
    for (unsigned int i = 0; i < invalidEdges.size(); i++) {
        invalidTwins.push_back(T.twin(invalidEdges[i]));
    }
}

void Delaunay::_getCleanupUpdateVertices(dcel::DCEL &T, 
                                         std::vector<dcel::HalfEdge> &invalidEdges,
                                         std::vector<dcel::Vertex> &vertices) {
    Vertex p1 = T.vertices[0];
    Vertex p2 = T.vertices[1];
    Vertex p3 = T.vertices[2];

    std::vector<bool> updateVertexTable(T.vertices.size(), false);
    for (unsigned int i = 0; i < invalidEdges.size(); i++) {
        Vertex v = T.origin(T.twin(invalidEdges[i]));
        if (v.id == p1.id || v.id == p2.id || v.id == p3.id) {
            continue;
        }
        updateVertexTable[v.id.ref] = true;
    }

    for (unsigned int i = 0; i < updateVertexTable.size(); i++) {
        if (updateVertexTable[i]) {
            vertices.push_back(T.vertices[i]);
        }
    }
}

void Delaunay::_updateCleanupVertices(dcel::DCEL &T, 
                                      std::vector<dcel::Vertex> &updateVertices,
                                      std::vector<bool> &invalidEdgeTable,
                                      std::vector<bool> &invalidFaceTable) {

    std::vector<Ref> incidentEdges;
    for (unsigned int i = 0; i < updateVertices.size(); i++) {
        Vertex v = updateVertices[i];

        incidentEdges.clear();
        T.getIncidentEdges(v, incidentEdges);

        HalfEdge h;
        for (unsigned int hidx = 0; hidx < incidentEdges.size(); hidx++) {
            h = T.getHalfEdge(incidentEdges[hidx]);
            int ref = h.id.ref;
            if (!invalidEdgeTable[ref] && invalidFaceTable[h.incidentFace.ref]) {
                h.incidentFace = Ref();
                T.updateHalfEdge(h);
            }
        }

        if (invalidEdgeTable[v.incidentEdge.ref]) {
            for (unsigned int hidx = 0; hidx < incidentEdges.size(); hidx++) {
                h = T.getHalfEdge(incidentEdges[hidx]);
                if (!invalidEdgeTable[h.id.ref]) {
                    v.incidentEdge = h.id;
                    T.updateVertex(v);
                }
            }
        }
    }
}

void Delaunay::_removeInvalidCleanupComponents(dcel::DCEL &T,
                                               std::vector<dcel::HalfEdge> &invalidEdges,
                                               std::vector<dcel::Face> &invalidFaces) {
    for (unsigned int i = 0; i < invalidEdges.size(); i++) {
        HalfEdge h = invalidEdges[i];
        HalfEdge twin = T.twin(h);

        Ref id = h.id;
        h = HalfEdge();
        h.id = id;
        T.updateHalfEdge(h);

        id = twin.id;
        twin = HalfEdge();
        twin.id = id;
        T.updateHalfEdge(twin);
    }

    for (unsigned int i = 0; i < invalidFaces.size(); i++) {
        Face f = invalidFaces[i];
        f.outerComponent = Ref();
        T.updateFace(f);
    }

    Vertex p1 = T.vertices[0];
    Vertex p2 = T.vertices[1];
    Vertex p3 = T.vertices[2];
    p1.incidentEdge = Ref();
    p2.incidentEdge = Ref();
    p3.incidentEdge = Ref();
    T.updateVertex(p1);
    T.updateVertex(p2);
    T.updateVertex(p3);
}

void Delaunay::_cleanup(dcel::DCEL &T) {

    std::vector<Face> invalidFaces;
    _getCleanupInvalidFaces(T, invalidFaces);
    std::vector<bool> invalidFaceTable(T.faces.size(), false);
    for (unsigned int i = 0; i < invalidFaces.size(); i++) {
        invalidFaceTable[invalidFaces[i].id.ref] = true;
    }

    std::vector<HalfEdge> invalidEdges;
    std::vector<HalfEdge> invalidTwinEdges;
    _getCleanupInvalidEdges(T, invalidEdges, invalidTwinEdges);

    std::vector<bool> invalidEdgeTable(T.edges.size(), false);
    for (unsigned int i = 0; i < invalidEdges.size(); i++) {
        invalidEdgeTable[invalidEdges[i].id.ref] = true;
    }
    for (unsigned int i = 0; i < invalidTwinEdges.size(); i++) {
        invalidEdgeTable[invalidTwinEdges[i].id.ref] = true;
    }

    std::vector<Vertex> updateVertices;
    _getCleanupUpdateVertices(T, invalidEdges, updateVertices);
    _updateCleanupVertices(T, updateVertices, invalidEdgeTable, invalidFaceTable);

    _removeInvalidCleanupComponents(T, invalidEdges, invalidFaces);
}
