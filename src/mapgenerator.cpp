#include "mapgenerator.h"

gen::MapGenerator::MapGenerator() {
    Extents2d defaultExtents(-20, -20, 20, 20);
    double defaultResolution = 0.4;

    _extents = defaultExtents;
    _resolution = defaultResolution;
}

gen::MapGenerator::MapGenerator(Extents2d extents, double resolution) :
                                _extents(extents), _resolution(resolution) {}

void gen::MapGenerator::initialize() {
    std::cout << "Generating samples within " << 
                 _extents.maxx - _extents.minx << " x " << 
                 _extents.maxy - _extents.miny << 
                 " area with radius " << _resolution << "." << std::endl;

    // Boundary vertices in the Voronoi diagram cannot be used as
    // grid nodes and will be excluded from the map grid. Extending the
    // bounds for the sample area will generate a Voronoi diagram that
    // extends past the map bounds, resulting in less excluded grid nodes.
    double pad = _samplePadFactor * _resolution;
    Extents2d sampleExtents(_extents.minx - pad, _extents.miny - pad,
                            _extents.maxx + pad, _extents.maxy + pad);

    std::vector<dcel::Point> samples;
    samples = PoissonDiscSampler::generateSamples(sampleExtents, 
                                                  _resolution, 25);

    std::cout << "Finished Generating " << samples.size() << 
                 " samples." << std::endl;
    std::cout << std::endl;

    std::cout << "Generating Voronoi diagram with " << samples.size() << 
                 " points." << std::endl;

    _voronoi = Voronoi::voronoi(samples);

    std::cout << "Finished generating Voronoi diagram." << std::endl;
    std::cout << "# Vertices:   " << _voronoi.vertices.size() << std::endl;
    std::cout << "# Half Edges: " << _voronoi.edges.size() << std::endl;
    std::cout << "# Faces:      " << _voronoi.faces.size() << std::endl;
    std::cout << std::endl;

    std::cout << "Initializing map vertices." << std::endl;
    _vertexMap = VertexMap(&_voronoi, _extents);
    std::cout << "Finished initializing map vertices." << std::endl;
    std::cout << "# Vertices: " << _vertexMap.vertices.size() << std::endl;
    std::cout << "# Interior: " << _vertexMap.interior.size() << std::endl;
    std::cout << "# Edge:     " << _vertexMap.edge.size() << std::endl;
    std::cout << std::endl;

    _heightMap = NodeMap<double>(&_vertexMap, 0);
    _isInitialized = true;
}

void gen::MapGenerator::normalize() {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }
    _heightMap.normalize();
}

void gen::MapGenerator::round() {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }
    _heightMap.round();
}

void gen::MapGenerator::relax() {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }
    _heightMap.relax();
}

void gen::MapGenerator::setSeaLevel(double level) {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }
    _heightMap.setLevel(level);
}

void gen::MapGenerator::setSeaLevelToMedian() {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }
    _heightMap.setLevelToMedian();
}

void gen::MapGenerator::addHill(double px, double py, double r, double height) {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }
    /*
        Create a rounded hill where height falls off smoothly
    */

    // coefficients for the Wyvill kernel function
    double coef1 = (4.0 / 9.0)*(1.0 / (r*r*r*r*r*r));
    double coef2 = (17.0 / 9.0)*(1.0 / (r*r*r*r));
    double coef3 = (22.0 / 9.0)*(1.0 / (r*r));
    double rsq = r*r;

    dcel::Point v;
    for (unsigned int i = 0; i < _heightMap.size(); i++) {
        v = _vertexMap.vertices[i].position;
        double dx = v.x - px;
        double dy = v.y - py;
        double dsq = dx*dx + dy*dy;
        if (dsq < rsq) {
            double kernel = 1.0 - coef1*dsq*dsq*dsq + coef2*dsq*dsq - coef3*dsq;
            double hval = _heightMap(i);
            _heightMap.set(i, hval + height*kernel);
        }
    }
}

void gen::MapGenerator::addCone(double px, double py, double radius, double height) {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }
    /*
        Create a cone where height falls off linearly
    */

    double invradius = 1.0 / radius;
    double rsq = radius * radius;
    dcel::Point v;
    for (unsigned int i = 0; i < _heightMap.size(); i++) {
        v = _vertexMap.vertices[i].position;
        double dx = v.x - px;
        double dy = v.y - py;
        double dsq = dx*dx + dy*dy;
        if (dsq < rsq) {
            double dist = sqrt(dsq);
            double kernel = 1.0 - dist * invradius;
            double hval = _heightMap(i);
            _heightMap.set(i, hval + height*kernel);
        }
    }
}

void gen::MapGenerator::addSlope(double px, double py, double dirx, double diry, 
                                 double radius, double height) {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }

    /*
        Add a slope that runs parallel to the line with direction (dirx, diry)
        and position (px, py). (dirx, diry) must be a unit vector. The slope
        varies in height linearly from 0 to height from the left side of the 
        direction vector to the right side of the direction vector.
    */

    dcel::Point v;
    for (unsigned int i = 0; i < _heightMap.size(); i++) {
        v = _vertexMap.vertices[i].position;
        
        double dx = px - v.x;
        double dy = py - v.y;
        double dot = dx*dirx + dy*diry;
        double distx = dx - dot*dirx;
        double disty = dy - dot*diry;
        double dist = sqrt(distx*distx + disty*disty);
        dist = fmin(dist, radius);

        double cross = dx*diry - dy*dirx;
        double min, max;
        if (cross < 0) {
            min = 0.5*height;
            max = 0;
        } else {
            min = 0.5*height;
            max = height;
        }

        double fieldval = min + (dist / radius)*(max - min);
        double hval = _heightMap(i);
        _heightMap.set(i, hval + fieldval);
    }
}

void gen::MapGenerator::erode()  {
    erode(_defaultErodeAmount);
}

void gen::MapGenerator::erode(double amount) {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }

    NodeMap<double> erosionMap(&_vertexMap, 0.0);
    _calculateErosionMap(erosionMap);

    for (unsigned int i = 0; i < _heightMap.size(); i++) {
        double currlevel = _heightMap(i);
        double newlevel = currlevel - amount * erosionMap(i);
        _heightMap.set(i, newlevel);
    }
}

void gen::MapGenerator::outputVoronoiDiagram(std::string filename) {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }
    using jsoncons::json;

    std::vector<std::vector<double> > faceVertices;
    faceVertices.reserve(_voronoi.faces.size());

    for (unsigned int i = 0; i < _voronoi.faces.size(); i++) {
        dcel::Face f = _voronoi.faces[i];

        if (f.outerComponent.ref == -1) {
            continue;
        }

        dcel::HalfEdge h = _voronoi.outerComponent(f);
        dcel::Ref startRef = h.id;

        std::vector<double> vertices;
        dcel::Vertex v;
        do {
            v = _voronoi.origin(h);
            vertices.push_back(v.position.x);
            vertices.push_back(v.position.y);
            h = _voronoi.next(h);
        } while (h.id != startRef);
        faceVertices.push_back(vertices);
    }
    
    json output;
    output["faces"] = faceVertices;
    output["extents"] = _getExtentsJSON();

    std::ofstream file(filename);
    file << output;
    file.close();
}

void gen::MapGenerator::outputVertices(std::string filename) {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }

    _outputVertices(_vertexMap.vertices, filename);
}

void gen::MapGenerator::outputEdgeVertices(std::string filename) {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }

    _outputVertices(_vertexMap.edge, filename);
}

void gen::MapGenerator::outputInteriorVertices(std::string filename) {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }

    _outputVertices(_vertexMap.interior, filename);
}

void gen::MapGenerator::outputHeightMap(std::string filename) {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }

    std::vector<double> facecolors = _computeFaceHeights(_heightMap);
    double min = facecolors[0];
    double max = facecolors[0];
    for (unsigned int i = 0; i < facecolors.size(); i++) {
        min = fmin(min, facecolors[i]);
        max = fmax(max, facecolors[i]);
    }

    for (unsigned int i = 0; i < facecolors.size(); i++) {
        facecolors[i] = (facecolors[i] - min) / (max - min);
    }

    jsoncons::json output;
    output["colors"] = facecolors;

    std::ofstream file(filename);
    file << output;
    file.close();
}

void gen::MapGenerator::outputDrawData(std::string filename) {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }

    std::vector<std::vector<double> > contourData;
    _getContourDrawData(contourData);

    std::vector<std::vector<double> > riverData;
    _getRiverDrawData(riverData);

    std::vector<double> slopeData;
    _getSlopeDrawData(slopeData);

    double width = _extents.maxx - _extents.minx;
    double height = _extents.maxy - _extents.miny;
    double aspectratio = width / height;

    jsoncons::json output;
    output["aspect_ratio"] = aspectratio;
    output["contour"] = contourData;
    output["river"] = riverData;
    output["slope"] = slopeData;

    std::ofstream file(filename);
    file << output;
    file.close();
}

jsoncons::json gen::MapGenerator::_getExtentsJSON() {
    jsoncons::json extents;
    extents["minx"] = _extents.minx;
    extents["miny"] = _extents.miny;
    extents["maxx"] = _extents.maxx;
    extents["maxy"] = _extents.maxy;

    return extents;
}

void gen::MapGenerator::_outputVertices(VertexList &verts, 
                                        std::string filename) {
    std::vector<double> coordinates;
    coordinates.reserve(2*_vertexMap.edge.size());
    for (unsigned int i = 0; i < verts.size(); i++) {
        coordinates.push_back(verts[i].position.x);
        coordinates.push_back(verts[i].position.y);
    }

    jsoncons::json output;
    output["vertices"] = coordinates;
    output["extents"] = _getExtentsJSON();

    std::ofstream file(filename);
    file << output;
    file.close();
}

std::vector<double> gen::MapGenerator::_computeFaceHeights(NodeMap<double> &heightMap) {
    std::vector<double> faceheights;
    faceheights.reserve(_voronoi.faces.size());

    dcel::Face f;
    dcel::Vertex v;
    std::vector<dcel::HalfEdge> edges;
    for (unsigned int i = 0; i < _voronoi.faces.size(); i++) {
        f = _voronoi.faces[i];

        edges.clear();
        _voronoi.getOuterComponents(f, edges);

        double sum = 0.0;
        for (unsigned int eidx = 0; eidx < edges.size(); eidx++) {
            v = _voronoi.origin(edges[eidx]);
            if (heightMap.isNode(v)) {
                sum += heightMap(v);
            }
        }
        double avg = sum / edges.size();
        faceheights.push_back(avg);
    }

    return faceheights;
}

std::vector<dcel::Point> gen::MapGenerator::_computeFacePositions() {
    std::vector<dcel::Point> positions;
    positions.reserve(_voronoi.faces.size());

    dcel::Face f;
    dcel::Point p;
    std::vector<dcel::HalfEdge> edges;
    for (unsigned int i = 0; i < _voronoi.faces.size(); i++) {
        f = _voronoi.faces[i];

        edges.clear();
        _voronoi.getOuterComponents(f, edges);

        double sumx = 0.0;
        double sumy = 0.0;
        for (unsigned int eidx = 0; eidx < edges.size(); eidx++) {
            p = _voronoi.origin(edges[eidx]).position;
            sumx += p.x;
            sumy += p.y;
        }
        positions.push_back(dcel::Point(sumx / edges.size(), 
                                        sumy / edges.size()));
    }

    return positions;
}

bool gen::MapGenerator::_isEdgeInMap(dcel::HalfEdge &h) {
    dcel::Vertex v1 = _voronoi.origin(h);
    dcel::Vertex v2 = _voronoi.origin(_voronoi.twin(h));

    return _vertexMap.isVertex(v1) && _vertexMap.isVertex(v2);
}

bool gen::MapGenerator::_isContourEdge(dcel::HalfEdge &h, 
                                       std::vector<double> &faceheights, 
                                       double isolevel) {
    dcel::Face f1 = _voronoi.incidentFace(h);
    dcel::Face f2 = _voronoi.incidentFace(_voronoi.twin(h));
    double iso1 = faceheights[f1.id.ref];
    double iso2 = faceheights[f2.id.ref];
    bool hasInside = iso1 < isolevel || iso2 < isolevel;
    bool hasOutside = iso1 >= isolevel || iso2 >= isolevel;

    return hasInside && hasOutside;
}

void gen::MapGenerator::_calculateErosionMap(NodeMap<double> &erosionMap) {
    _fillDepressions();

    NodeMap<double> fluxMap(&_vertexMap, 0.0);
    _calculateFluxMap(fluxMap);

    NodeMap<double> slopeMap(&_vertexMap, 0.0);
    _calculateSlopeMap(slopeMap);

    for (unsigned int i = 0; i < erosionMap.size(); i++) {
        double flux = fluxMap(i);
        double slope = slopeMap(i);
        double river = _erosionRiverFactor * sqrt(flux) * slope;
        double creep = _erosionCreepFactor * slope * slope;
        double erosion = fmin(river + creep, _maxErosionRate);
        erosionMap.set(i, erosion);
    }

    erosionMap.normalize();
}

void gen::MapGenerator::_fillDepressions() {
    NodeMap<double> finalHeightMap(&_vertexMap, _heightMap.max());
    dcel::Vertex v;
    for (unsigned int i = 0; i < _vertexMap.edge.size(); i++) {
        v = _vertexMap.edge[i];
        finalHeightMap.set(v, _heightMap(v));
    }

    double eps = 1e-5;
    std::vector<double> nbs;
    for (;;) {
        bool heightUpdated = false;
        for (unsigned int i = 0; i < _heightMap.size(); i++) {
            if (_heightMap(i) == finalHeightMap(i)) {
                continue;
            }

            nbs.clear();
            finalHeightMap.getNeighbours(i, nbs);
            for (unsigned int nidx = 0; nidx < nbs.size(); nidx++) {
                if (_heightMap(i) >= nbs[nidx] + eps) {
                    finalHeightMap.set(i, _heightMap(i));
                    heightUpdated = true;
                    break;
                }

                double hval = nbs[nidx] + eps;
                if ((finalHeightMap(i) > hval) && (hval > _heightMap(i))) {
                    finalHeightMap.set(i, hval);
                    heightUpdated = true;
                }
            }
        }

        if (!heightUpdated) {
            break;
        }
    }

    _heightMap = finalHeightMap;
}

void gen::MapGenerator::_calculateFlowMap(NodeMap<int> &flowMap) {
    dcel::Vertex v, n;
    VertexList nbs;
    for (unsigned int i = 0; i < _vertexMap.interior.size(); i++) {
        v = _vertexMap.interior[i];

        nbs.clear();
        _vertexMap.getNeighbours(v, nbs);
        dcel::Vertex minVertex;
        double minHeight = _heightMap(v);
        for (unsigned int nidx = 0; nidx < nbs.size(); nidx++) {
            n = nbs[nidx];
            if (!_vertexMap.isVertex(n)) {
                continue;
            }

            if (_heightMap(n) < minHeight) {
                minHeight = _heightMap(n);
                minVertex = n;
            }
        }

        flowMap.set(v, flowMap.getNodeIndex(minVertex));
    }
}

void gen::MapGenerator::_calculateFluxMap(NodeMap<double> &fluxMap) {
    NodeMap<int> flowMap(&_vertexMap, -1);
    _calculateFlowMap(flowMap);

    for (unsigned int i = 0; i < flowMap.size(); i++) {
        int next = i;
        while (next != -1) {
            fluxMap.set(next, fluxMap(next) + 1.0);
            next = flowMap(next);
        }
    }

    double maxFlux = _calculateFluxCap(fluxMap);
    for (unsigned int i = 0; i < fluxMap.size(); i++) {
        double f = fluxMap(i);
        f = fmin(maxFlux, f);
        f /= maxFlux;
        fluxMap.set(i, f);
    }

    _fluxMap = fluxMap;
    _flowMap = flowMap;
}

double gen::MapGenerator::_calculateFluxCap(NodeMap<double> &fluxMap) {
    double max = fluxMap.max();

    int nbins = 1000;
    int bins[1000];
    for (int i = 0; i < nbins; i++) {
        bins[i] = 0;
    }

    double step = (double)max / (double)nbins;
    double invstep = 1.0 / step;
    for (unsigned int i = 0; i < fluxMap.size(); i++) {
        double f = fluxMap(i);
        int binidx = floor(f * invstep);
        bins[binidx]++;
    }

    double acc = 0.0;
    double maxflux = 0.0;
    for (int i = 0; i < nbins; i++) {
        double pct = (double)bins[i] / (double)fluxMap.size();
        acc += pct;
        if (acc > _fluxCapPercentile) {
            maxflux = (i + 1)*step;
            break;
        }
    }

    return maxflux;
}

void gen::MapGenerator::_calculateSlopeMap(NodeMap<double> &slopeMap) {
    for (unsigned int i = 0; i < slopeMap.size(); i++) {
        slopeMap.set(i, _calculateSlope(i));
    }
}

double gen::MapGenerator::_calculateSlope(int i) {
    dcel::Vertex v = _vertexMap.vertices[i];
    if (!_vertexMap.isInterior(v)) {
        return 0.0;
    }

    double nx, ny, nz;
    _calculateVertexNormal(i, &nx, &ny, &nz);

    double slope = sqrt(nx*nx + ny*ny);

    return slope;
}

void gen::MapGenerator::_getContourDrawData(std::vector<std::vector<double> > &data) {
    std::vector<VertexList> paths;
    _getContourPaths(paths);

    double invwidth = 1.0 / (_extents.maxx - _extents.minx);
    double invheight = 1.0 / (_extents.maxy - _extents.miny);
    for (unsigned int j = 0; j < paths.size(); j++) {
        std::vector<double> drawPath;
        drawPath.reserve(2*paths[j].size());
        for (unsigned int i = 0; i < paths[j].size(); i++) {
            dcel::Vertex v = paths[j][i];
            double nx = (v.position.x - _extents.minx) * invwidth;
            double ny = (v.position.y - _extents.miny) * invheight;
            drawPath.push_back(nx);
            drawPath.push_back(ny);
        }
        data.push_back(drawPath);
    }
}

void gen::MapGenerator::_getContourPaths(std::vector<VertexList> &paths) {
    std::vector<bool> isLandFace;
    _getLandFaces(isLandFace);

    std::vector<int> adjacentEdgeCounts(_vertexMap.vertices.size(), 0);
    std::vector<bool> isEdgeVisited(_voronoi.edges.size(), false);
    dcel::HalfEdge h;
    dcel::Vertex v1, v2;
    for (unsigned int i = 0; i < _voronoi.edges.size(); i++) {
        h = _voronoi.edges[i];
        if (!_isEdgeInMap(h) || isEdgeVisited[h.id.ref]) { 
            continue; 
        }

        if (!_isContourEdge(h, isLandFace)) { 
            continue; 
        }

        v1 = _voronoi.origin(h);
        v2 = _voronoi.origin(_voronoi.twin(h));
        int idx1 = _vertexMap.getVertexIndex(v1);
        int idx2 = _vertexMap.getVertexIndex(v2);
        adjacentEdgeCounts[idx1]++;
        adjacentEdgeCounts[idx2]++;

        isEdgeVisited[h.id.ref] = true;
        isEdgeVisited[h.twin.ref] = true;
    }

    std::vector<bool> isEndVertex(_vertexMap.vertices.size(), false);
    std::vector<bool> isContourVertex(_vertexMap.vertices.size(), false);
    for (unsigned int i = 0; i < adjacentEdgeCounts.size(); i++) {
        if (adjacentEdgeCounts[i] == 1) {
            isEndVertex[i] = true;
            isContourVertex[i] = true;
        } else if (adjacentEdgeCounts[i] == 2) {
            isContourVertex[i] = true;
        }
    }

    std::vector<bool> isVertexInContour(_vertexMap.vertices.size(), false);
    for (unsigned int i = 0; i < isEndVertex.size(); i++) {
        if (!isEndVertex[i] || isVertexInContour[i]) {
            continue;
        }

        VertexList path;
        _getContourPath(i, isContourVertex, isEndVertex, isLandFace,
                           isVertexInContour, path);
        paths.push_back(path);
    }

    for (unsigned int i = 0; i < isContourVertex.size(); i++) {
        if (!isContourVertex[i] || isVertexInContour[i]) {
            continue;
        }

        VertexList path;
        _getContourPath(i, isContourVertex, isEndVertex, isLandFace,
                           isVertexInContour, path);
        paths.push_back(path);
    }
}

void gen::MapGenerator::_getLandFaces(std::vector<bool> &isLandFace) {
    std::vector<double> faceHeights;
    _getFaceHeights(faceHeights);

    isLandFace.clear();
    isLandFace.reserve(_voronoi.faces.size());
    for (unsigned int i = 0; i < faceHeights.size(); i++) {
        isLandFace.push_back(_isLand(faceHeights[i]));
    }

    _cleanupLandFaces(isLandFace);    
}

void gen::MapGenerator::_getFaceHeights(std::vector<double> &faceHeights) {
    faceHeights.clear();
    faceHeights.reserve(_voronoi.faces.size());

    dcel::Face f;
    dcel::Vertex v;
    std::vector<dcel::HalfEdge> edges;
    for (unsigned int i = 0; i < _voronoi.faces.size(); i++) {
        f = _voronoi.faces[i];

        edges.clear();
        _voronoi.getOuterComponents(f, edges);

        double sum = 0.0;
        for (unsigned int eidx = 0; eidx < edges.size(); eidx++) {
            v = _voronoi.origin(edges[eidx]);
            if (_heightMap.isNode(v)) {
                sum += _heightMap(v);
            }
        }
        double avg = sum / edges.size();
        faceHeights.push_back(avg);
    }
}

bool gen::MapGenerator::_isLand(double isolevel) {
    return isolevel >= _isolevel;
}

void gen::MapGenerator::_cleanupLandFaces(std::vector<bool> &isLandFace) {
    std::vector<std::vector<int> > islands;
    std::vector<bool> isFaceProcessed(_voronoi.faces.size(), false);
    std::vector<int> connectedFaces;
    for (unsigned int i = 0; i < isLandFace.size(); i++) {
        if (isFaceProcessed[i]) {
            continue;
        }

        connectedFaces.clear();
        _getConnectedFaces(i, isLandFace, isFaceProcessed, connectedFaces);
        islands.push_back(connectedFaces);
    }

    for (unsigned int j = 0; j < islands.size(); j++) {
        if (islands[j].size() >= _minIslandFaceThreshold) {
            continue;
        }

        for (unsigned int i = 0; i < islands[j].size(); i++) {
            int fidx = islands[j][i];
            isLandFace[fidx] = !isLandFace[fidx];
        }
    }
}

void gen::MapGenerator::_getConnectedFaces(int seed, 
                                           std::vector<bool> &isLandFace,
                                           std::vector<bool> &isFaceProcessed,
                                           std::vector<int> &faces) {
    std::vector<int> queue;
    queue.push_back(seed);
    isFaceProcessed[seed] = true;

    dcel::HalfEdge h;
    std::vector<dcel::HalfEdge> outerComponents;
    while (!queue.empty()) {
        int fidx = queue.back();
        queue.pop_back();

        bool faceType = isLandFace[fidx];
        outerComponents.clear();
        _voronoi.getOuterComponents(_voronoi.faces[fidx], outerComponents);
        for (unsigned int nidx = 0; nidx < outerComponents.size(); nidx++) {
            h = _voronoi.twin(outerComponents[nidx]);
            int nfidx = h.incidentFace.ref;
            if ((nfidx != -1) && (isLandFace[nfidx] == faceType) && !isFaceProcessed[nfidx]) {
                queue.push_back(nfidx);
                isFaceProcessed[nfidx] = true;
            }
        }

        faces.push_back(fidx);
    }
}


bool gen::MapGenerator::_isContourEdge(dcel::HalfEdge &h, 
                                       std::vector<bool> &isLandFace) {
    dcel::Face f1 = _voronoi.incidentFace(h);
    dcel::Face f2 = _voronoi.incidentFace(_voronoi.twin(h));
    return (isLandFace[f1.id.ref] && !isLandFace[f2.id.ref]) ||
           (isLandFace[f2.id.ref] && !isLandFace[f1.id.ref]);
}

bool gen::MapGenerator::_isContourEdge(dcel::Vertex &v1, dcel::Vertex &v2, 
                                       std::vector<bool> &isLandFace) {
    std::vector<dcel::HalfEdge> edges;
    edges.reserve(3);
    _voronoi.getIncidentEdges(v1, edges);

    dcel::HalfEdge h;
    dcel::Vertex v;
    for (unsigned int i = 0; i < edges.size(); i++) {
        h = edges[i];
        v = _voronoi.origin(_voronoi.twin(h));
        if (v.id.ref == v2.id.ref) {
            return _isContourEdge(h, isLandFace);
        }
    }

    return false;
}

void gen::MapGenerator::_getContourPath(int seed, std::vector<bool> &isContourVertex, 
                                                  std::vector<bool> &isEndVertex, 
                                                  std::vector<bool> &isLandFace,
                                                  std::vector<bool> &isVertexInContour, 
                                                  VertexList &path) {
    dcel::Vertex v = _vertexMap.vertices[seed];
    dcel::Vertex lastVertex = v;

    std::vector<dcel::Vertex> nbs;
    for (;;) {
        path.push_back(v);
        isVertexInContour[_vertexMap.getVertexIndex(v)] = true;
        
        nbs.clear();
        _vertexMap.getNeighbours(v, nbs);
        bool isFound = false;
        for (unsigned int i = 0; i < nbs.size(); i++) {
            dcel::Vertex n = nbs[i];
            int nidx = _vertexMap.getVertexIndex(n);
            if (n.id.ref != lastVertex.id.ref && isContourVertex[nidx] && 
                    _isContourEdge(v, n, isLandFace)) {
                lastVertex = v;
                v = n;
                isFound = true;
                break;
            }
        }

        if (!isFound) {
            break;
        }

        int vidx = _vertexMap.getVertexIndex(v);
        if (isEndVertex[vidx] || isVertexInContour[vidx]) {
            path.push_back(v);
            isVertexInContour[vidx] = true;
            break;
        }
    }
}

void gen::MapGenerator::_getRiverDrawData(std::vector<std::vector<double> > &data) {
    std::vector<VertexList> riverPaths;
    _getRiverPaths(riverPaths);

    for (unsigned int i = 0; i < riverPaths.size(); i++) {
        riverPaths[i] = _smoothPath(riverPaths[i], _riverSmoothingFactor);
    }

    double invwidth = 1.0 / (_extents.maxx - _extents.minx);
    double invheight = 1.0 / (_extents.maxy - _extents.miny);
    for (unsigned int j = 0; j < riverPaths.size(); j++) {
        std::vector<double> drawPath;
        drawPath.reserve(2*riverPaths[j].size());
        for (unsigned int i = 0; i < riverPaths[j].size(); i++) {
            dcel::Vertex v = riverPaths[j][i];
            double nx = (v.position.x - _extents.minx) * invwidth;
            double ny = (v.position.y - _extents.miny) * invheight;
            drawPath.push_back(nx);
            drawPath.push_back(ny);
        }
        data.push_back(drawPath);
    }
}

void gen::MapGenerator::_getRiverPaths(std::vector<VertexList> &riverPaths) {
    VertexList riverVertices;
    _getRiverVertices(riverVertices);

    std::vector<bool> isVertexInRiver(_vertexMap.vertices.size(), false);
    for (unsigned int i = 0; i < riverVertices.size(); i++) {
        int idx = _vertexMap.getVertexIndex(riverVertices[i]);
        isVertexInRiver[idx] = true;
    }

    VertexList fixedVertices;
    _getFixedRiverVertices(riverVertices, fixedVertices);
    std::vector<bool> isFixedVertex(_vertexMap.vertices.size(), false);
    for (unsigned int i = 0; i < fixedVertices.size(); i++) {
        int vidx = _vertexMap.getVertexIndex(fixedVertices[i]);
        isFixedVertex[vidx] = true;
    }

    VertexList path;
    for (unsigned int i = 0; i < isFixedVertex.size(); i++) {
        if (!isFixedVertex[i]) {
            continue;
        }

        path.clear();
        int next = i;
        while (isVertexInRiver[next]) {
            path.push_back(_vertexMap.vertices[next]);
            next = _flowMap(next);

            if (next == -1) { break; }
            if (isFixedVertex[next]) {
                path.push_back(_vertexMap.vertices[next]);
                break;
            }
        }

        if (path.size() >= 2) {
            riverPaths.push_back(path);
        }
    }
}

void gen::MapGenerator::_getRiverVertices(VertexList &vertices) {
    std::vector<bool> isVertexAdded(_vertexMap.vertices.size(), false);

    std::vector<bool> isLandFace;
    _getLandFaces(isLandFace);

    dcel::Vertex v;
    VertexList pathVertices;
    for (unsigned int i = 0; i < _vertexMap.vertices.size(); i++) {
        v = _vertexMap.vertices[i];
        if (_fluxMap(v) < _riverFluxThreshold || _isCoastVertex(i, isLandFace)) {
            continue;
        }

        int next = _flowMap.getNodeIndex(v);
        pathVertices.clear();
        while (next != -1) {
            pathVertices.push_back(_vertexMap.vertices[next]);
            if (_isCoastVertex(next, isLandFace)) { 
                break; 
            }
            if (!_isLandVertex(next, isLandFace)) { 
                pathVertices.clear();
                break; 
            }
            next = _flowMap(next);
        }

        if (pathVertices.empty()) { continue; }

        for (unsigned int i = 0; i < pathVertices.size(); i++) {
            int idx = _vertexMap.getVertexIndex(pathVertices[i]);
            if (!isVertexAdded[idx]) {
                vertices.push_back(pathVertices[i]);
                isVertexAdded[idx] = true;
            }
        }
    }
}

bool gen::MapGenerator::_isLandVertex(int vidx, std::vector<bool> &isLandFace) {
    std::vector<dcel::Face> faces;
    faces.reserve(6);
    _voronoi.getIncidentFaces(_vertexMap.vertices[vidx], faces);

    for (unsigned int i = 0; i < faces.size(); i++) {
        if (isLandFace[faces[i].id.ref]) {
            return true;
        }
    }

    return false;
}

bool gen::MapGenerator::_isCoastVertex(int vidx, std::vector<bool> &isLandFace) {
    std::vector<dcel::Face> faces;
    faces.reserve(6);
    _voronoi.getIncidentFaces(_vertexMap.vertices[vidx], faces);

    bool hasLand = false;
    bool hasSea = false;
    for (unsigned int i = 0; i < faces.size(); i++) {
        if (isLandFace[faces[i].id.ref]) {
            hasLand = true;
        } else {
            hasSea = true;
        }
    }

    return hasLand && hasSea;
}

void gen::MapGenerator::_getFixedRiverVertices(VertexList &riverVertices, 
                                               VertexList &fixedVertices) {
    std::vector<bool> isVertexInRiver(_vertexMap.vertices.size(), false);
    for (unsigned int i = 0; i < riverVertices.size(); i++) {
        int idx = _vertexMap.getVertexIndex(riverVertices[i]);
        isVertexInRiver[idx] = true;
    }

    std::vector<bool> isEdgeProcessed(_voronoi.edges.size(), false);
    std::vector<int> adjacentEdgeCounts(_vertexMap.vertices.size(), 0);

    dcel::HalfEdge h;
    dcel::Vertex p1, p2;
    for (unsigned int i = 0; i < _voronoi.edges.size(); i++) {
        h = _voronoi.edges[i];
        if (!_isEdgeInMap(h) || isEdgeProcessed[h.id.ref]) { 
            continue; 
        }

        p1 = _voronoi.origin(h);
        p2 = _voronoi.origin(_voronoi.twin(h));
        int idx1 = _vertexMap.getVertexIndex(p1);
        int idx2 = _vertexMap.getVertexIndex(p2);

        if (!isVertexInRiver[idx1] || !isVertexInRiver[idx2]) {
            continue;
        }

        adjacentEdgeCounts[idx1]++;
        adjacentEdgeCounts[idx2]++;

        isEdgeProcessed[h.id.ref] = true;
        isEdgeProcessed[h.twin.ref] = true;
    }

    for (unsigned int i = 0; i < adjacentEdgeCounts.size(); i++) {
        if (adjacentEdgeCounts[i] == 1 || adjacentEdgeCounts[i] == 3) {
            fixedVertices.push_back(_vertexMap.vertices[i]);
        }
    }
}

gen::MapGenerator::VertexList gen::MapGenerator::_smoothPath(VertexList &path,
                                                             double factor) {
    if (path.size() < 2) {
        return path;
    }

    VertexList smoothedPath = path;
    for (unsigned int i = 1; i < path.size() - 1; i++) {
        dcel::Point v0 = path[i-1].position;
        dcel::Point v1 = path[i].position;
        dcel::Point v2 = path[i+1].position;

        v1.x = (1-factor)*v1.x + factor*0.5*(v0.x + v2.x);
        v1.y = (1-factor)*v1.y + factor*0.5*(v0.y + v2.y);

        smoothedPath[i].position = v1;
    }

    return smoothedPath;
}

void gen::MapGenerator::_getSlopeDrawData(std::vector<double> &data) {
    std::vector<Segment> slopeSegments;
    _getSlopeSegments(slopeSegments);

    double invwidth = 1.0 / (_extents.maxx - _extents.minx);
    double invheight = 1.0 / (_extents.maxy - _extents.miny);
    for (unsigned int i = 0; i < slopeSegments.size(); i++) {
        Segment s = slopeSegments[i];

        double nx1 = (s.p1.x - _extents.minx) * invwidth;
        double ny1 = (s.p1.y - _extents.miny) * invheight;
        double nx2 = (s.p2.x - _extents.minx) * invwidth;
        double ny2 = (s.p2.y - _extents.miny) * invheight;

        data.push_back(nx1);
        data.push_back(ny1);
        data.push_back(nx2);
        data.push_back(ny2);
    }
}

void gen::MapGenerator::_getSlopeSegments(std::vector<Segment> &segments) {
    NodeMap<double> slopeMap(&_vertexMap, 0.0);
    _calculateHorizontalSlopeMap(slopeMap);

    NodeMap<double> nearSlopeMap(&_vertexMap, 0.0);
    _calculateVerticalSlopeMap(nearSlopeMap);

    std::vector<double> faceSlopes = _computeFaceHeights(slopeMap);
    std::vector<double> nearSlopes = _computeFaceHeights(nearSlopeMap);
    std::vector<dcel::Point> facePositions = _computeFacePositions();

    std::vector<bool> isLandFace;
    _getLandFaces(isLandFace);

    for (unsigned int i = 0; i < faceSlopes.size(); i++) {
        double slope = faceSlopes[i];
        if (!isLandFace[i] || fabs(slope) < _minSlopeThreshold) {
            continue;
        }

        double factor = (fabs(slope) - _minSlope) / (_maxSlope - _minSlope);
        factor = fmin(1.0, factor);
        factor = fmax(0.0, factor);

        double angle = _minSlopeAngle + factor * (_maxSlopeAngle - _minSlopeAngle);
        angle = slope < 0 ? angle : -angle;

        double dirx =  cos(angle);
        double diry =  sin(angle);
        
        double minlength = _minSlopeLength*_resolution;
        double maxlength = _maxSlopeLength*_resolution;
        double minv = _minVerticalSlope;
        double maxv = _maxVerticalSlope;
        double vslope = nearSlopes[i];
        double nf = (vslope - minv) / (maxv - minv);
        nf = fmin(1.0, nf);
        nf = fmax(0.0, nf);
        double length = minlength + nf * (maxlength - minlength);

        Segment s;
        s.p1 = facePositions[i];
        s.p2 = dcel::Point(s.p1.x + dirx*length, s.p1.y + diry*length);
        segments.push_back(s);
    }
}

void gen::MapGenerator::_calculateHorizontalSlopeMap(NodeMap<double> &slopeMap) {
    for (unsigned int i = 0; i < slopeMap.size(); i++) {
        slopeMap.set(i, _calculateHorizontalSlope(i));
    }
}

double gen::MapGenerator::_calculateHorizontalSlope(int i) {
    dcel::Vertex v = _vertexMap.vertices[i];
    if (!_vertexMap.isInterior(v)) {
        return 0.0;
    }

    double nx, ny, nz;
    _calculateVertexNormal(i, &nx, &ny, &nz);

    return nx;
}

void gen::MapGenerator::_calculateVerticalSlopeMap(NodeMap<double> &slopeMap) {
    for (unsigned int i = 0; i < slopeMap.size(); i++) {
        slopeMap.set(i, _calculateVerticalSlope(i));
    }
}

double gen::MapGenerator::_calculateVerticalSlope(int i) {
    dcel::Vertex v = _vertexMap.vertices[i];
    if (!_vertexMap.isInterior(v)) {
        return 0.0;
    }

    double nx, ny, nz;
    _calculateVertexNormal(i, &nx, &ny, &nz);

    return ny;
}

void gen::MapGenerator::_calculateVertexNormal(int vidx, 
                                                 double *nx, double *ny, double *nz) {
    dcel::Vertex v = _vertexMap.vertices[vidx];
    VertexList nbs;
    nbs.reserve(3);
    _vertexMap.getNeighbours(v, nbs);
    if (nbs.size() != 3) {
        return;
    }

    dcel::Point p0 = nbs[0].position;
    dcel::Point p1 = nbs[1].position;
    dcel::Point p2 = nbs[2].position;

    double v0x = p1.x - p0.x;
    double v0y = p1.y - p0.y;
    double v0z = _heightMap(nbs[1]) - _heightMap(nbs[0]);
    double v1x = p2.x - p0.x;
    double v1y = p2.y - p0.y;
    double v1z = _heightMap(nbs[2]) - _heightMap(nbs[0]);

    double vnx = v0y*v1z - v0z*v1y;
    double vny = v0z*v1x - v0x*v1z;
    double vnz = v0x*v1y - v0y*v1x;
    double invlen = 1.0 / sqrt(vnx*vnx + vny*vny + vnz*vnz);
    *nx = vnx * invlen;
    *ny = vny * invlen;
    *nz = vnz * invlen;
}