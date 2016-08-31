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
    
    /* 
        Normalize height map values to range [0, 1]
    */
    double min = std::numeric_limits<double>::infinity();
    double max = -std::numeric_limits<double>::infinity();
    for (unsigned int i = 0; i < _heightMap.size(); i++) {
        min = fmin(min, _heightMap(i));
        max = fmax(max, _heightMap(i));
    }

    for (unsigned int i = 0; i < _heightMap.size(); i++) {
        double val = _heightMap(i);
        double normalized = (val - min) / (max - min);
        _heightMap.set(i, normalized);
    }
}

void gen::MapGenerator::round() {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }

    /* 
        Normalize height map and square root the values
    */
    normalize();
    for (unsigned int i = 0; i < _heightMap.size(); i++) {
        double rounded = sqrt(_heightMap(i));
        _heightMap.set(i, rounded);
    }
}

void gen::MapGenerator::relax() {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }

    /* 
        Replace height with average of its neighbours
    */
    std::vector<double> averages;
    averages.reserve(_heightMap.size());

    std::vector<double> nbs;
    for (unsigned int i = 0; i < _heightMap.size(); i++) {
        nbs.clear();
        _heightMap.getNeighbours(i, nbs);
        if (nbs.size() == 0) {
            continue;
        }

        double sum = 0.0;
        for (unsigned int nidx = 0; nidx < nbs.size(); nidx++) {
            sum += nbs[nidx];
        }
        averages.push_back(sum / nbs.size());
    }

    for (unsigned int i = 0; i < _heightMap.size(); i++) {
        _heightMap.set(i, averages[i]);
    }
}

void gen::MapGenerator::setSeaLevel(double level) {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }

    /* 
        Translate height map so that level is at 0.5
    */
    for (unsigned int i = 0; i < _heightMap.size(); i++) {
        double newval = _heightMap(i) - (level - 0.5);
        _heightMap.set(i, newval);
    }
}

void gen::MapGenerator::setSeaLevelToMedian() {
    if (!_isInitialized) {
        throw std::runtime_error("MapGenerator must be initialized.");
    }

    std::vector<double> values;
    values.reserve(_heightMap.size());
    for (unsigned int i = 0; i < _heightMap.size(); i++) {
        values.push_back(_heightMap(i));
    }

    std::sort(values.begin(), values.end());
    int mididx = values.size() / 2;
    double median;
    if (values.size() % 2 == 0) {
        median = 0.5 * (values[mididx - 1] + values[mididx]);
    } else {
        median = values[mididx];
    }

    setSeaLevel(median);
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
        _heightMap.set(i, _heightMap(i) - amount * erosionMap(i));
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

    jsoncons::json output;
    output["colors"] = facecolors;

    std::ofstream file(filename);
    file << output;
    file.close();
}

void gen::MapGenerator::outputContour(std::string filename, double isolevel) {
    jsoncons::json output;
    output["contour"] = _computeContour(isolevel);
    output["extents"] = _getExtentsJSON();

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

void gen::MapGenerator::_outputVertices(std::vector<dcel::Vertex> &verts, 
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

std::vector<double> gen::MapGenerator::_computeContour(double isolevel) {
    std::vector<double> contour;
    std::vector<double> faceheights = _computeFaceHeights(_heightMap);
    std::vector<bool> isEdgeInContour(_voronoi.edges.size(), false);

    dcel::HalfEdge h;
    dcel::Point p1, p2;
    for (unsigned int i = 0; i < _voronoi.edges.size(); i++) {
        h = _voronoi.edges[i];
        if (!_isEdgeInMap(h) || isEdgeInContour[h.id.ref]) {
            continue;
        }

        if (!_isContourEdge(h, faceheights, isolevel)) {
            continue;
        }

        p1 = _voronoi.origin(h).position;
        p2 = _voronoi.origin(_voronoi.twin(h)).position;

        contour.push_back(p1.x);
        contour.push_back(p1.y);
        contour.push_back(p2.x);
        contour.push_back(p2.y);

        isEdgeInContour[h.id.ref] = true;
        isEdgeInContour[h.twin.ref] = true;
    }

    return contour;
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

    double max = erosionMap.max();
    for (unsigned int i = 0; i < erosionMap.size(); i++) {
        erosionMap.set(i, erosionMap(i) / max);
    }
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
    std::vector<dcel::Vertex> nbs;
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

    std::vector<dcel::Vertex> nbs;
    nbs.reserve(3);
    _vertexMap.getNeighbours(v, nbs);
    if (nbs.size() != 3) {
        return 0.0;
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

    double nx = v0y*v1z - v0z*v1y;
    double ny = v0z*v1x - v0x*v1z;
    double nz = v0x*v1y - v0y*v1x;
    double invlen = 1.0 / sqrt(nx*nx + ny*ny + nz*nz);
    nx *= invlen;
    ny *= invlen;

    double slope = sqrt(nx*nx + ny*ny);

    return slope;
}