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
                 " area with radius " << _resolution << std::endl;

    std::vector<dcel::Point> samples;
    samples = PoissonDiscSampler::generateSamples(_extents, _resolution, 25);

    std::cout << "Finished Generating " << samples.size() << 
                 " samples." << std::endl;

    std::cout << "Generating Voronoi diagram with " << samples.size() << 
                 " points." << std::endl;

    dcel::DCEL V = Voronoi::voronoi(samples);

    std::cout << "Finished generating Voronoi diagram." << std::endl;
    std::cout << "# Vertices:   " << V.vertices.size() << std::endl;
    std::cout << "# Half Edges: " << V.edges.size() << std::endl;
    std::cout << "# Faces:      " << V.faces.size() << std::endl;

    _voronoi = V;
    _isInitialized = true;
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

    std::ofstream file(filename);
    file << output;
    file.close();
}
