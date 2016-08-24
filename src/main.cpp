#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>

#include "dcel.h"
#include "voronoi.h"
#include "poissondiscsampler.h"
#include "extents2d.h"

#include "jsoncons/json.hpp"

void outputDCEL(dcel::DCEL &T, std::string filename) {
    using jsoncons::json;

    std::vector<std::vector<double> > faceVertices;
    faceVertices.reserve(T.faces.size());

    for (unsigned int i = 0; i < T.faces.size(); i++) {
        dcel::Face f = T.faces[i];

        if (f.outerComponent.ref == -1) {
            continue;
        }

        dcel::HalfEdge h = T.outerComponent(f);
        dcel::Ref startRef = h.id;


        std::vector<double> vertices;
        dcel::Vertex v;
        do {
            v = T.origin(h);
            vertices.push_back(v.position.x);
            vertices.push_back(v.position.y);
            h = T.next(h);
        } while (h.id != startRef);
        faceVertices.push_back(vertices);
    }
    
    json output;
    output["faces"] = faceVertices;

    std::ofstream file(filename);
    file << output;
    file.close();
}

int main() {
    //srand(time(NULL));

    using namespace dcel;

    double minx = -20.0;
    double maxx = 20.0;
    double miny = -20.0;
    double maxy = 20.0;
    Extents2d bounds(minx, miny, maxx, maxy);
    double radius = 0.25;
    int k = 20;

    std::cout << "Generating samples within " << 
                 maxx - minx << " x " << maxy - miny << 
                 " area with radius " << radius << std::endl;

    std::vector<Point> samples;
    samples = PoissonDiscSampler::generateSamples(bounds, radius, k);

    std::cout << "Finished Generating " << samples.size() << 
                 " samples." << std::endl;

    std::cout << "Generating Voronoi diagram with " << samples.size() << 
                 " points." << std::endl;

    dcel::DCEL T = Voronoi::voronoi(samples);

    std::cout << "Finished generating Voronoi diagram." << std::endl;
    std::cout << "# Vertices:   " << T.vertices.size() << std::endl;
    std::cout << "# Half Edges: " << T.edges.size() << std::endl;
    std::cout << "# Faces:      " << T.faces.size() << std::endl;

    outputDCEL(T, "output.json");

    return 0;
}
