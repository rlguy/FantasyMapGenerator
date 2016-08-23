#include <stdio.h>
#include <iostream>
#include <vector>

#include <stdlib.h>
#include <time.h>
#include <fstream>

#include "dcel.h"
#include "voronoi.h"

#include "jsoncons/json.hpp"

double randomDouble(double min, double max) {
    return min + (double)rand() / ((double)RAND_MAX / (max - min));
}

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
    srand(time(NULL));

    // Generate an n point Voronoi diagram and write the face vertices
    // to a JSON file.
    int n = 10000;
    double minx = -20.0;
    double maxx = 20.0;
    double miny = -20.0;
    double maxy = 20.0;
    std::vector<dcel::Point> points;
    points.reserve(n);
    for (int i = 0; i < n; i++) {
        dcel::Point p(randomDouble(minx, maxx), randomDouble(miny, maxy));
        points.push_back(p);
    }

    std::cout << "Generating Voronoi diagram with " << n << " points." << std::endl;

    dcel::DCEL T = Voronoi::voronoi(points);

    std::cout << "Finished generating Voronoi diagram." << std::endl;
    std::cout << "# Vertices:   " << T.vertices.size() << std::endl;
    std::cout << "# Half Edges: " << T.edges.size() << std::endl;
    std::cout << "# Faces:      " << T.faces.size() << std::endl;

    outputDCEL(T, "output.json");

    return 0;
}
