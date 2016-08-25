#include <stdio.h>
#include <iostream>

#include "mapgenerator.h"

int main() {
    //srand(time(NULL));

    gen::MapGenerator map(Extents2d(-20, -10, 20, 10), 0.4);
    map.initialize();
    map.outputVoronoiDiagram("voronoi.json");
    map.outputEdgeVertices("edgevertices.json");
    map.outputInteriorVertices("interiorvertices.json");

    return 0;
}
