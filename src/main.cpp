#include <stdio.h>
#include <iostream>

#include "mapgenerator.h"

int main() {
    //srand(time(NULL));

    gen::MapGenerator map;
    map.initialize();
    map.outputVoronoiDiagram("voronoi.json");

    return 0;
}
