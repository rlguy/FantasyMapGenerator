#include <stdio.h>
#include <iostream>

#include "mapgenerator.h"

double randomDouble(double min, double max) {
    return min + (double)rand() / ((double)RAND_MAX / (max - min));
}

dcel::Point randomPoint(Extents2d &extents) {
    double px = randomDouble(extents.minx, extents.maxx);
    px = randomDouble(extents.minx, extents.maxx);
    double py = randomDouble(extents.miny, extents.maxy);

    return dcel::Point(px, py);
}

dcel::Point randomDirection() {
    double angle = randomDouble(0, 2*3.14159);
    return dcel::Point(sin(angle), cos(angle));
}

int main() {
    auto seed = time(NULL);
    std::cout << "SEED " << (unsigned int)seed << std::endl;
    srand(seed);

    Extents2d extents(0, 0, 1.7777*20.0, 20.0);
    gen::MapGenerator map(extents, 0.09);
    map.initialize();
    
    int n = randomDouble(60, 120);
    double minr = 1.0;
    double maxr = 6.0;
    for (int i = 0; i < n; i++) {
        dcel::Point p = randomPoint(extents);
        double r = randomDouble(minr, maxr);
        double strength = randomDouble(0.5, 1.5);
        if (randomDouble(0, 1) > 0.5) {
            map.addHill(p.x, p.y, r, strength);
        } else {
            map.addCone(p.x, p.y, r, strength);
        }
    }

    if (randomDouble(0, 1) > 0.5) {
        dcel::Point p = randomPoint(extents);
        double r = randomDouble(6.0, 12.0);
        double strength = randomDouble(1.0, 3.0);
        map.addCone(p.x, p.y, r, strength);
    }

    if (randomDouble(0, 1) > 0.5) {
        dcel::Point dir = randomDirection();
        dcel::Point lp = randomPoint(extents);
        double slopewidth = randomDouble(0.5, 5.0);
        double strength = randomDouble(2.0, 3.0);
        map.addSlope(lp.x, lp.y, dir.x, dir.y, slopewidth, strength);
    }
    
    if (randomDouble(0, 1) > 0.5) {
        map.normalize();
    } else {
        map.round();
    }

    if (randomDouble(0, 1) > 0.5) {
        map.relax();
    }

    map.setSeaLevelToMedian();
    for (int i = 0; i < 3; i++) {
        double e = randomDouble(0.05, 0.12);
        map.erode(e);
    }
    map.setSeaLevelToMedian();

    for (int i = 0; i < 5; i++) {
        map.addCity();
    }

    map.outputDrawData("output.json");

    return 0;
}
