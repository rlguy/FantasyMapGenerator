#ifndef EXTENTS2D_H
#define EXTENTS2D_H

#include "dcel.h"

struct Extents2d {
    double minx = 0.0;
    double miny = 0.0;
    double maxx = 0.0;
    double maxy = 0.0;

    Extents2d() {}
    Extents2d(double mina, double minb, double maxa, double maxb) : 
                    minx(mina), miny(minb), maxx(maxa), maxy(maxb) {}

    bool containsPoint(dcel::Point p) {
        return containsPoint(p.x, p.y);
    }
    
    bool containsPoint(double x, double y) {
        return x >= minx && x < maxx && y >= miny && y < maxy;
    }
};

#endif