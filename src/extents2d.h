#ifndef EXTENTS2D_H
#define EXTENTS2D_H

struct Extents2d {
    double minx = 0.0;
    double miny = 0.0;
    double maxx = 0.0;
    double maxy = 0.0;

    Extents2d() {}
    Extents2d(double mina, double minb, double maxa, double maxb) : 
                    minx(mina), miny(minb), maxx(maxa), maxy(maxb) {}

    bool containsPoint(double x, double y) {
        return x >= minx && x < maxx && y >= miny && y < maxy;
    }
};

#endif