#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <math.h>

#include "dcel.h"

namespace Geometry {

inline bool lineIntersection(dcel::Point &p, dcel::Point &r, 
                              dcel::Point &q, dcel::Point &s, 
                              dcel::Point *intersect) {
    double cross = r.x*s.y - r.y*s.x;

    double eps = 1e-9;
    if (fabs(cross) < eps) {
        // parallel case
        return false;
    }

    double vx = q.x - p.x;
    double vy = q.y - p.y;
    double t = (vx*s.y - vy*s.x) / cross;

    intersect->x = p.x + t*r.x;
    intersect->y = p.y + t*r.y;

    return true;
}

inline bool lineSegmentIntersection(dcel::Point &A, dcel::Point &B, 
                                    dcel::Point &C, dcel::Point &D) {
    bool c1 = (D.y-A.y)*(C.x-A.x) > (C.y-A.y)*(D.x-A.x);
    bool c2 = (D.y-B.y)*(C.x-B.x) > (C.y-B.y)*(D.x-B.x);
    bool c3 = (C.y-A.y)*(B.x-A.x) > (B.y-A.y)*(C.x-A.x);
    bool c4 = (D.y-A.y)*(B.x-A.x) > (B.y-A.y)*(D.x-A.x);

    return (c1 != c2) && (c3 != c4);
}

}

#endif