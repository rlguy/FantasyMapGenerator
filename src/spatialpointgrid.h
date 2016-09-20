#ifndef SPATIALPOINTGRID_H
#define SPATIALPOINTGRID_H

#include <stdio.h>
#include <iostream>
#include <vector>

#include "dcel.h"
#include "extents2d.h"

namespace gen {

class SpatialPointGrid {

public:
    SpatialPointGrid();
    SpatialPointGrid(std::vector<dcel::Point> &points, double dx);

    int getPointCount(Extents2d extents);

private:
    void _initializeGrid(std::vector<dcel::Point> &points);
    Extents2d _getPointSetExtents(std::vector<dcel::Point> &points);
    void _insertPointsIntoGrid(std::vector<dcel::Point> &points);
    int _flattenIndex(int i, int j);

    double _dx;
    Extents2d _extents;
    double _width;
    double _height;
    double _isize;
    double _jsize;
    dcel::Point _offset;

    std::vector<std::vector<dcel::Point> > _grid;
};

}

#endif