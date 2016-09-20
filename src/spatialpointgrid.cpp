#include "spatialpointgrid.h"

gen::SpatialPointGrid::SpatialPointGrid() {
}

gen::SpatialPointGrid::SpatialPointGrid(std::vector<dcel::Point> &points, 
                                        double dx) : _dx(dx) {
    _initializeGrid(points);
}

int gen::SpatialPointGrid::getPointCount(Extents2d extents) {
    Extents2d temp = extents;
    temp.minx -= _offset.x;
    temp.miny -= _offset.y;
    temp.maxx -= _offset.x;
    temp.maxy -= _offset.y;

    double invdx = 1.0 / _dx;
    int mini = fmax(0, floor(temp.minx * invdx));
    int minj = fmax(0, floor(temp.miny * invdx));
    int maxi = fmin(_isize - 1, floor(temp.maxx * invdx));
    int maxj = fmin(_jsize - 1, floor(temp.maxy * invdx));

    int count = 0;
    for (int j = minj; j <= maxj; j++) {
        for (int i = mini; i <= maxi; i++) {
            int grididx = _flattenIndex(i, j);
            for (unsigned int g = 0; g < _grid[grididx].size(); g++) {
                if (extents.containsPoint(_grid[grididx][g])) {
                    count++;
                }
            }
        }
    }

    return count;
}

void gen::SpatialPointGrid::_initializeGrid(std::vector<dcel::Point> &points) {
    _extents = _getPointSetExtents(points);
    _width = _extents.maxx - _extents.minx;
    _height = _extents.maxy - _extents.miny;
    _offset = dcel::Point(_extents.minx, _extents.miny);
    _isize = (int)ceil(_width / _dx);
    _jsize = (int)ceil(_height / _dx);
    _grid = std::vector<std::vector<dcel::Point> >(_isize*_jsize);

    _insertPointsIntoGrid(points);
}

Extents2d gen::SpatialPointGrid::_getPointSetExtents(
                                    std::vector<dcel::Point> &points) {

    dcel::Point p = points[0];
    Extents2d e(p.x, p.y, p.x, p.y);
    for (unsigned int i = 0; i < points.size(); i++) {
        p = points[i];
        e.minx = fmin(p.x, e.minx);
        e.miny = fmin(p.y, e.miny);
        e.maxx = fmax(p.x, e.maxx);
        e.maxy = fmax(p.y, e.maxy);
    }

    return e;
}

void gen::SpatialPointGrid::_insertPointsIntoGrid(
                                std::vector<dcel::Point> &points) {
    double invdx = 1.0 / _dx;
    for (unsigned int idx = 0; idx < points.size(); idx++) {
        int i = (int)floor((points[idx].x - _offset.x) * invdx);
        int j = (int)floor((points[idx].y - _offset.y) * invdx);
        int grididx = _flattenIndex(i, j);
        _grid[grididx].push_back(points[idx]);
    }
}

int gen::SpatialPointGrid::_flattenIndex(int i, int j) {
    return i + _isize*j;
}