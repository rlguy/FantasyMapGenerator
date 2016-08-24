#ifndef POISSONDISCSAMPLER_H
#define POISSONDISCSAMPLER_H

#include <vector>
#include <cmath>
#include <stdlib.h>
#include <time.h>

#include "dcel.h"

namespace PoissonDiscSampler {
    using namespace dcel;

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

    struct GridIndex {
        int i = 0;
        int j = 0;

        GridIndex() {}
        GridIndex(int ii, int jj) : i(ii), j(jj) {}
    };

    struct SampleGrid {
        Extents2d bounds;
        int width = 0;
        int height = 0;
        double dx = 0.0;
        std::vector<int> grid;

        SampleGrid(Extents2d extents, double cellsize) {
            bounds = extents;
            dx = cellsize;
            double bw = bounds.maxx - bounds.minx;
            double bh = bounds.maxy - bounds.miny;
            width = ceil(bw / cellsize);
            height = ceil(bh / cellsize);
            grid = std::vector<int>(width*height, -1);
        }

        int getFlatIndex(int i, int j) {
            return i + j*width;
        }

        int getSample(GridIndex g) {
            return getSample(g.i, g.j);
        }

        int getSample(int i, int j) {
            if (i < 0 || i > width || j < 0 || j > height) {
                throw std::range_error("getSample");
            }
            return grid[getFlatIndex(i, j)];
        }

        void setSample(GridIndex g, int s) {
            setSample(g.i, g.j, s);
        }

        void setSample(int i, int j, int s) {
            if (i < 0 || i > width || j < 0 || j > height) {
                throw std::range_error("setSample");
            }
            grid[getFlatIndex(i, j)] = s;
        }

        GridIndex getCell(Point p) {
            return getCell(p.x, p.y);
        }

        GridIndex getCell(double x, double y) {
            x -= bounds.minx;
            y -= bounds.miny;
            return GridIndex(floor(x / dx), floor(y / dx));
        }
    };

    std::vector<Point> generateSamples(Extents2d bounds, double r, int k);

    double _randomDouble(double min, double max);
    int _randomRange(int min, int max);
    Point _randomPoint(Extents2d &extents);
    Point _randomDiscPoint(Point &center, double r);
    bool _findDiscPoint(Point &center, double r, int k, 
                        std::vector<Point> &points, SampleGrid &grid, Point *p);
    bool _isSampleValid(Point &p, double r, 
                        std::vector<Point> &points, SampleGrid &grid);

}

#endif