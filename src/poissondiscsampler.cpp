#include "poissondiscsampler.h"

/*
    Poisson disc sampling method based on "Fast Poisson Disk Sampling in 
    Arbitrary Dimensions" by Robert Bridson.

    https://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf
*/

std::vector<dcel::Point> PoissonDiscSampler::generateSamples(Extents2d bounds, 
                                                             double r, int k) {

    double dx = r / sqrt(2);
    SampleGrid grid(bounds, dx);

    Point seed = _randomPoint(bounds);
    std::vector<Point> points;
    points.push_back(seed);

    std::vector<int> activeList;
    activeList.push_back(0);

    GridIndex g = grid.getCell(seed);
    grid.setSample(g, 0);

    while (!activeList.empty()) {
        int randidx = _randomRange(0, activeList.size());
        int pidx = activeList[randidx];
        Point p = points[pidx];

        Point newPoint;
        bool isFound = _findDiscPoint(p, r, k, points, grid, &newPoint);

        if (!isFound) {
            activeList.erase(activeList.begin() + randidx);
            continue;
        }

        int newidx = points.size();
        activeList.push_back(newidx);
        points.push_back(newPoint);

        GridIndex g = grid.getCell(newPoint);
        grid.setSample(g, newidx);
    }
    
    return points;
}

double PoissonDiscSampler::_randomDouble(double min, double max) {
    return min + (double)rand() / ((double)RAND_MAX / (max - min));
}

int PoissonDiscSampler::_randomRange(int min, int max) {
    return min + (rand() % (int)(max - min));
}

dcel::Point PoissonDiscSampler::_randomPoint(Extents2d &extents) {
    double px = _randomDouble(extents.minx, extents.maxx);
    px = _randomDouble(extents.minx, extents.maxx);
    double py = _randomDouble(extents.miny, extents.maxy);

    return Point(px, py);
}

dcel::Point PoissonDiscSampler::_randomDiscPoint(dcel::Point &center, double r) {
    double angle = _randomDouble(0, 2*3.141592653);
    double nx = sin(angle);
    double ny = cos(angle);
    double rl = _randomDouble(r, 2*r);

    return Point(center.x + nx*rl, center.y + ny*rl);
}

bool PoissonDiscSampler::_findDiscPoint(dcel::Point &center, double r, int k, 
                                        std::vector<dcel::Point> &points, 
                                        SampleGrid &grid, dcel::Point *p) {
    for (int i = 0; i < k; i++) {
        Point sample = _randomDiscPoint(center, r);
        if (!grid.bounds.containsPoint(sample.x, sample.y)) {
            continue;
        }

        if (_isSampleValid(sample, r, points, grid)) {
            *p = sample;
            return true;
        }
    }

    return false;
}

bool PoissonDiscSampler::_isSampleValid(dcel::Point &p, double r, 
                                        std::vector<dcel::Point> &points, 
                                        SampleGrid &grid) {
    GridIndex g = grid.getCell(p);
    int sampleid = grid.getSample(g);
    if (sampleid != -1) {
        return false;
    }

    int mini = fmax(g.i - 2, 0);
    int minj = fmax(g.j - 2, 0);
    int maxi = fmin(g.i + 2, grid.width - 1);
    int maxj = fmin(g.j + 2, grid.height - 1);

    double rsq = r*r;
    for (int j = minj; j <= maxj; j++) {
        for (int i = mini; i <= maxi; i++) {
            int sampleid = grid.getSample(i, j);
            if (sampleid == -1) {
                continue;
            }

            Point o = points[sampleid];
            double dx = p.x - o.x;
            double dy = p.y - o.y;
            double distsq = dx*dx + dy*dy;
            if (distsq < rsq) {
                return false;
            }
        }
    }

    return true;
}
