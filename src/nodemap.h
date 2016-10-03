#ifndef NODEMAP_H
#define NODEMAP_H

#include "vertexmap.h"
#include "dcel.h"

namespace gen {

template <class T>
class NodeMap {

public:
    NodeMap() {}
    NodeMap(VertexMap *vertexMap) : _vertexMap(vertexMap) {
        _initializeNodes();
    }

    NodeMap(VertexMap *vertexMap, T fillval) : _vertexMap(vertexMap) {
        _initializeNodes();
        fill(fillval);
    }

    unsigned int size() {
        return _vertexMap->size();
    }

    T min() {
        T minval = _nodes[0];
        for (unsigned int i = 0; i < _nodes.size(); i++) {
            if (_nodes[i] < minval) {
                minval = _nodes[i];
            }
        }
        return minval;
    }

    T max() {
        T maxval = _nodes[0];
        for (unsigned int i = 0; i < _nodes.size(); i++) {
            if (_nodes[i] > maxval) {
                maxval = _nodes[i];
            }
        }
        return maxval;
    }

    int getNodeIndex(dcel::Vertex &v) {
        return _vertexMap->getVertexIndex(v);
    }

    void fill(T fillval) {
        for (unsigned int i = 0; i < _nodes.size(); i++) {
            _nodes[i] = fillval;
        }
    }

    T operator()(int idx) {
        if (!_isInRange(idx)) {
            throw std::range_error("Index out of range:" + idx);
        }
        return _nodes[idx];
    }

    T operator()(dcel::Vertex v) {
        int idx = getNodeIndex(v);
        if (idx == -1) {
            throw std::range_error("Vertex is not in NodeMap.");
        }
        return _nodes[idx];
    }

    T* getPointer(int idx) {
        if (!_isInRange(idx)) {
            throw std::range_error("Index out of range:" + idx);
        }
        return &_nodes[idx];
    }

    T* getPointer(dcel::Vertex v) {
        int idx = getNodeIndex(v);
        if (idx == -1) {
            throw std::range_error("Vertex is not in NodeMap.");
        }
        return &_nodes[idx];
    }

    void set(int idx, T val) {
        if (!_isInRange(idx)) {
            throw std::range_error("Index out of range:" + idx);
        }
        _nodes[idx] = val;
    }

    void set(dcel::Vertex &v, T val) {
        int idx = getNodeIndex(v);
        if (idx == -1 || !_isInRange(idx)) {
            throw std::range_error("Vertex is not in NodeMap.");
        }
        _nodes[idx] = val;
    }

    void getNeighbours(int idx, std::vector<T> &nbs) {
        if (!_isInRange(idx)) {
            throw std::range_error("Index out of range:" + idx);
        }
        std::vector<int> indices;
        indices.reserve(3);
        _vertexMap->getNeighbourIndices(_vertexMap->vertices[idx], indices);

        for (unsigned int i = 0; i < indices.size(); i++) {
            nbs.push_back(_nodes[indices[i]]);
        }
    }

    void getNeighbours(dcel::Vertex &v, std::vector<T> &nbs) {
        int idx = getNodeIndex(v);
        if (idx == -1) {
            throw std::range_error("Vertex is not in NodeMap.");
        }
        std::vector<int> indices;
        indices.reserve(3);
        _vertexMap->getNeighbourIndices(v, indices);

        for (unsigned int i = 0; i < indices.size(); i++) {
            nbs.push_back(_nodes[indices[i]]);
        }
    }

    bool isNode(int idx) {
        return _isInRange(idx);
    }

    bool isNode(dcel::Vertex &v) {
        int idx = getNodeIndex(v);
        return idx != -1;
    }

    bool isEdge(int idx) {
        if (!_isInRange(idx)) {
            return false;
        }
        return _vertexMap->isEdge(_vertexMap->vertices[idx]);
    }

    bool isEdge(dcel::Vertex &v) {
        return _vertexMap->isEdge(v);
    }

    bool isInterior(int idx) {
        if (!_isInRange(idx)) {
            return false;
        }
        return _vertexMap->isInterior(_vertexMap->vertices[idx]);
    }

    bool isInterior(dcel::Vertex &v) {
        return _vertexMap->isInterior(v);
    }

    // Normalize height map values to range [0, 1]
    void normalize() {
        double min = std::numeric_limits<double>::infinity();
        double max = -std::numeric_limits<double>::infinity();
        for (unsigned int i = 0; i < size(); i++) {
            min = fmin(min, (*this)(i));
            max = fmax(max, (*this)(i));
        }

        for (unsigned int i = 0; i < size(); i++) {
            double val = (*this)(i);
            double normalized = (val - min) / (max - min);
            set(i, normalized);
        }
    }

    // Normalize height map and square root the values
    void round() {
        normalize();
        for (unsigned int i = 0; i < size(); i++) {
            double rounded = sqrt((*this)(i));
            set(i, rounded);
        }
    }

    //  Replace height with average of its neighbours
    void relax() {
        std::vector<double> averages;
        averages.reserve(size());

        std::vector<double> nbs;
        for (unsigned int i = 0; i < size(); i++) {
            nbs.clear();
            getNeighbours(i, nbs);
            if (nbs.size() == 0) {
                continue;
            }

            double sum = 0.0;
            for (unsigned int nidx = 0; nidx < nbs.size(); nidx++) {
                sum += nbs[nidx];
            }
            averages.push_back(sum / nbs.size());
        }

        for (unsigned int i = 0; i < size(); i++) {
            set(i, averages[i]);
        }
    }

    // Translate height map so that level is at zero
    void setLevel(double level) {
        for (unsigned int i = 0; i < size(); i++) {
            double newval = (*this)(i) - level;
            set(i, newval);
        }
    }

    void setLevelToMedian() {
        std::vector<double> values;
        values.reserve(size());
        for (unsigned int i = 0; i < size(); i++) {
            values.push_back((*this)(i));
        }

        std::sort(values.begin(), values.end());
        int mididx = values.size() / 2;
        double median;
        if (values.size() % 2 == 0) {
            median = 0.5 * (values[mididx - 1] + values[mididx]);
        } else {
            median = values[mididx];
        }

        setLevel(median);
    }
    
private:
    void _initializeNodes() {
        _nodes = std::vector<T>(size(), T());
    }

    bool _isInRange(int idx) {
        return idx >= 0 && idx < (int)_nodes.size();
    }

    VertexMap *_vertexMap;
    std::vector<T> _nodes;
};

}

#endif