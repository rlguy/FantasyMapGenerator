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

    int getNodeIndex(dcel::Vertex &v) {
        return _vertexMap->getVertexIndex(v);
    }

    void fill(T fillval) {
        std::fill(_nodes.begin(), _nodes.end(), fillval);
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

    void set(int idx, T val) {
        if (!_isInRange(idx)) {
            throw std::range_error("Index out of range:" + idx);
        }
        _nodes[idx] = val;
    }

    void set(dcel::Vertex &v, T val) {
        int idx = getNodeIndex(v);
        if (idx == -1) {
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

private:
    void _initializeNodes() {
        _nodes = std::vector<T>(size());
    }

    bool _isInRange(int idx) {
        return idx >= 0 || idx < (int)_nodes.size();
    }

    VertexMap *_vertexMap;
    std::vector<T> _nodes;
};

}

#endif