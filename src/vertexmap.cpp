#include "vertexmap.h"

gen::VertexMap::VertexMap() {}

gen::VertexMap::VertexMap(dcel::DCEL *V, Extents2d extents) :
                        _dcel(V), _extents(extents) {
    vertices.reserve(_dcel->vertices.size());
    interior.reserve(_dcel->vertices.size());
    _vertexTypes.reserve(_dcel->vertices.size());
    _vertexIdToMapIndex = std::vector<int>(_dcel->vertices.size(), -1);

    dcel::Vertex v;
    std::vector<dcel::Vertex> neighbours;
    for (unsigned int i = 0; i < _dcel->vertices.size(); i++) {
        v = _dcel->vertices[i];
        if (!extents.containsPoint(v.position) || _isBoundaryVertex(v)) {
            continue;
        }

        vertices.push_back(v);
        _vertexIdToMapIndex[v.id.ref] = vertices.size() - 1;

        neighbours.clear();
        getNeighbours(v, neighbours);
        if (neighbours.size() < 3) {
            edge.push_back(v);
            _vertexTypes.push_back(VertexType::edge);
        } else {
            interior.push_back(v);
            _vertexTypes.push_back(VertexType::interior);
        }
    }
}

unsigned int gen::VertexMap::size() {
    return vertices.size();
}

void gen::VertexMap::getNeighbours(dcel::Vertex v, 
                                   std::vector<dcel::Vertex> &nbs) {
    dcel::HalfEdge h = _dcel->incidentEdge(v);
    dcel::Ref startRef = h.id;

    dcel::HalfEdge twin;
    dcel::Vertex n;
    do {
        twin = _dcel->twin(h);
        n = _dcel->origin(twin);
        if (_extents.containsPoint(n.position)) {
            nbs.push_back(n);
        }
        h = _dcel->next(twin);
    } while (h.id != startRef);
}

void gen::VertexMap::getNeighbourIndices(dcel::Vertex v, std::vector<int> &nbs) {
    dcel::HalfEdge h = _dcel->incidentEdge(v);
    dcel::Ref startRef = h.id;

    dcel::HalfEdge twin;
    dcel::Vertex n;
    do {
        twin = _dcel->twin(h);
        n = _dcel->origin(twin);
        if (_extents.containsPoint(n.position)) {
            nbs.push_back(getVertexIndex(n));
        }
        h = _dcel->next(twin);
    } while (h.id != startRef);
}

int gen::VertexMap::getVertexIndex(dcel::Vertex &v) {
    if (!_isInRange(v.id.ref)) {
        throw std::range_error("Index out of range:" + v.id.ref);
    }
    return _vertexIdToMapIndex[v.id.ref];
}

bool gen::VertexMap::isVertex(dcel::Vertex &v) {
    if (!_isInRange(v.id.ref)) {
        return false;
    }
    return getVertexIndex(v) != -1;
}

bool gen::VertexMap::isEdge(dcel::Vertex &v) {
    if (!_isInRange(v.id.ref)) {
        return false;
    }
    return _vertexTypes[getVertexIndex(v)] == VertexType::edge;
}

bool gen::VertexMap::isInterior(dcel::Vertex &v) {
    if (!_isInRange(v.id.ref)) {
        return false;
    }
    return _vertexTypes[getVertexIndex(v)] == VertexType::interior;
}

bool gen::VertexMap::_isBoundaryVertex(dcel::Vertex &v) {
    if (v.incidentEdge.ref == -1) {
        return true;
    }

    dcel::HalfEdge h = _dcel->incidentEdge(v);
    dcel::Ref startid = h.id;

    do {
        if (h.incidentFace.ref == -1) {
            return true;
        }

        h = _dcel->twin(h);
        if (h.next.ref == -1) {
            return true;
        }

        h = _dcel->next(h);
    } while (h.id != startid);

    return false;
}

bool gen::VertexMap::_isInRange(int id) {
    return id >= 0 || id < (int)_vertexIdToMapIndex.size();
}
