#ifndef DCEL_H
#define DCEL_H

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <stdexcept>
#include <exception>

namespace dcel {

class Point {
public:
	Point() {}
	Point(double px, double py) : x(px), y(py) {}

	double x = 0.0;
	double y = 0.0;
};

class Ref {
public:
	Ref() {}
	Ref(int id) : ref(id) {};

	bool operator == (const Ref &r) {
		return ref == r.ref;
	}

	bool operator != (const Ref &r) {
		return ref != r.ref;
	}

	int ref = -1;
};

class Face;
class HalfEdge {
public:
	HalfEdge() {}

	Ref origin;
	Ref twin;
	Ref incidentFace;
	Ref next;
	Ref prev;
	Ref id;
};

class Face {
public:
	Face() {}

	Ref outerComponent;
	std::vector<Ref> innerComponents;
	Ref id;
};

class Vertex {
public: 
	Vertex() {}
	Vertex(double x, double y) : position(x, y) {}

	Point position;
	Ref incidentEdge;
	Ref id;
};

class DCEL {
public:
	DCEL() {
	}

	inline Vertex createVertex(Point p) {
		return createVertex(p.x, p.y);
	}

	inline Vertex createVertex(double px, double py) {
		Vertex vert(px, py);
	    vert.id = Ref(vertices.size());
		vertices.push_back(vert);
		return vert;
	}

	inline HalfEdge createHalfEdge() {
		HalfEdge edge;
		edge.id = Ref(edges.size());
		edges.push_back(edge);
		return edge;
	}

	inline Face createFace() {
		Face face;
		face.id = Ref(faces.size());
		faces.push_back(face);
		return face;
	}

	inline Vertex getVertex(Ref &id) {
		return getVertex(id.ref);
	}

	inline Vertex getVertex(int id) {
		if (!_isVertexInRange(id)) {
			throw std::range_error("Vertex out of range: " + id);
		}
		return vertices[id];
	}

	inline HalfEdge getHalfEdge(Ref &id) {
		return getHalfEdge(id.ref);
	}

	inline HalfEdge getHalfEdge(int id) {
		if (!_isHalfEdgeInRange(id)) {
			throw std::range_error("HalfEdge out of range: " + id);
		}
		return edges[id];
	}

	inline Face getFace(Ref &id) {
		return getFace(id.ref);
	}

	inline Face getFace(int id) {
		if (!_isFaceInRange(id)) {
			throw std::range_error("Face out of range: " + id);
		}
		return faces[id];
	}

	inline void updateVertex(Vertex &v) {
		if (!_isVertexInRange(v)) {
			throw std::range_error("Vertex out of range: " + v.id.ref);
		}
		vertices[v.id.ref] = v;
	}

	inline void updateHalfEdge(HalfEdge &e) {
		if (!_isHalfEdgeInRange(e)) {
			throw std::range_error("HalfEdge out of range: " + e.id.ref);
		}
		edges[e.id.ref] = e;
	}

	inline void updateFace(Face &f) {
		if (!_isFaceInRange(f)) {
			throw std::range_error("Face out of range: " + f.id.ref);
		}
		faces[f.id.ref] = f;
	}

	inline Vertex origin(HalfEdge h) {
		if (!_isVertexInRange(h.origin)) {
			throw std::range_error("HalfEdge origin out of range: " + h.origin.ref);
		}
		return vertices[h.origin.ref];
	}

	inline HalfEdge twin(HalfEdge h) {
		if (!_isHalfEdgeInRange(h.twin)) {
			throw std::range_error("HalfEdge twin out of range: " + h.twin.ref);
		}
		return edges[h.twin.ref];
	}

	inline Face incidentFace(HalfEdge h) {
		if (!_isFaceInRange(h.incidentFace)) {
			throw std::range_error("HalfEdge incident face out of range: " + h.incidentFace.ref);
		}
		return faces[h.incidentFace.ref];
	}

	inline HalfEdge next(HalfEdge h) {
		if (!_isHalfEdgeInRange(h.next)) {
			throw std::range_error("HalfEdge next of range: " + h.next.ref);
		}
		return edges[h.next.ref];
	}

	inline HalfEdge prev(HalfEdge h) {
		if (!_isHalfEdgeInRange(h.prev)) {
			throw std::range_error("HalfEdge prev of range: " + h.prev.ref);
		}
		return edges[h.prev.ref];
	}

	inline HalfEdge outerComponent(Face f) {
		if (!_isHalfEdgeInRange(f.outerComponent)) {
			throw std::range_error("HalfEdge outer component out of range: " + f.outerComponent.ref);
		}
		return edges[f.outerComponent.ref];
	}

	inline HalfEdge incidentEdge(Vertex v) {
		if (!_isHalfEdgeInRange(v.incidentEdge)) {
			throw std::range_error("Vertex incident edge out of range: " + v.incidentEdge.ref);
		}
		return edges[v.incidentEdge.ref];
	}

	inline bool _isVertexInRange(Vertex &v) {
		return v.id.ref >= 0 && v.id.ref < (int)vertices.size();
	}

	inline bool _isVertexInRange(Ref &id) {
		return id.ref >= 0 && id.ref < (int)vertices.size();
	}

	inline bool _isVertexInRange(int id) {
		return id >= 0 && id < (int)vertices.size();
	}

	inline bool _isHalfEdgeInRange(HalfEdge &h) {
		return h.id.ref >= 0 && h.id.ref < (int)edges.size();
	}

	inline bool _isHalfEdgeInRange(Ref &id) {
		return id.ref >= 0 && id.ref < (int)edges.size();
	}

	inline bool _isHalfEdgeInRange(int id) {
		return id >= 0 && id < (int)edges.size();
	}

	inline bool _isFaceInRange(Face &f) {
		return f.id.ref >= 0 && f.id.ref < (int)faces.size();
	}

	inline bool _isFaceInRange(Ref &id) {
		return id.ref >= 0 && id.ref < (int)faces.size();
	}

	inline bool _isFaceInRange(int id) {
		return id >= 0 && id < (int)faces.size();
	}

	void getOuterComponents(Face f, std::vector<HalfEdge> &edges);
	void getOuterComponents(Face f, std::vector<Ref> &edges);
	void getIncidentEdges(Vertex v, std::vector<HalfEdge> &edges);
	void getIncidentEdges(Vertex v, std::vector<Ref> &edges);
	void getIncidentFaces(Vertex v, std::vector<Face> &faces);
	void getIncidentFaces(Vertex v, std::vector<Ref> &faces);
	bool isBoundary(HalfEdge h);

	std::vector<Vertex> vertices;
	std::vector<HalfEdge> edges;
	std::vector<Face> faces;
};

}

#endif