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
		if (id < 0 || id >= (int)vertices.size()) {
			throw std::range_error("");
		}
		return vertices[id];
	}

	inline HalfEdge getHalfEdge(Ref &id) {
		return getHalfEdge(id.ref);
	}

	inline HalfEdge getHalfEdge(int id) {
		if (id < 0 || id >= (int)edges.size()) {
			throw std::range_error("");
		}
		return edges[id];
	}

	inline Face getFace(Ref &id) {
		return getFace(id.ref);
	}

	inline Face getFace(int id) {
		if (id < 0 || id >= (int)faces.size()) {
			throw std::range_error("");
		}
		return faces[id];
	}

	inline void updateVertex(Vertex &v) {
		if (v.id.ref < 0 || v.id.ref >= (int)vertices.size()) {
			throw std::range_error("");
		}
		vertices[v.id.ref] = v;
	}

	inline void updateHalfEdge(HalfEdge &e) {
		if (e.id.ref < 0 || e.id.ref >= (int)edges.size()) {
			throw std::range_error("");
		}
		edges[e.id.ref] = e;
	}

	inline void updateFace(Face &f) {
		if (f.id.ref < 0 || f.id.ref >= (int)faces.size()) {
			throw std::range_error("");
		}
		faces[f.id.ref] = f;
	}

	inline Vertex origin(HalfEdge h) {
		if (h.origin.ref < 0 || h.origin.ref >= (int)vertices.size()) {
			throw std::range_error("");
		}
		return vertices[h.origin.ref];
	}

	inline HalfEdge twin(HalfEdge h) {
		if (h.twin.ref < 0 || h.twin.ref >= (int)edges.size()) {
			throw std::range_error("");
		}
		return edges[h.twin.ref];
	}

	inline Face incidentFace(HalfEdge h) {
		if (h.incidentFace.ref < 0 || h.incidentFace.ref >= (int)faces.size()) {
			throw std::range_error("");
		}
		return faces[h.incidentFace.ref];
	}

	inline HalfEdge next(HalfEdge h) {
		if (h.next.ref < 0 || h.next.ref >= (int)edges.size()) {
			throw std::range_error("");
		}
		return edges[h.next.ref];
	}

	inline HalfEdge prev(HalfEdge h) {
		if (h.prev.ref < 0 || h.prev.ref >= (int)edges.size()) {
			throw std::range_error("");
		}
		return edges[h.prev.ref];
	}

	inline HalfEdge outerComponent(Face f) {
		if (f.outerComponent.ref < 0 || f.outerComponent.ref >= (int)edges.size()) {
			throw std::range_error("");
		}
		return edges[f.outerComponent.ref];
	}

	inline std::vector<HalfEdge> innerComponents(Face f) {
		std::vector<HalfEdge> inner;
		for (unsigned int i = 0; i < f.innerComponents.size(); i++) {
			Ref r = f.innerComponents[i];
			if (r.ref < 0 || r.ref >= (int)f.innerComponents.size()) {
				throw std::range_error("");
			}
			inner.push_back(edges[f.innerComponents[i].ref]);
		}
		return inner;
	}

	inline HalfEdge incidentEdge(Vertex v) {
		if (v.incidentEdge.ref < 0 || v.incidentEdge.ref >= (int)edges.size()) {
			throw std::range_error("");
		}
		return edges[v.incidentEdge.ref];
	}

	void getOuterComponents(Face f, std::vector<HalfEdge> &edges);
	void getOuterComponents(Face f, std::vector<Ref> &edges);
	void getIncidentEdges(Vertex v, std::vector<HalfEdge> &edges);
	void getIncidentEdges(Vertex v, std::vector<Ref> &edges);
	void getIncidentFaces(Vertex v, std::vector<Face> &faces);
	void getIncidentFaces(Vertex v, std::vector<Ref> &faces);
	bool isBoundary(HalfEdge h);
	void removeEdge(HalfEdge h);

	std::vector<Vertex> vertices;
	std::vector<HalfEdge> edges;
	std::vector<Face> faces;
};

}

#endif