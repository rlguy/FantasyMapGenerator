#include "dcel.h"

void dcel::DCEL::getOuterComponents(Face f, std::vector<HalfEdge> &edges) {
	HalfEdge h = outerComponent(f);
	Ref startid = h.id;

	do {
		edges.push_back(h);
		h = next(h);
	} while (h.id != startid);
}

void dcel::DCEL::getOuterComponents(Face f, std::vector<Ref> &edges) {
	HalfEdge h = outerComponent(f);
	Ref startid = h.id;

	do {
		edges.push_back(h.id);
		h = next(h);
	} while (h.id != startid);
}

void dcel::DCEL::getIncidentEdges(Vertex v, std::vector<HalfEdge> &edges) {
	HalfEdge h = incidentEdge(v);
	Ref startid = h.id;

	do {
		edges.push_back(h);
		h = next(twin(h));
	} while (h.id != startid);
}

void dcel::DCEL::getIncidentEdges(Vertex v, std::vector<Ref> &edges) {
	HalfEdge h = incidentEdge(v);
	Ref startid = h.id;

	do {
		edges.push_back(h.id);
		h = next(twin(h));
	} while (h.id != startid);
}

void dcel::DCEL::getIncidentFaces(Vertex v, std::vector<Face> &faces) {
	HalfEdge h = incidentEdge(v);
	Ref startid = h.id;

	do {
		if (!isBoundary(h)) {
			faces.push_back(incidentFace(h));
		}
		h = next(twin(h));
	} while (h.id != startid);
}

void dcel::DCEL::getIncidentFaces(Vertex v, std::vector<Ref> &faces) {
	HalfEdge h = incidentEdge(v);
	Ref startid = h.id;

	do {
		if (!isBoundary(h)) {
			faces.push_back(incidentFace(h).id);
		}
		h = next(twin(h));
	} while (h.id != startid);
}

bool dcel::DCEL::isBoundary(HalfEdge h) {
	return h.incidentFace.ref == -1;
}
