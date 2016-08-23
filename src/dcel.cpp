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

void dcel::DCEL::removeEdge(HalfEdge h) {
	if (h.origin.ref == -1) {
		return;
	}

	HalfEdge eij = h;
	HalfEdge eji = twin(h);
	if (isBoundary(eji)) {
		eij = eji;
		eji = h;
	}

	Ref faceid = eij.incidentFace;

	Vertex vi = origin(eij);
	Vertex vj = origin(eji);

	HalfEdge eri = prev(eij);
	HalfEdge ein = next(eji);
	HalfEdge emj = prev(eji);
	HalfEdge ejk = next(eij);

	std::vector<Ref> loopedges;
	if (!isBoundary(eji)) {
		getOuterComponents(incidentFace(eji), loopedges);
	}

	// update components
	if (vi.incidentEdge == eij.id) {
		vi.incidentEdge = ein.id;
		updateVertex(vi);
	}

	if (vj.incidentEdge == eji.id) {
		vj.incidentEdge = ejk.id;
		updateVertex(vj);
	}

	eri.next = ein.id;
	updateHalfEdge(eri);

	ein.prev = eri.id;
	updateHalfEdge(ein);

	emj.next = ejk.id;
	updateHalfEdge(emj);

	ejk.prev = emj.id;
	updateHalfEdge(ejk);

	if (!isBoundary(eij)) {
		Face f = incidentFace(eij);
		if (f.outerComponent == eij.id) {
			f.outerComponent = ejk.id;
			updateFace(f);
		}
	}

	std::cout << "start " << faceid.ref << std::endl;
	for (unsigned int i = 0; i < loopedges.size(); i++) {
		HalfEdge h = getHalfEdge(loopedges[i]);
		if (h.id == eij.id || h.id == eji.id) {
			continue;
		}

		h.incidentFace = faceid;
		if (incidentFace(h).outerComponent.ref == -1) {
			std::cout << "1 Error: " << faceid.ref << i << std::endl;
		}

		updateHalfEdge(h);
	}
	std::cout << "end " << faceid.ref << std::endl;

	// remove components
	if (!isBoundary(eji)) {
		Face f = incidentFace(eji);
		f.outerComponent = Ref();
		updateFace(f);

		std::cout << "faceid: " << f.id.ref << std::endl;
	}

	eij.origin = Ref();
	eij.twin = Ref();
	eij.incidentFace = Ref();
	eij.next = Ref();
	eij.prev = Ref();
	updateHalfEdge(eij);

	eji.origin = Ref();
	eji.twin = Ref();
	eji.incidentFace = Ref();
	eji.next = Ref();
	eji.prev = Ref();
	updateHalfEdge(eji);

	for (unsigned int i = 0; i < edges.size(); i++) {
		HalfEdge h = edges[i];
		if (h.incidentFace.ref != -1 && incidentFace(h).outerComponent.ref == -1) {
			std::cout << "Error: " << i << std::endl;
		}
	}
}