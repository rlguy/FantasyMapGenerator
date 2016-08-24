# Fantasy Map Generator (Work in Progress)

This program is an implementation of a fantasy map generator based on the methods described in Martin O'Leary's "Generating fantasy map" notes (https://mewo2.com/notes/terrain/).

# Progress

## Generating Irregular Grids

A Poisson Disc Sampler generates a set of random points with the property that no two points are within a set radius of eachother.
![alt tag](http://rlguy.com/map_generation/images/uniform_vs_poisson_sampling.jpg)

The set of points are triangulated in a Delaunay triangulation. The triangulation is stored in a doubly connected edge list (DCEL) data structure.
![alt tag](http://rlguy.com/map_generation/images/uniform_vs_poisson_delaunay.jpg)

The dual of the Delaunay triangulation is computed to produce a Voronoi diagram, which is also stored as a DCEL.
![alt tag](http://rlguy.com/map_generation/images/uniform_vs_poisson_voronoi.jpg)

Each vertex in the Delaunay triangulation becomes a face in the Voronoi diagram, and each triangle in the Delaunay triangulation becomes a vertex in the Voronoi diagram. A triangle is transformed into a vertex by fitting a circle to the three triangle vertices and setting the circle's center as the position of a Voronoi vertex. The following image displays the relationship between a Delaunay triangulation and a Voronoi diagram.

![alt tag](http://rlguy.com/map_generation/images/voronoi_delaunay_overlay.jpg)

The vertices in the Voronoi diagram will be used as the nodes in an irregular grid. Note that each node has exactly three neighbours.