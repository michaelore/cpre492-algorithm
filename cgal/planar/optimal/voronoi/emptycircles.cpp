// CprE 491 May14-31 optimal prototype
// Author: Michael Ore

#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <queue>
#include <utility>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Point_2<K>             Point;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;

// This program solves the largest empty circle problem, excluding circles with
// centers on the convex hull edge. Those are harder to find.
int main() {
    std::ifstream in("points.txt");
    std::istream_iterator<Point> begin(in);
    std::istream_iterator<Point> end;

    // Create a Delaunay triangulation from the points in "points.txt"
    Delaunay dt;
    dt.insert(begin, end);

    // Priority queue containing pointers to faces on the triangulation
    // Sorted by decreasing squared radius (same as sorting by radius)
    std::priority_queue<std::pair<double, Delaunay::Face_handle> > faces;

    // Iterate over all faces of the triangulation
    Delaunay::Finite_faces_iterator it;
    for (it = dt.finite_faces_begin(); it != dt.finite_faces_end(); ++it) {
        // dt.dual(...) takes an element of the Delaunay triangulation, and
        // returns the corresponding element of the Voronoi diagram. For faces,
        // dt.dual(...) returns a vertex on the Voronoi diagram aka circumcenter
        Point circumcenter = dt.dual(it);
        // dt.locate(...) returns the face that the given point is inside of.
        Delaunay::Face_handle other_face = dt.locate(circumcenter);
        // If the point is not outside the convex hull
        if (other_face != NULL && !dt.is_infinite(other_face)) {
            // Calculate squared radius from squared distance between the
            // circumcenter and an arbitrary vertex of the face, then push onto
            // the priority queue.
            faces.push(std::make_pair(CGAL::to_double(squared_distance(circumcenter, it->vertex(0)->point())), it));
        }
    }
    // Print the circle centers and radii. The largest circle comes first.
    while (!faces.empty()) {
        std::cout << dt.dual(faces.top().second) << std::endl;
        std::cout << sqrt(faces.top().first) << std::endl;
        std::cout << std::endl;
        faces.pop();
    }
    return 0;
}
