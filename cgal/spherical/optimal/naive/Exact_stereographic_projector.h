
#ifndef EAXCT_STEREOGRAPHIC_PROJECTOR_H
#define EAXCT_STEREOGRAPHIC_PROJECTOR_H

#include <CGAL/Object.h>

namespace CGAL {
    template <class K>
    class Exact_stereographic_projector {
        typedef typename K::Point_2  Point_2;
        typedef typename K::Point_3  Point_3;
        typedef typename K::Line_3   Line_3;
        typedef typename K::Plane_3  Plane_3;
        typedef typename K::Sphere_3 Sphere_3;

    public:
        Exact_stereographic_projector(Sphere_3 &s, Point_3 &p) :
            sphere(s) {
                antipodal = sphere.center()-(p-sphere.center());
                plane = Plane_3(p, sphere.center()-p);
            }

        Point_2 operator()(Point_3 &p3) { return projection(p3); }

    private:
        Point_3 antipodal;
        Plane_3 plane;
        Sphere_3 sphere;

        Point_2 projection(Point_3 &p3) {
            Line_3 line(antipodal, p3);
            Object o_projection = intersection(line, plane);
            Point_2 projection = plane.to_2d(object_cast<Point_3>(o_projection));
            Point_2 adjusted(-projection.x(), -projection.y());
            return adjusted;
        }
    };
}
#endif // EAXCT_STEREOGRAPHIC_PROJECTOR_H
