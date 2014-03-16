
#ifndef STEREOGRAPHIC_PROJECTOR_H
#define STEREOGRAPHIC_PROJECTOR_H

#include <CGAL/Object.h>

namespace CGAL {
    template <class K>
    class Stereographic_projector {
        typedef typename K::Point_2  Point_2;
        typedef typename K::Point_3  Point_3;
        typedef typename K::Line_3   Line_3;
        typedef typename K::Plane_3  Plane_3;
        typedef typename K::Sphere_3 Sphere_3;

    public:
        Stereographic_projector(Sphere_3 &s, Point_3 &p) :
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
            Object projection = intersection(line, plane);
            return plane.to_2d(object_cast<Point_3>(projection));
        }
    };
}
#endif // STEREOGRAPHIC_PROJECTOR_H
