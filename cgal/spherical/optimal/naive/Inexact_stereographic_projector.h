
#ifndef INEAXCT_STEREOGRAPHIC_PROJECTOR_H
#define INEAXCT_STEREOGRAPHIC_PROJECTOR_H

#define TAU 6.2831853071

#include <CGAL/Object.h>

namespace CGAL {
    template <class K>
    class Inexact_stereographic_projector {
        typedef typename K::Point_2  Point_2;

    public:
        Inexact_stereographic_projector(double lat, double lon, double radius) {
            phi1 = lat*TAU/360;
            lam0 = lon*TAU/360;
            R = radius;
        }

        template <class P3>
        Point_2 operator()(P3 &p3) { return projection(p3); }

    private:
        double phi1;
        double lam0;
        double R;

        template <class P3>
        Point_2 projection(P3 &p3) {
            Point_2 sphere_coord = spherical(p3);
            double phi = sphere_coord.x()*TAU/360;
            double lam = sphere_coord.y()*TAU/360;
            double k = 2*R/(1 + sin(phi1)*sin(phi)+cos(phi1)*cos(phi)*cos(lam-lam0));
            double x = k*cos(phi)*sin(lam-lam0);
            double y = k*(cos(phi1)*sin(phi)-sin(phi1)*cos(phi)*cos(lam-lam0));
            Point_2 result(x, y);
            return result;
        }

        template <class P3>
        Point_2 spherical(P3 &p) {
            double x = CGAL::to_double(p.x());
            double y = CGAL::to_double(p.y());
            double z = CGAL::to_double(p.z());
            double r = sqrt(x*x+y*y+z*z);
            double lon = atan2(y, x)*360/TAU;
            double lat = 90-(acos(z/r)*360/TAU);
            Point_2 result(lat, lon);
            return result;
        }
    };
}
#endif // INEAXCT_STEREOGRAPHIC_PROJECTOR_H
