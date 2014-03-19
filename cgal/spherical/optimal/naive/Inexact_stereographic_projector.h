
#ifndef INEAXCT_STEREOGRAPHIC_PROJECTOR_H
#define INEAXCT_STEREOGRAPHIC_PROJECTOR_H

#define TAU 6.2831853071

#include <CGAL/Object.h>

namespace CGAL {
    template <class K>
    class Inexact_stereographic_projector {

    public:
        Inexact_stereographic_projector(double lat, double lon, double radius) {
            phi1 = lat*TAU/360;
            lam0 = lon*TAU/360;
            R = radius;
        }

        template <class P3>
        Point_2<K> operator()(P3 &p3) { return projection(p3); }

        template <class P3, class P3E>
        void project_and_display(P3E center, P3 p1, P3 p2, P3 p3) {
            Point_2<K> points2[3];
            points2[0] = projection(p1);
            points2[1] = projection(p2);
            points2[2] = projection(p3);
            Circle_2<K> circle(points2[0], points2[1], points2[2]);
            Point_2<K> true_center = projection(center);
            std::cout << to_double(true_center.x()) << "\t" << to_double(true_center.y()) << "\t";
            std::cout << to_double(circle.center().x()) << "\t" << to_double(circle.center().y()) << "\t";
            std::cout << sqrt(to_double((circle.center()-points2[0]).squared_length()));
            std::cout << std::endl;
        }

        template <class P3, class P3E>
        void project_and_display(P3E center, P3 p1, P3 p2) {
            Point_2<K> points2[2];
            points2[0] = projection(p1);
            points2[1] = projection(p2);
            Circle_2<K> circle(points2[0], points2[1]);
            Point_2<K> true_center = projection(center);
            std::cout << to_double(true_center.x()) << "\t" << to_double(true_center.y()) << "\t";
            std::cout << to_double(circle.center().x()) << "\t" << to_double(circle.center().y()) << "\t";
            std::cout << sqrt(to_double((circle.center()-points2[0]).squared_length()));
            std::cout << std::endl;
        }

    private:
        double phi1;
        double lam0;
        double R;

        template <class P3>
        Point_2<K> projection(P3 &p3) {
            Point_2<K> sphere_coord = spherical(p3);
            double phi = sphere_coord.x()*TAU/360;
            double lam = sphere_coord.y()*TAU/360;
            double k = 2*R/(1 + sin(phi1)*sin(phi)+cos(phi1)*cos(phi)*cos(lam-lam0));
            double x = k*cos(phi)*sin(lam-lam0);
            double y = k*(cos(phi1)*sin(phi)-sin(phi1)*cos(phi)*cos(lam-lam0));
            Point_2<K> result(x, y);
            return result;
        }

        template <class P3>
        Point_2<K> spherical(P3 &p) {
            double x = CGAL::to_double(p.x());
            double y = CGAL::to_double(p.y());
            double z = CGAL::to_double(p.z());
            double r = sqrt(x*x+y*y+z*z);
            double lon = atan2(y, x)*360/TAU;
            double lat = 90-(acos(z/r)*360/TAU);
            Point_2<K> result(lat, lon);
            return result;
        }
    };
}
#endif // INEAXCT_STEREOGRAPHIC_PROJECTOR_H
