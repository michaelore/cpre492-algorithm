
#ifndef SOLUTION_H
#define SOLUTION_H

#define TAU 6.2831853071

#include <vector>
#include <utility>
#include <CGAL/Object.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace CGAL {
    template <class K>
    class Solution {
        typedef typename K::Root_of_2 KE;
        typedef typename K::Circle_3 Circle_3;
        typedef typename K::Point_3 Point_3;
        typedef Cartesian<KE> E;
        typedef typename E::Point_3 Point_3E;
        typedef typename K::FT FT;
        typedef typename E::FT EFT;
        typedef CGAL::Exact_predicates_inexact_constructions_kernel D;
        typedef typename D::Point_2 Point_2D;
        typedef typename D::Circle_2 Circle_2D;

    public:
        Circle_3 circle;
        Point_3E center;
        std::vector<Point_3> points;

        Solution(Circle_3 &circ, Point_3E &cent, std::vector<Point_3 > &ps) :
            circle(circ), center(cent), points(ps) {}

        EFT sq_dist() {
            Point_3E extend(EFT(circle.center().x()), EFT(circle.center().y()), EFT(circle.center().z()));
            return (center-extend).squared_length();
        }

        template <class P3>
        Point_2D spherical(P3 &p) {
            double x = to_double(p.x());
            double y = to_double(p.y());
            double z = to_double(p.z());
            double r = sqrt(x*x+y*y+z*z);
            double lon = atan2(y, x)*360/TAU;
            double lat = 90-(acos(z/r)*360/TAU);
            Point_2D result(lat, lon);
            return result;
        }

        void project_and_display() {
            Point_2D true_center = spherical(center);
            std::cout << to_double(true_center.x()) << "\t" << to_double(true_center.y()) << "\t";
            for (int i = 0; i < points.size(); i++) {
                Point_2D point = spherical(points[i]);
                std::cout << to_double(point.x()) << "\t" << to_double(point.y()) << "\t";
            }
            std::cout << std::endl;
        }

        template <class STEREO>
        void project_and_display(STEREO stereo) {
            Circle_2D circle;
            if (points.size() == 2) {
                circle = Circle_2D(stereo(points[0]), stereo(points[1]));
            } else if (points.size() == 3) {
                circle = Circle_2D(stereo(points[0]), stereo(points[1]), stereo(points[2]));
            }
            Point_2D true_center = stereo(center);
            std::cout << to_double(true_center.x()) << "\t" << to_double(true_center.y()) << "\t";
            std::cout << to_double(circle.center().x()) << "\t" << to_double(circle.center().y()) << "\t";
            std::cout << sqrt(to_double((circle.center()-stereo(points[0])).squared_length()));
            std::cout << std::endl;
        }

    };
}
#endif // SOLUTION_H
