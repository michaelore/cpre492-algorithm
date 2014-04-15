
#ifndef SOLUTION_H
#define SOLUTION_H

#include <vector>
#include <utility>
#include <CGAL/Object.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include "Utility.h"

namespace CGAL {
    class Solution {

    public:
        Point_3<K> center;
        double rad, eps;

        Solution(Point_3<K> pcenter, double prad, double peps) :
            center(pcenter), rad(prad), eps(peps) {}

        void project_and_display(std::ostream &output) {
            Point_2<K> coordinates = spherical(center);
            output << coordinates.x() << "\t" << coordinates.y() << "\t";
            output << rad << "\t" << eps << "\t";
            output << std::endl;
        }

        /*
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
        */

    };
}
#endif // SOLUTION_H
