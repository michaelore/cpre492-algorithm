#include <iostream>
#include <fstream>

#include <CGAL/Combination_enumerator.h>

#include "Solution.h"
#include "Utility.h"

#include <vector>
#include <utility>
#include <cmath>

//#define USE_STEREOGRAPHIC

using namespace CGAL;

bool compare_radii(Solution a, Solution b) {
    return a.rad > b.rad;
}

class Exact_max {
private:
    int k;
    int ncircs;
    int nrefines;
    std::istream *input;
    std::ostream *output;
    std::ostream *error;

    std::vector<Point_3<K> > coordinates;
    std::vector<Solution> solutions;
    Point_2<K> proj_point;
public:
    Exact_max(int pk, int pncircs, int pnrefines, std::istream &pinput, std::ostream &poutput, std::ostream &perror) :
            k(pk), ncircs(pncircs), nrefines(pnrefines) {
        input = &pinput;
        output = &poutput;
        error = &perror;
        std::istream_iterator<Point_2<K> > iend;
        for (std::istream_iterator<Point_2<K> > it(*input); it != iend; it++) {
            coordinates.push_back(cartesian(*it));
            proj_point = *it;
        }
    }

    void exec_output() {
        Combination_enumerator<std::vector<Point_3<K> >::iterator> set3(3, coordinates.begin(), coordinates.end());
        for (; !set3.finished(); set3++) {
            Plane_3<K> divide(*set3[0], *set3[1], *set3[2]);
            if (divide.is_degenerate()) {
                *error << "WARNING: Skipping combination of 3 because the dividing plane is degenerate" << std::endl;
                continue;
            }
            Plane_3<K> bisect1 = bisecting_plane(*set3[0], *set3[1]);
            Plane_3<K> bisect2 = bisecting_plane(*set3[1], *set3[2]);
            Line_3<K> central;
            Object intersect = intersection(bisect1, bisect2);
            if (!assign(central, intersect)) {
                *error << "WARNING: Skipping combination of 3 because of trouble finding the circumcenters" << std::endl;
                continue;
            }
            Point_3<K> center0 = normalize(ORIGIN+central.to_vector());
            Point_3<K> pos_center;
            Point_3<K> neg_center;
            if (divide.has_on_positive_side(center0)) {
                pos_center = center0;
                neg_center = reverse(center0);
            } else {
                pos_center = reverse(center0);
                neg_center = center0;
            }
            int positives = 0;
            int negatives = 0;
            for (int i = 0; i < coordinates.size(); i++) {
                Point_3<K> p = coordinates[i];
                if (divide.has_on_positive_side(p)) {
                    positives++;
                } else if (divide.has_on_negative_side(p)) {
                    negatives++;
                } else {
                    *error << "DEBUG: Coordinate directly on dividing plane in combination of 3" << std::endl;
                }
            }
            if (positives == k) {
                Solution pos_sol(pos_center, great_circle_dist(pos_center, *set3[0]), EPS);
                solutions.push_back(pos_sol);
            }
            if (negatives == k) {
                Solution neg_sol(neg_center, great_circle_dist(neg_center, *set3[0]), EPS);
                solutions.push_back(neg_sol);
            }
        }
        Combination_enumerator<std::vector<Point_3<K> >::iterator> set2(2, coordinates.begin(), coordinates.end());
        for (; !set2.finished(); set2++) {
            Point_3<K> mid = midpoint(*set2[0], *set2[1]);
            if ((mid-ORIGIN).squared_length() < EPS) {
                *error << "WARNING: Skipping combination of 2 because they are ~antipodal" << std::endl;
            }
            Plane_3<K> divide(*set3[0], mid-ORIGIN);
            if (divide.is_degenerate()) {
                *error << "WARNING: Skipping combination of 2 because the dividing plane is degenerate" << std::endl;
                continue;
            }
            Point_3<K> center0 = normalize(mid);
            Point_3<K> pos_center;
            Point_3<K> neg_center;
            if (divide.has_on_positive_side(center0)) {
                pos_center = center0;
                neg_center = reverse(center0);
            } else {
                pos_center = reverse(center0);
                neg_center = center0;
            }
            int positives = 0;
            int negatives = 0;
            for (int i = 0; i < coordinates.size(); i++) {
                Point_3<K> p = coordinates[i];
                if (divide.has_on_positive_side(p)) {
                    positives++;
                } else if (divide.has_on_negative_side(p)) {
                    negatives++;
                } else {
                    *error << "DEBUG: Coordinate directly on dividing plane in combination of 2" << std::endl;
                }
            }
            if (positives == k) {
                Solution pos_sol(pos_center, great_circle_dist(pos_center, *set2[0]), EPS);
                solutions.push_back(pos_sol);
            }
            if (negatives == k) {
                Solution neg_sol(neg_center, great_circle_dist(neg_center, *set2[0]), EPS);
                solutions.push_back(neg_sol);
            }
        }
        std::sort(solutions.begin(), solutions.end(), compare_radii);
        int nsolutions = solutions.size();
        int limit = std::min(nsolutions, ncircs);
        for (int i = 0; i < limit; i++) {
            #ifdef USE_STEREOGRAPHIC
            //solutions[i].project_and_display(stereo);
            #else
            solutions[i].project_and_display(*output);
            #endif
        }
    }
};

namespace CGAL {
    void run_exact_max(int k, int ncircs, int nrefines, std::istream &input, std::ostream &output, std::ostream &error) {
        Exact_max em(k, ncircs, nrefines, input, output, error);
        em.exec_output();
    }
}
