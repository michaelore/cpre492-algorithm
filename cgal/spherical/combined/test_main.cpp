#include <fstream>
#include <iostream>
#include <utility>
#include <random>
#include "Exact_max.h"
#include "Grid_sampling.h"
#include "Utility.h"

//#define USE_STEREOGRAPHIC

using namespace std;

default_random_engine gen;

bool test_cartesian_inverse() {
    uniform_real_distribution<double> latgen(-90, 90);
    uniform_real_distribution<double> longen(-180, 180);
    bool success = true;
    for (int i = 0; i < 10000; i++) {
        CGAL::Point_2<K> test_point(latgen(gen), longen(gen));
        CGAL::Point_2<K> result = spherical(cartesian(test_point));
        if ((test_point-result).squared_length() > EPSILON) {
            success = false;
        }
    }
    return success;
}

bool test_stereo_inverse() {
    uniform_real_distribution<double> latgen(-90, 90);
    uniform_real_distribution<double> longen(-180, 180);
    bool success = true;
    for (int i = 0; i < 10000; i++) {
        CGAL::Point_2<K> test_spherical(latgen(gen), longen(gen));
        CGAL::Point_3<K> test_point = cartesian(test_spherical);
        double clat = latgen(gen);
        double clon = longen(gen);
        CGAL::Point_3<K> result = inv_stereographic(stereographic(test_point, clat, clon, RADIUS), clat, clon, RADIUS);
        if ((test_point-result).squared_length() > EPSILON) {
            success = false;
            cout << test_point << "\t" << result << endl;
        } else {
        }
    }
    return success;
}

int main(int argc, char *argv[]) {
    cout << "test_cartesian_inverse(): " << (test_cartesian_inverse() ? "TRUE" : "FALSE")  << endl;
    cout << "test_stereo_inverse(): " << (test_stereo_inverse() ? "TRUE" : "FALSE")  << endl;
    return 0;
}
