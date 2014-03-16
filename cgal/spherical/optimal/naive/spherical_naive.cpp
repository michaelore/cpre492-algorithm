#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Circular_arc_point_3.h>
#include <CGAL/Circular_arc_3.h>
#include <CGAL/Combination_enumerator.h>

#include "Exact_stereographic_projector.h"
#include "Inexact_stereographic_projector.h"
#include "Angle_range.h"

#include <vector>
#include <utility>
#include <cmath>

#define TAU 6.2831853071

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_spherical_kernel_3 S;
typedef CGAL::Point_3<S>             Point_3;
typedef CGAL::Vector_3<S>            Vector_3;
typedef CGAL::Line_3<S>              Line_3;
typedef CGAL::Plane_3<S>             Plane_3;
typedef CGAL::Sphere_3<S>            Sphere_3;
typedef CGAL::Circle_3<S>            Circle_3;
typedef CGAL::Circular_arc_3<S> Arc_3;
typedef CGAL::Circular_arc_point_3<S> Arc_Point_3;
typedef CGAL::Exact_stereographic_projector<S> Exact_Stereo_Projector;
typedef CGAL::Inexact_stereographic_projector<K> Inexact_Stereo_Projector;

using namespace std;

Sphere_3 UNIT_SPHERE(CGAL::ORIGIN, 1);

Point_3 cartesian(CGAL::Point_2<K> point2) {
    K::FT lat = point2.x()*TAU/360;
    K::FT lon = point2.y()*TAU/360;
    K::FT x = cos(lon)*sin(lat);
    K::FT y = sin(lon)*sin(lat);
    K::FT z = cos(lat);
    Point_3 result(x, y, z);
    return result;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: sphericalnaive <k> <filename>" << endl;
        return 1;
    }
    int k = atoi(argv[1]);
    istream_iterator<CGAL::Point_2<K> > iend;
    vector<Point_3> coordinates;
    CGAL::Point_2<K> proj_point;
    ifstream ifs(argv[2], ifstream::in);
    for (istream_iterator<CGAL::Point_2<K> > it(ifs); it != iend; it++) {
        coordinates.push_back(cartesian(*it));
        proj_point = *it;
    }
    ifs.close();
    //Exact_Stereo_Projector stereo(UNIT_SPHERE, cartesian(proj_point));
    Inexact_Stereo_Projector stereo(proj_point.x(), proj_point.y());

    vector<vector<CGAL::Point_2<K> > > circlesets(coordinates.size());
    CGAL::Combination_enumerator<vector<Point_3>::iterator> set3(3, coordinates.begin(), coordinates.end());
    for (; !set3.finished(); set3++) {
        Plane_3 divide = Plane_3(*set3[0], *set3[1], *set3[2]);
        Circle_3 circle = Circle_3(*set3[0], *set3[1], *set3[2]);
        int positives = 0;
        int negatives = 0;
        vector<Point_3> zeroes;
        for (int i = 0; i < coordinates.size(); i++) {
            Point_3 p = coordinates[i];
            if (divide.has_on_positive_side(p)) {
                positives++;
            } else if (divide.has_on_negative_side(p)) {
                negatives++;
            } else {
                zeroes.push_back(p);
            }
        }
        for (int i = 0; i < zeroes.size(); i++) {
            int retain = zeroes.size()-i;
            bool escapable = false;
            CGAL::Combination_enumerator<vector<Point_3>::iterator> selected(retain, zeroes.begin(), zeroes.end());
            for (; !selected.finished(); selected++) {
                CGAL::Angle_range<S> range(circle);
                for (int j = 0; j < retain; j++) {
                    range.restrict_with_point(*selected[j]);
                }
                if (!range.is_empty()) {
                    escapable = true;
                }
            }
            if (!escapable && (positives + i == k || negatives + i == k)) {
                CGAL::Point_2<K> points2[3];
                points2[0] = stereo(*set3[0]);
                points2[1] = stereo(*set3[1]);
                points2[2] = stereo(*set3[2]);
                for (int i = 0; i < 3; i++) {
                    cout << CGAL::to_double(points2[i].x()) << "\t" << CGAL::to_double(points2[i].y()) << "\t";
                }
                CGAL::Circle_2<K> circle(points2[0], points2[1], points2[2]);
                cout << CGAL::to_double(circle.center().x()) << "\t" << CGAL::to_double(circle.center().y()) << "\t";
                cout << sqrt(CGAL::to_double((circle.center()-points2[0]).squared_length()));
                cout << endl;
            }
        }
    }

    return 0;
}
