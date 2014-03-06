#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Circular_arc_point_3.h>
#include <CGAL/Circular_arc_3.h>
#include <CGAL/Combination_enumerator.h>

#include <Stereographic_projector.h>

#include <vector>
#include <utility>
#include <cmath>

#define TAU 6.283185

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
typedef CGAL::Stereographic_projector<S> Stereo_Projector;

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

template <class P3>
CGAL::Point_2<K> spherical(P3 p) {
    double x = CGAL::to_double(p.x());
    double y = CGAL::to_double(p.y());
    double z = CGAL::to_double(p.z());
    double r = sqrt(x*x+y*y+z*z);
    double lon = atan2(y, x)*360/TAU;
    double lat = acos(z/r)*360/TAU;
    CGAL::Point_2<K> result(lat, lon);
    return result;
}

CGAL::Point_2<K> projection(Point_3 p3) {
    CGAL::Point_2<K> sphere_coord = spherical(p3);
    double phi = sphere_coord.x()*TAU/360;
    double lam = sphere_coord.y()*TAU/360;
    double lam0 = 0;
    double phi1 = 0;
    double R = 1;
    double k = 2*R/(1 + sin(phi1)*sin(phi)+cos(phi1)*cos(phi)*cos(lam-lam0));
    double x = k*cos(phi)*sin(lam-lam0);
    double y = k*(cos(phi1)*sin(phi)-sin(phi1)*cos(phi)*cos(lam-lam0));
    CGAL::Point_2<K> result(x, y);
    return result;
}

Arc_3 get_opposing_arc(Circle_3 circle, Point_3 point) {
    Plane_3 orthogonal_plane(circle.center(), point-circle.center());

    vector<CGAL::Object> orthogonal_points;
    CGAL::intersection(circle, orthogonal_plane, back_inserter(orthogonal_points));
    Arc_Point_3 orthogonal_one = CGAL::object_cast<pair<Arc_Point_3, unsigned> >(orthogonal_points[0]).first;
    Arc_Point_3 orthogonal_two = CGAL::object_cast<pair<Arc_Point_3, unsigned> >(orthogonal_points[1]).first;

    Arc_3 arc_one(circle, orthogonal_one, orthogonal_two);
    Arc_3 arc_two(circle, orthogonal_two, orthogonal_one);
    if (S().has_on_3_object()(arc_one, point)) {
        return arc_two;
    } else if (S().has_on_3_object()(arc_two, point)) {
        return arc_one;
    } else {
        cout << "Something is very wrong!1" << endl;
        return arc_one;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        return 1;
    }
    int k = atoi(argv[1]);
    istream_iterator<CGAL::Point_2<K> > iend;
    vector<Point_3> coordinates;
    Point_3 proj_point;
    ifstream ifs(argv[2], ifstream::in);
    for (istream_iterator<CGAL::Point_2<K> > it(ifs); it != iend; it++) {
        coordinates.push_back(cartesian(*it));
        proj_point = cartesian(*it);
    }
    ifs.close();
    Stereo_Projector stereo(UNIT_SPHERE, proj_point);

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
                bool empty_range = false;
                Arc_3 full_range(circle);
                CGAL::Object range = CGAL::make_object(full_range);
                for (int j = 0; j < retain; j++) {
                    Arc_3 opp_arc = get_opposing_arc(circle, *selected[j]);
                    vector<CGAL::Object> restricted;
                    if (const Arc_3* casted_range = CGAL::object_cast<Arc_3>(&range)) {
                        CGAL::intersection(*casted_range, opp_arc, back_inserter(restricted));
                        if (restricted.size() > 1) {
                            cout << "Something is very wrong!2" << endl;
                        } else if (restricted.size() == 1) {
                            range = restricted[0];
                        } else {
                            empty_range = true;
                            break;
                        }
                    } else if (const pair<Arc_Point_3, unsigned>* casted_range_with_mult = CGAL::object_cast<pair<Arc_Point_3, unsigned> >(&range)) {
                        if (!S().has_on_3_object()(opp_arc, casted_range_with_mult->first)) {
                            empty_range = true;
                            break;
                        }
                    } else {
                        //
                    }
                }
                if (!empty_range) {
                    escapable = true;
                    break;
                }
            }
            if (!escapable && (positives + i == k || negatives + i == k)) {
                CGAL::Point_2<S> points2[3];
                points2[0] = stereo(*set3[0]);
                points2[1] = stereo(*set3[1]);
                points2[2] = stereo(*set3[2]);
                for (int i = 0; i < 3; i++) {
                    cout << CGAL::to_double(points2[i].x()) << "\t" << CGAL::to_double(points2[i].y()) << "\t";
                }
                CGAL::Circle_2<S> circle(points2[0], points2[1], points2[2]);
                cout << CGAL::to_double(circle.center().x()) << "\t" << CGAL::to_double(circle.center().y()) << "\t";
                cout << sqrt(CGAL::to_double((circle.center()-points2[0]).squared_length()));
                cout << endl;
            }
        }
    }

    return 0;
}
