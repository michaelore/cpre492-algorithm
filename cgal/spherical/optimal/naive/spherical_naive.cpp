#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Circular_arc_point_3.h>
#include <CGAL/Circular_arc_3.h>
#include <CGAL/Combination_enumerator.h>

#include <vector>
#include <utility>
#include <cmath>

#define TAU 6.283185

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_spherical_kernel_3 S;
typedef CGAL::Point_2<K>             Point_2;
typedef CGAL::Point_3<S>             Point_3;
typedef CGAL::Vector_3<S>            Vector_3;
typedef CGAL::Plane_3<S>             Plane_3;
typedef CGAL::Sphere_3<S>            Sphere_3;
typedef CGAL::Circle_3<S>            Circle_3;
typedef CGAL::Circular_arc_3<S> Arc_3;
typedef CGAL::Circular_arc_point_3<S> Arc_Point_3;

using namespace std;

Sphere_3 UNIT_SPHERE(CGAL::ORIGIN, 1);

Point_3 cartesian(Point_2 point2) {
    K::FT lat = point2.x()*TAU/360;
    K::FT lon = point2.y()*TAU/360;
    K::FT x = cos(lon)*sin(lat);
    K::FT y = sin(lon)*sin(lat);
    K::FT z = cos(lat);
    Point_3 result(x, y, z);
    return result;
}

Point_2 spherical(Arc_Point_3 ap) {
    double x = CGAL::to_double(ap.x());
    double y = CGAL::to_double(ap.y());
    double z = CGAL::to_double(ap.z());
    double lon = atan2(y, x)*360/TAU;
    double lat = acos(z)*360/TAU;
    Point_2 result(lat, lon);
    return result;
}

Point_2 spherical(Point_3 p) {
    double x = CGAL::to_double(p.x());
    double y = CGAL::to_double(p.y());
    double z = CGAL::to_double(p.z());
    double r = sqrt(x*x+y*y+z*z);
    double lon = atan2(y, x)*360/TAU;
    double lat = acos(z/r)*360/TAU;
    Point_2 result(lat, lon);
    return result;
}

void print_circumcenters(Point_3 a, Point_3 b, Point_3 c) {
    Plane_3 plane_ab = CGAL::bisector(a, b);
    Plane_3 plane_bc = CGAL::bisector(b, c);
    vector<CGAL::Object> circumcenters;
    CGAL::intersection(UNIT_SPHERE, plane_ab, plane_bc, back_inserter(circumcenters));
    for (int i = 0; i < circumcenters.size(); i++) {
        Arc_Point_3 ap = CGAL::object_cast<pair<Arc_Point_3, unsigned> >(circumcenters[i]).first;
        cout << "\t" << spherical(ap) << endl;
    }
    cout << "\t" << spherical(CGAL::circumcenter(a, b, c)) << endl;
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
    //} else {
    } else if (S().has_on_3_object()(arc_two, point)) {
        return arc_one;
    } else {
        cout << "Something is very wrong!1" << endl;
        return arc_one;
    }
}

int main() {
    istream_iterator<Point_2> iend;
    vector<Point_3> coordinates;
    for (istream_iterator<Point_2> it(cin); it != iend; it++) {
        coordinates.push_back(cartesian(*it));
    }

    vector<vector<Point_2> > circlesets(coordinates.size());
    CGAL::Combination_enumerator<vector<Point_3>::iterator> set3(3, coordinates.begin(), coordinates.end());
    for (; !set3.finished(); set3++) {
        cout << "Combination: {";
        for (int i = 0; i < 3; i++) {
            cout << spherical(*set3[i]) << ", ";
        }
        cout << "}" << endl;
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
            cout << "\tLocal max with " << i << " point(s) removed: " << (escapable ? "No" : "Yes") << endl;
            //cout << "\t" << escapable << endl;
        }
        cout << "\tLocal max with " << zeroes.size() << " point(s) removed: " << (true ? "No" : "Yes") << endl;
        cout << "\tLesser point count:  " << min(positives, negatives) << endl;
        cout << "\tGreater point count: " << max(positives, negatives) << endl;
        
        //print_circumcenters(*set3[0], *set3[1], *set3[2]);
    }

    return 0;
}
