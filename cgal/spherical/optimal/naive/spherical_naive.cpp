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
#define RADIUS 6371009

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_spherical_kernel_3 S;
typedef CGAL::Cartesian<S::Root_of_2> E;
typedef CGAL::Exact_stereographic_projector<S> Exact_Stereo_Projector;
typedef CGAL::Inexact_stereographic_projector<K> Inexact_Stereo_Projector;

using namespace CGAL;

Sphere_3<S> UNIT_SPHERE(ORIGIN, 1);

Point_3<S> cartesian(Point_2<K> point2) {
    K::FT lat = point2.x()*TAU/360;
    K::FT lon = point2.y()*TAU/360;
    K::FT x = cos(lon)*sin(lat);
    K::FT y = sin(lon)*sin(lat);
    K::FT z = cos(lat);
    Point_3<S> result(x, y, z);
    return result;
}

Plane_3<S> bisecting_plane(Point_3<S> p1, Point_3<S> p2) {
    Plane_3<S> result(midpoint(p1, p2), p2-p1);
    return result;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: sphericalnaive <k> <filename>" << std::endl;
        return 1;
    }
    int k = atoi(argv[1]);
    std::istream_iterator<Point_2<K> > iend;
    std::vector<Point_3<S> > coordinates;
    Point_2<K> proj_point;
    std::ifstream ifs(argv[2], std::ifstream::in);
    for (std::istream_iterator<Point_2<K> > it(ifs); it != iend; it++) {
        coordinates.push_back(cartesian(*it));
        proj_point = *it;
    }
    ifs.close();
    //Exact_Stereo_Projector stereo(UNIT_SPHERE, cartesian(proj_point));
    Inexact_Stereo_Projector stereo(proj_point.x(), proj_point.y(), RADIUS);

    Combination_enumerator<std::vector<Point_3<S> >::iterator> set3(3, coordinates.begin(), coordinates.end());
    for (; !set3.finished(); set3++) {
        Plane_3<S> divide(*set3[0], *set3[1], *set3[2]);
        Plane_3<E> e_divide(E::FT(divide.a()), E::FT(divide.b()), E::FT(divide.c()), E::FT(divide.d()));
        Circle_3<S> circle(*set3[0], *set3[1], *set3[2]);
        int positives = 0;
        int negatives = 0;
        std::vector<Point_3<S> > zeroes;
        for (int i = 0; i < coordinates.size(); i++) {
            Point_3<S> p = coordinates[i];
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
            Combination_enumerator<std::vector<Point_3<S> >::iterator> selected(retain, zeroes.begin(), zeroes.end());
            for (; !selected.finished(); selected++) {
                Angle_range<S> range(circle);
                for (int j = 0; j < retain; j++) {
                    range.restrict_with_point(*selected[j]);
                }
                if (!range.empty()) {
                    escapable = true;
                }
            }
            if (escapable) {
                continue;
            }
            std::vector<Object> centers;
            Plane_3<S> bisect1 = bisecting_plane(*set3[0], *set3[1]);
            Plane_3<S> bisect2 = bisecting_plane(*set3[1], *set3[2]);
            intersection(UNIT_SPHERE, bisect1, bisect2, back_inserter(centers));
            Point_3<E> pos_center;
            Point_3<E> neg_center;
            for (int j = 0; j < centers.size(); j++) {
                std::pair<Circular_arc_point_3<S>, unsigned> app;
                assign(app, centers[j]);
                Point_3<E> center(E::FT(app.first.x()), E::FT(app.first.y()), E::FT(app.first.z()));
                if (e_divide.has_on_positive_side(center)) {
                    pos_center = center;
                } else if (e_divide.has_on_negative_side(center)) {
                    neg_center = center;
                } else {
                    //
                }
            }
            if (positives + i == k) {
                Point_2<K> points2[3];
                points2[0] = stereo(*set3[0]);
                points2[1] = stereo(*set3[1]);
                points2[2] = stereo(*set3[2]);
                Circle_2<K> circle(points2[0], points2[1], points2[2]);
                Point_2<K> true_center = stereo(pos_center);
                std::cout << to_double(true_center.x()) << "\t" << to_double(true_center.y()) << "\t";
                std::cout << to_double(circle.center().x()) << "\t" << to_double(circle.center().y()) << "\t";
                std::cout << sqrt(to_double((circle.center()-points2[0]).squared_length()));
                std::cout << std::endl;
            }
            if (negatives + i == k) {
                Point_2<K> points2[3];
                points2[0] = stereo(*set3[0]);
                points2[1] = stereo(*set3[1]);
                points2[2] = stereo(*set3[2]);
                Circle_2<K> circle(points2[0], points2[1], points2[2]);
                Point_2<K> true_center = stereo(neg_center);
                std::cout << to_double(true_center.x()) << "\t" << to_double(true_center.y()) << "\t";
                std::cout << to_double(circle.center().x()) << "\t" << to_double(circle.center().y()) << "\t";
                std::cout << sqrt(to_double((circle.center()-points2[0]).squared_length()));
                std::cout << std::endl;
            }
        }
    }
    return 0;
}
