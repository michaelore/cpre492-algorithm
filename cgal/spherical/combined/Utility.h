
#ifndef UTILITY_H
#define UTILITY_H

#define TAU 6.2831853071
#define RADIUS 6371009
#define EPSILON 1e-8

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Sphere_3.h>

namespace CGAL {

    typedef Simple_cartesian<double>                K;

    inline double clmp(double val) {
        if (val > 1) {
            return 1;
        } else if (val < -1) {
            return -1;
        } else {
            return val;
        }
    }

    inline double chord_to_arc(double c) {
        double s = 2*asin(c/2);
        return RADIUS*s;
    }

    inline double great_circle_dist(Point_3<K> p1, Point_3<K> p2) {
        return chord_to_arc(sqrt((p1-p2).squared_length()));
    }

    Sphere_3<K> UNIT_SPHERE(ORIGIN, 1);

    Vector_3<K> normalize(const Vector_3<K> &v) {
        return sqrt(1/v.squared_length())*v;
    }

    Point_3<K> normalize(const Point_3<K> &p) {
        Vector_3<K> v = p-ORIGIN;
        return ORIGIN+normalize(v);
    }

    Vector_3<K> reverse(const Vector_3<K> &v) {
        return -v;
    }

    Point_3<K> reverse(const Point_3<K> &p) {
        Vector_3<K> v = p-ORIGIN;
        return ORIGIN+reverse(v);
    }

    Plane_3<K> bisecting_plane(Point_3<K> p1, Point_3<K> p2) {
        Plane_3<K> result(midpoint(p1, p2), p2-p1);
        return result;
    }

    // Converts a 2D point in spherical coordinates (lat, lon) into a 3D point in cartesian coordinates (x, y, z) on a unit sphere
    Point_3<K> cartesian(const Point_2<K> &p) {
        double lat = to_double(p.x())*TAU/360;
        double colat = TAU/4-lat;
        double lon = to_double(p.y())*TAU/360;
        double x = cos(lon)*sin(colat);
        double y = sin(lon)*sin(colat);
        double z = cos(colat);
        Point_3<K> result(x, y, z);
        return result;
    }

    // Converts a 3D point in cartesian coordinates (x, y, z) into a 2D point in spherical coordinates (lat, lon)
    Point_2<K> spherical(const Point_3<K> &p) {
        double x = to_double(p.x());
        double y = to_double(p.y());
        double z = to_double(p.z());
        double r = sqrt(x*x+y*y+z*z);
        double lon = atan2(y, x)*360/TAU;
        double colat = acos(z/r)*360/TAU;
        double lat = 90-colat;
        Point_2<K> result(lat, lon);
        return result;
    }

    // Converts a 3D point in cartesian coordinates (x, y, z) into a 2D point with a stereographic projection
    // (clat, clon) are the latitude and longitude of the "center" of the stereographic projection, where there is no distortion
    // radius is just the radius of the earth, used to determine scale. The RADIUS constant defined at the top of Utility.h should always be used (average radius of earth in meters)
    Point_2<K> stereographic(const Point_3<K> &p3, double clat, double clon, double radius) {
        double phi1 = clat*TAU/360;
        double lam0 = clon*TAU/360;
        double R = radius;
        Point_2<K> sphere_coord = spherical(p3);
        double phi = sphere_coord.x()*TAU/360;
        double lam = sphere_coord.y()*TAU/360;
        double k = 2*R/(1 + sin(phi1)*sin(phi)+cos(phi1)*cos(phi)*cos(lam-lam0));
        double x = k*cos(phi)*sin(lam-lam0);
        double y = k*(cos(phi1)*sin(phi)-sin(phi1)*cos(phi)*cos(lam-lam0));
        Point_2<K> result(x, y);
        return result;
    }

    // Inverse of stereographic(...), WARNING: I believe it is incorrect
    Point_3<K> inv_stereographic(const Point_2<K> &p2, double clat, double clon, double radius) {
        double phi1 = clat*TAU/360;
        double lam0 = clon*TAU/360;
        double R = radius;
        double x = p2.x();
        double y = p2.y();
        double rho = sqrt(x*x+y*y);
        double c = 2*atan2(rho, 2*R);
        double phi = asin(cos(c)*sin(phi1)+y*sin(c)*cos(phi1)/rho);
        double lam = lam0 + atan2(x*sin(c), rho*cos(phi1)*cos(c) - y*sin(phi1)*sin(c));
        Point_2<K> sphere_coord(phi, lam);
        Point_3<K> result = cartesian(sphere_coord);
        return result;
    }

    // functions for stereographic projections of circles are missing
}

#endif // UTILITY_H
