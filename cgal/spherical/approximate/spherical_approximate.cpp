#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Subdivision_method_3.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/number_utils.h>
#include <CGAL/circulator.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include "Angle_well.h"
#include "Inexact_stereographic_projector.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>

#define TAU 6.2831853071
#define RADIUS 6371009

#define USE_STEREOGRAPHIC

using namespace CGAL;

typedef Simple_cartesian<double>                K;
typedef Plane_3<K>                              Plane;
typedef Polyhedron_3<K>                         Polyhedron;
typedef Vector_3<K>                             Vector;
typedef Polyhedron::Vertex_handle               Vertex_handle;
typedef Polyhedron::Vertex_iterator             Vertex_iterator;
typedef Polyhedron::Halfedge_around_vertex_circulator Halfedge_circulator;
typedef Container_from_circulator<Halfedge_circulator> Halfedge_container;
typedef Halfedge_container::iterator            Halfedge_iterator;
typedef Polyhedron::Halfedge_handle             Halfedge_handle;
typedef Search_traits_3<K>                      TreeTraits;
typedef Orthogonal_incremental_neighbor_search<TreeTraits> NN_search;
typedef NN_search::iterator                     NN_iterator;
typedef NN_search::Tree                         Tree;
typedef Aff_transformation_3<K>                 Transform_3;
typedef Inexact_stereographic_projector<K>      Inexact_Stereo_Projector;

Sphere_3<K> UNIT_SPHERE(ORIGIN, 1);

Point_3<K> cartesian(Point_2<K> point2) {
    K::FT colat = (90-point2.x())*TAU/360;
    K::FT lon = point2.y()*TAU/360;
    K::FT x = cos(lon)*sin(colat);
    K::FT y = sin(lon)*sin(colat);
    K::FT z = cos(colat);
    Point_3<K> result(x, y, z);
    return result;
}

inline double clmp(double val) {
    if (val > 1) {
        return 1;
    } else if (val < -1) {
        return -1;
    } else {
        return val;
    }
}

double dist_convert(double c2) {
    double s = acos(clmp(1-c2/2));
    return RADIUS*s;
}

Vector normalize(const Vector &v) {
    return sqrt(1/v.squared_length())*v;
}

Point_3<K> normalize(const Point_3<K> &p) {
    Vector v = p-ORIGIN;
    return ORIGIN+normalize(v);
}

Vector reverse(const Vector &v) {
    return -v;
}

Point_3<K> reverse(const Point_3<K> &p) {
    Vector v = p-ORIGIN;
    return ORIGIN+reverse(v);
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: spherical_approximate <k> <ncircs> <nrefines> <filename>" << std::endl;
        return 1;
    }
    int k = atoi(argv[1]);
    int ncircs = atoi(argv[2]);
    int nrefines = atoi(argv[3]);
    std::istream_iterator<Point_2<K> > iend;
    Tree tree;
    Point_2<K> proj_point;
    std::ifstream ifs(argv[4], std::ifstream::in);
    for (std::istream_iterator<Point_2<K> > it(ifs); it != iend; it++) {
        tree.insert(cartesian(*it));
        proj_point = *it;
    }
    ifs.close();
    Inexact_Stereo_Projector stereo(proj_point.x(), proj_point.y(), RADIUS);
    Polyhedron P;
    std::ifstream polygon_reader("icosahedron.off", std::ifstream::in);
    polygon_reader >> P;
    polygon_reader.close();

    Subdivision_method_3::Loop_subdivision(P, nrefines);

    std::transform(P.points_begin(), P.points_end(), P.points_begin(), (Point_3<K>(*)(const Point_3<K>&))normalize);

    for (Vertex_iterator v_iter = P.vertices_begin(); v_iter != P.vertices_end(); v_iter++) {
        //std::cerr << (v_iter->point()-ORIGIN).squared_length() << std::endl;
        Halfedge_container h_cont(v_iter->vertex_begin());
        K::FT max_dist = 0;
        Point_3<K> max_center;
        //std::cerr << "Vertex:" << std::endl;
        for (Halfedge_iterator h_iter = h_cont.begin(); h_iter != h_cont.end(); h_iter++) {
            Point_3<K> p1 = h_iter->vertex()->point();
            Point_3<K> p2 = h_iter->next()->vertex()->point();
            Point_3<K> p3 = h_iter->next()->next()->vertex()->point();
            Point_3<K> circ_center = normalize(circumcenter(p1, p2, p3));
            K::FT dist = dist_convert((v_iter->point()-circ_center).squared_length());
            if (max_dist < dist) {
                max_center = circ_center;
            }
            max_dist = std::max(max_dist, dist);
            //std::cerr << "  dist: " << dist << std::endl;
        }
        //std::cerr << "  max_dist: " << max_dist << std::endl;
        NN_search NN(tree, v_iter->point());
        std::vector<std::pair<Point_3<K>, K::FT> > near_neighbors;
        int i = 0;
        NN_iterator it = NN.begin();
        for (; i < k+1 && it != NN.end(); i++, it++) {
            near_neighbors.push_back(std::make_pair(it->first, dist_convert(it->second)));
        }
        if (i < k+1) {
            //
        } else {
            Point_3<K> k1_point = near_neighbors[(k+1)-1].first;
            K::FT k1_dist = near_neighbors[(k+1)-1].second;
            K::FT k1_min_rad = min(abs(k1_dist-max_dist), k1_dist);
            K::FT k1_max_rad = k1_dist+max_dist;
            std::vector<Point_3<K> > border;
            for (int j = 1; j <= k; j++) {
                K::FT j_dist = near_neighbors[j-1].second;
                K::FT j_max_rad = j_dist+max_dist;
                if (j_max_rad >= k1_min_rad) {
                    border.push_back(near_neighbors[j-1].first);
                }
            }
            int extra_on_border = border.size();
            border.push_back(k1_point);
            for (; it != NN.end(); i++, it++) {
                K::FT it_dist = dist_convert(it->second);
                K::FT it_min_rad = min(abs(it_dist-max_dist), it_dist);
                if (it_min_rad <= k1_max_rad) {
                    border.push_back(it->first);
                } else {
                    break;
                }
            }
            //std::cerr << "  border.size(): " << border.size() << std::endl;
            std::vector<double> angles;
            Vector reference = v_iter->point()-ORIGIN;
            Plane tangent(v_iter->point(), reference);
            Vector primary = normalize(v_iter->point()-tangent.projection(border[0]));
            for (int i = 0; i < border.size(); i++) {
                Vector secondary = normalize(v_iter->point()-tangent.projection(border[i]));
                Vector cross = cross_product(primary, secondary);
                //std::cerr << to_double(primary*secondary) << std::endl;
                double angle = acos(clmp(to_double(primary*secondary)));
                if (reference*cross < 0) {
                    angle *= -1;
                }
                if (isnan(angle)) {
                    std::cerr << "Warning: NaN detected" << std::endl;
                } else {
                    angles.push_back(angle);
                }
                //std::cout << angle << std::endl;
            }
            long remov = Angle_well::removability(angles);
            if (remov > extra_on_border) {
                //std::cout << border.size() << "\t" << remov << "\t" << k1_dist << std::endl;
                Point_3<K> infinity = reverse(cartesian(proj_point));
                Point_2<K> true_center = stereo(v_iter->point());
                Point_2<K> circle_center;
                if ((v_iter->point()-infinity).squared_length() < (v_iter->point()-near_neighbors[(k+1)-1].first).squared_length()) {
                    circle_center = stereo(reverse(v_iter->point()));
                } else {
                    circle_center = true_center;
                }
                Point_2<K> k1_proj = stereo(k1_point);
                K::FT rad = sqrt((circle_center-k1_proj).squared_length());
                std::cout << true_center.x() << "\t" << true_center.y() << "\t";
                std::cout << circle_center.x() << "\t" << circle_center.y() << "\t";
                std::cout << rad << "\t";
                std::cout << sqrt((true_center-stereo(max_center)).squared_length());
                std::cout << std::endl;
            }
        }
    }

    return 0;
}

