#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Subdivision_method_3.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/circulator.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include "Angle_well.h"
#include "Utility.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>

//#define USE_STEREOGRAPHIC

using namespace CGAL;

typedef Polyhedron_3<K>::Vertex_handle               Vertex_handle;
typedef Polyhedron_3<K>::Vertex_iterator             Vertex_iterator;
typedef Polyhedron_3<K>::Halfedge_around_vertex_circulator Halfedge_circulator;
typedef Container_from_circulator<Halfedge_circulator> Halfedge_container;
typedef Halfedge_container::iterator            Halfedge_iterator;
typedef Polyhedron_3<K>::Halfedge_handle             Halfedge_handle;
typedef Search_traits_3<K>                      TreeTraits;
typedef Orthogonal_incremental_neighbor_search<TreeTraits> NN_search;
typedef NN_search::iterator                     NN_iterator;
typedef NN_search::Tree                         Tree;

class Grid_sampling {
private:
    int k;
    int ncircs;
    int nrefines;
    std::istream *input;
    std::ostream *output;
    std::ostream *error;

    Tree tree;
    Point_2<K> proj_point;
    Polyhedron_3<K> geodesic_grid;

    double distance_func(Point_3<K> p) {
        NN_search NN(tree, p);
        int i = 0;
        NN_iterator it = NN.begin();
        double k1_dist = 0;
        for (; i < k+1 && it != NN.end(); i++, it++) {
            k1_dist = chord_to_arc(sqrt(it->second));
        }
        if (i < k+1) {
            k1_dist = chord_to_arc(std::sqrt(2*RADIUS));
        }
        return k1_dist;
    }

public:
    Grid_sampling(int pk, int pncircs, int pnrefines, std::istream &pinput, std::ostream &poutput, std::ostream &perror) :
            k(pk), ncircs(pncircs), nrefines(pnrefines) {
        input = &pinput;
        output = &poutput;
        error = &perror;
        std::istream_iterator<Point_2<K> > iend;
        for (std::istream_iterator<Point_2<K> > it(*input); it != iend; it++) {
            tree.insert(cartesian(*it));
            proj_point = *it;
        }
        std::ifstream polygon_reader("icosahedron.off", std::ifstream::in);
        polygon_reader >> geodesic_grid;
        polygon_reader.close();

        Subdivision_method_3::Loop_subdivision(geodesic_grid, nrefines);

        std::transform(geodesic_grid.points_begin(), geodesic_grid.points_end(), geodesic_grid.points_begin(), (Point_3<K>(*)(const Point_3<K>&))normalize);
    }

    void exec_output() {
        for (Vertex_iterator v_iter = geodesic_grid.vertices_begin(); v_iter != geodesic_grid.vertices_end(); v_iter++) {
            Point_3<K> point = v_iter->point();
            Halfedge_container h_cont(v_iter->vertex_begin());
            double max_dist = 0;
            for (Halfedge_iterator h_iter = h_cont.begin(); h_iter != h_cont.end(); h_iter++) {
                Point_3<K> p1 = h_iter->vertex()->point();
                Point_3<K> p2 = h_iter->next()->vertex()->point();
                Point_3<K> p3 = h_iter->next()->next()->vertex()->point();
                Point_3<K> circ_center = normalize(circumcenter(p1, p2, p3));
                double dist = great_circle_dist(point, circ_center);
                max_dist = std::max(max_dist, dist);
            }
            double k1_dist = distance_func(point);
            #ifdef USE_STEREOGRAPHIC
            /*
            Point_3<K> infinity = reverse(cartesian(proj_point));
            Point_2<K> true_center = stereo(v_iter->point());
            Point_2<K> circle_center;
            if ((v_iter->point()-infinity).squared_length() < (v_iter->point()-near_neighbors[(k+1)-1].first).squared_length()) {
                circle_center = stereo(reverse(v_iter->point()));
            } else {
                circle_center = true_center;
            }
            Point_2<K> k1_proj = stereo(k1_point);
            double rad = sqrt((circle_center-k1_proj).squared_length());
            output << true_center.x() << "\t" << true_center.y() << "\t";
            output << circle_center.x() << "\t" << circle_center.y() << "\t";
            output << rad << "\t";
            output << sqrt((true_center-stereo(max_center)).squared_length());
            output << std::endl;
            */
            #else
            Point_2<K> coords = spherical(point);
            *output << coords.x() << "\t" << coords.y() << "\t";
            *output << k1_dist << "\t" << max_dist << "\t";
            *output << std::endl;
            #endif
        }
    }
};

namespace CGAL {
    void run_grid_sampling(int k, int ncircs, int nrefines, std::istream &input, std::ostream &output, std::ostream &error) {
        Grid_sampling gs(k, ncircs, nrefines, input, output, error);
        gs.exec_output();
    }
}
