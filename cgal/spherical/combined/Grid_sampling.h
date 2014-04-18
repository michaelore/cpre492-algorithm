#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Subdivision_method_3.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/circulator.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include "KCircle.h"
#include "Utility.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>

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
    std::istream *input;
    std::ostream *output;
    std::ostream *error;

    Tree tree;
    Point_2<K> proj_point;
    std::vector<KCircle> kcircles;

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

    void calc_kcircle_radius(KCircle &kcirc) {
        kcirc.rad = distance_func(kcirc.center);
    }

public:
    Grid_sampling(int pk, int pncircs, std::istream &pinput, std::ostream &poutput, std::ostream &perror) :
            k(pk), ncircs(pncircs) {
        input = &pinput;
        output = &poutput;
        error = &perror;
        std::istream_iterator<Point_2<K> > iend;
        for (std::istream_iterator<Point_2<K> > it(*input); it != iend; it++) {
            tree.insert(cartesian(*it));
            proj_point = *it;
        }
    }

    void sample_with_geodesic_grid(int nrefines) {
        Polyhedron_3<K> geodesic_grid;
        std::ifstream polygon_reader("icosahedron.off", std::ifstream::in);
        polygon_reader >> geodesic_grid;
        polygon_reader.close();

        Subdivision_method_3::Loop_subdivision(geodesic_grid, nrefines);

        std::transform(geodesic_grid.points_begin(), geodesic_grid.points_end(), geodesic_grid.points_begin(), (Point_3<K>(*)(const Point_3<K>&))normalize);
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
            KCircle partial_circ(point, -1, max_dist);
            kcircles.push_back(partial_circ);
        }
    }

    void run_sequential() {
        for (int i = 0; i < kcircles.size(); i++) {
            calc_kcircle_radius(kcircles[i]);
        }
    }

    void output_results() {
        std::sort(kcircles.begin(), kcircles.end(), KCircle::compare_radii);
        int nkcircles = kcircles.size();
        int limit = std::min(nkcircles, ncircs);
        for (int i = 0; i < limit; i++) {
            kcircles[i].project_and_display(*output);
        }
    }
};

namespace CGAL {
    void run_grid_sampling(int k, int ncircs, int nrefines, std::istream &input, std::ostream &output, std::ostream &error) {
        Grid_sampling gs(k, ncircs, input, output, error);
        gs.sample_with_geodesic_grid(nrefines);
        gs.run_sequential();
        gs.output_results();
    }
}
