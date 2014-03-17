
#ifndef ANGLE_RANGE_H
#define ANGLE_RANGE_H

#include <vector>
#include <utility>
#include <CGAL/Object.h>

namespace CGAL {
    template <class K>
    class Angle_range {
        typedef typename K::Point_2  Point_2;
        typedef typename K::Point_3  Point_3;
        typedef typename K::Line_3   Line_3;
        typedef typename K::Plane_3  Plane_3;
        typedef typename K::Circle_3 Circle_3;
        typedef typename K::Sphere_3 Sphere_3;
        typedef typename K::Circular_arc_3 Arc_3;
        typedef typename K::Circular_arc_point_3 Arc_Point_3;

    public:
        Angle_range(Circle_3 &c) : circle(c) {
                Arc_3 full_range = Arc_3(c);
                available_sets.push_back(make_object(full_range));
        }

        void restrict_with_point(Point_3 &p) {
            Arc_3 opp_arc = get_opposing_arc(circle, p);
            std::vector<Object> restricted;
            std::pair<Arc_Point_3, unsigned> app;
            Arc_3 arc;
            for (int i = 0; i < available_sets.size(); i++) {
                if (assign(app, available_sets[i])) {
                    if (K().has_on_3_object()(opp_arc, app.first)) {
                        restricted.push_back(make_object(app.first));
                    }
                } else if (assign(arc, available_sets[i])) {
                    intersection(opp_arc, arc, back_inserter(restricted));
                }
            }
            available_sets = restricted;
        }

        bool is_empty() {
            return available_sets.empty();
        }

    private:
        Circle_3 circle;
        std::vector<Object> available_sets;

        Arc_3 get_opposing_arc(Circle_3 circle, Point_3 point) {
            Plane_3 orthogonal_plane(circle.center(), point-circle.center());

            std::vector<Object> orthogonal_points;
            intersection(circle, orthogonal_plane, back_inserter(orthogonal_points));
            Arc_Point_3 orthogonal_one = object_cast<std::pair<Arc_Point_3, unsigned> >(orthogonal_points[0]).first;
            Arc_Point_3 orthogonal_two = object_cast<std::pair<Arc_Point_3, unsigned> >(orthogonal_points[1]).first;

            Arc_3 arc_one(circle, orthogonal_one, orthogonal_two);
            Arc_3 arc_two(circle, orthogonal_two, orthogonal_one);
            if (K().has_on_3_object()(arc_one, point)) {
                return arc_two;
            } else if (K().has_on_3_object()(arc_two, point)) {
                return arc_one;
            } else {
                std::cerr << "Something is very wrong!1" << std::endl;
                return arc_one;
            }
        }
    };
}
#endif // ANGLE_RANGE_H
