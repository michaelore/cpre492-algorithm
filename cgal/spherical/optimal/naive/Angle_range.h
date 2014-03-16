
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
                range = CGAL::make_object(full_range);
                empty_range = false;
        }

        void restrict_with_point(Point_3 &p) {
            if (!empty_range) {
                Arc_3 opp_arc = get_opposing_arc(circle, p);
                std::vector<CGAL::Object> restricted;
                if (const Arc_3* casted_range = CGAL::object_cast<Arc_3>(&range)) {
                    CGAL::intersection(*casted_range, opp_arc, back_inserter(restricted));
                    if (restricted.size() > 1) {
                        std::cerr << "Something is very wrong!2" << std::endl;
                    } else if (restricted.size() == 1) {
                        range = restricted[0];
                    } else {
                        empty_range = true;
                    }
                } else if (const std::pair<Arc_Point_3, unsigned>* casted_range_with_mult = CGAL::object_cast<std::pair<Arc_Point_3, unsigned> >(&range)) {
                    if (!K().has_on_3_object()(opp_arc, casted_range_with_mult->first)) {
                        empty_range = true;
                    }
                } else {
                    //
                }
            }
        }

        bool is_empty() {
            return empty_range;
        }

    private:
        Circle_3 circle;
        Object range;
        bool empty_range;

        Arc_3 get_opposing_arc(Circle_3 circle, Point_3 point) {
            Plane_3 orthogonal_plane(circle.center(), point-circle.center());

            std::vector<CGAL::Object> orthogonal_points;
            CGAL::intersection(circle, orthogonal_plane, back_inserter(orthogonal_points));
            Arc_Point_3 orthogonal_one = CGAL::object_cast<std::pair<Arc_Point_3, unsigned> >(orthogonal_points[0]).first;
            Arc_Point_3 orthogonal_two = CGAL::object_cast<std::pair<Arc_Point_3, unsigned> >(orthogonal_points[1]).first;

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
