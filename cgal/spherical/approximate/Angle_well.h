
#ifndef ANGLE_WELL_H
#define ANGLE_WELL_H

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>
#include <vector>

#define TAU 6.2831853071

inline double wrapAngle( double angle )
{
        return angle - TAU * floor( angle / TAU );
}

namespace CGAL {

    class Angle_well {

    public:

        Angle_well() {}

        void add_intervals(std::vector<double> &starts, std::vector<double> &ends, std::vector<long> &weights) {
            for (int i = 0; i < starts.size(); i++) {
                double start = wrapAngle(starts[i]);
                double end = wrapAngle(ends[i]);
                int weight = weights[i];
                angles[start].second += weight;
                angles[end].first += weight;
                if (end < start) {
                    base += weight;
                }
            }
        }

        long removability() {
            long minval = base;
            long cumsum = base;
            for (std::map<double, std::pair<long, long> >::iterator it=angles.begin(); it!=angles.end(); it++) {
                cumsum -= it->second.first;
                minval = std::min(minval, cumsum);
                //std::cout << cumsum << " ";
                cumsum += it->second.second;
            }
            //std::cout << std::endl;
            return minval;
        }

        static long removability(std::vector<double> &starts, std::vector<double> &ends, std::vector<long> &weights) {
            Angle_well a_well;
            a_well.add_intervals(starts, ends, weights);
            return a_well.removability();
        }

        static long removability(std::vector<double> &angles) {
            std::vector<double> starts;
            std::vector<double> ends;
            std::vector<long> weights;
            for (int i = 0; i < angles.size(); i++) {
                starts.push_back(angles[i] - TAU/4);
                ends.push_back(angles[i] + TAU/4);
                weights.push_back(1);
            }
            Angle_well a_well;
            a_well.add_intervals(starts, ends, weights);
            return a_well.removability();
        }

    private:
        std::map<double, std::pair<long, long> > angles;
        double base;

    };
}
#endif // ANGLE_WELL_H
