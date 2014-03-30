#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Subdivision_method_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <iostream>

using namespace CGAL;

typedef Simple_cartesian<double>               Kernel;
typedef Kernel::Point_3                        Point_3;
typedef Polyhedron_3<Kernel>                   Polyhedron;

int main(int argc, char *argv[]) {
    Polyhedron P;
    std::cin >> P;
    int n;
    if (argc > 1) {
        n = atoi(argv[1]);
    } else {
        n = 1;
    }

    Subdivision_method_3::CatmullClark_subdivision(P, n);

    std::cout << P;
    return 0;
}

