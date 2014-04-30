#include <iostream>
#include <fstream>
#include "Exact_max.h"

//#define DEBUG

using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: exact_main <k> <ncircs> <nrefines> <filename>" << std::endl;
        return 1;
    }
    int k = atoi(argv[1]);
    int ncircs = atoi(argv[2]);
    int nrefines = atoi(argv[3]);
    std::ifstream ifs(argv[4], std::ifstream::in);
    CGAL::run_exact_max(k, ncircs, nrefines, ifs, cout, cerr);
    ifs.close();
    return 0;
}
