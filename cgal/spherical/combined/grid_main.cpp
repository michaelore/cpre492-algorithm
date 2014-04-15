#include <fstream>
#include <iostream>
#include <utility>
#include "Grid_sampling.h"

//#define USE_STEREOGRAPHIC

using namespace std;

int main(int argc, char *argv[]) {
    if (argc < 5) {
        cerr << "Usage: grid_main <k> <ncircs> <nrefines> <filename>" << endl;
        return 1;
    }
    int k = atoi(argv[1]);
    int ncircs = atoi(argv[2]);
    int nrefines = atoi(argv[3]);
    ifstream ifs(argv[4], ifstream::in);
    CGAL::run_grid_sampling(k, ncircs, nrefines, ifs, cout, cerr);
    ifs.close();
    return 0;
}

