#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <climits>

using namespace std;
 
#include <omp.h>
#define NUM_RANDS 1000000000

int NUM_THREADS;

// Function that will compute the parallel data
double getSampleData() {
     srand (time(NULL));
     double inside = 0.0;
     for(int i = 0; i <= NUM_RANDS / NUM_THREADS; i++) {
             double x = (rand() % INT_MAX) / (double)(INT_MAX);
             double y = (rand() % INT_MAX) / (double)(INT_MAX);
             double radius = sqrt(pow(x, 2.0)+pow(y, 2.0));
             if(radius <= 1) inside++;
     }
     return inside;
}

int main(int argc, char *argv[])
{
  // Create a time for program start
  time_t beginning;
  time(&beginning);
  
  // Thread id and number of threads respectfully
  int th_id, nthreads;
  
  // Total number of distances located within the circle
  double insideCircle = 0.0;
  
  // Get and set number of threads to run in parallel
  NUM_THREADS=omp_get_max_threads();
  omp_set_num_threads(NUM_THREADS);
  
  // Main function used by every thread
  #pragma omp parallel private(th_id) shared(nthreads)
  {
    // Create a time for program finished
    time_t end;
    
    // Get the thread number and return value from the function call
    th_id = omp_get_thread_num();
    double ret = getSampleData();
    
    // Define what to do in critical region (No other thread can access a variable in this region at the same time
    #pragma omp critical
    {
      insideCircle += ret;
      time(&end);
      cout << "Time to finish thread " << th_id << ": " << difftime(end, beginning) << '\n';
    }
    // Wait for all threads to complete
    #pragma omp barrier
    
    // Master thread runs when all threads have finished
    #pragma omp master
    {
      nthreads = omp_get_num_threads();
      cout << "There are " << nthreads << " threads" << '\n';
    }
  }
  // Compute result
  double ratio = 4.0 * (insideCircle / (double)(NUM_RANDS));
  cout << "Ratio = " << ratio;
  
  return 0;
}
