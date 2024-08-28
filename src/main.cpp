#include "Simulator.h"

void simulate()
{
    Simulator sim;
    if (sim.Init())
        sim.Start();
    sim.ShutDown();
}


void neighborhood_search_test()
{
    // REVIEW: Neighborhood search functions
    s32 n = 100;
    std::vector<f64> runtimes(n, 0.0);
    for (s32 i = 0; i < n; i++)
    {
        auto start = std::chrono::system_clock::now();

        //test::naive();                 // p=5000 n=100 mean 41.5364ms
        //test::hash_map();              // p=5000 n=100 mean 4.3266ms;  p=25000 n=100 mean 42.6990ms
        //test::parallel_hash_map();     // p=5000 n=100 mean 26.3149ms; p=25000 n=100 mean 35.9189ms
        //test::unordered_map();         // p=5000 n=100 mean 4.6511ms;  p=25000 n=100 mean 48.1093ms
        //test::precompute_neighbor();   // p=5000 n=100 mean 8.3499ms 

        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration<f64, std::milli>(end - start);
        runtimes.push_back(elapsed.count());
    }
    f64 sum = 0.0;
    for (s32 i = 0; i < runtimes.size(); i++)
        sum += runtimes[i];
    f64 mean_elapsed_time = sum / n;
    std::printf("n=%d mean %.4fms\n", n, mean_elapsed_time);
}


//#define TEST

int main()
{
#ifdef TEST
    neighborhood_search_test();
#else
    simulate();
#endif
    return 0;
}