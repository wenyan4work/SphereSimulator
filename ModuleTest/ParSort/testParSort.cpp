#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

#include "sctl/sctl.hpp"

template <class T>
void testParSort(const int length) {
    std::vector<T> data(length);
    std::srand(0);
    std::generate(data.begin(), data.end(), std::rand); // Using the C function std::rand()

    // sort
    sctl::omp_par::merge_sort(data.begin(), data.end());

    // print
    std::cout <<"result:------------------------"<<std::endl;;
    for (auto &d : data) {
        std::cout << d << " ";
    }
    std::cout <<"-------------------------------"<<std::endl;;
}

int main(int argc, char **argv) {

    int length = 100;
    if (argc > 1) {
        length = atoi(argv[1]);
    }

    testParSort<int>(length);
    testParSort<size_t>(length);
    testParSort<float>(length);
    testParSort<double>(length);

    return 0;
}