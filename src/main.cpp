#include <iostream>

#include "kdtree.h"
#include "fileutils.h"

int main(int argc, char **argv) 
{
    const std::string filename =
        argc > 1 ? argv[1] : "../benchmark/benchmark1.csv";

    // Read the dataset as a 1D array, dim consecutive items of dt are a data point
    auto [dataPoints, dim] = parkdtree::utils::readDataset<double>(filename);

    parkdtree::KDTree<double> tree(dataPoints, dim);
    std::cout << tree;
    
}
