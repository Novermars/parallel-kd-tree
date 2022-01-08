#pragma once

#include <iosfwd>
#include <chrono>
#include <omp.h>
#include <cmath>

#include "knode.h"
#include "datapoint.h"
#include "fileutils.h"
#include "utils.h"


namespace parkdtree
{
template <typename T>
class KDTree
{
    KNode<T>* d_root;
    int d_dimension;
    int d_size;

public:
    KDTree(std::vector<T> const& data, int dimension);
    ~KDTree();
    KNode<T>* getRoot();
    template <typename TS>
    friend std::ostream& operator<<(std::ostream& os, const KDTree<TS>& kdTree);

private:
    void generateKDTree(std::vector<T> const& data);
    template <typename Iterator>
    void buildTree(Iterator start, int size, 
                   int depth, int region_width, int region_start_index, 
                   int branch_starting_index, std::vector<bool>& initialized,
                   std::vector<DataPoint<T>>& splitsTree);
    std::ostream& print_node_values(std::ostream &os,
                                    const KNode<T> &node) const;
    std::ostream &print_tree(std::ostream &os, const KNode<T> &node,
                         const std::string &prefix, bool isLeft) const;
};

template <typename T>
KDTree<T>::KDTree(std::vector<T> const& data, int dimension)
:   
    d_dimension{dimension}
{
    d_size = data.size() / dimension;
    auto start = std::chrono::high_resolution_clock::now();
    generateKDTree(data);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "It took: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << '\n';
}

template <typename T>
KDTree<T>::~KDTree()
{
    // All the nodes will get deleted recursively 
    delete d_root;
}

template <typename T>
void KDTree<T>::generateKDTree(std::vector<T> const& data)
{
    std::vector<DataPoint<T>> dataPoints(d_size, DataPoint<T>(d_dimension));
    for (int dp = 0; dp < d_size; ++dp)
    {
        for (int idx = 0; idx < d_dimension; ++idx)
        {
            dataPoints[dp][idx] = data[dp * d_dimension + idx];
        }
    }

    int splitsTreeSize = parkdtree::utils::biggerPowerSumOfTwo(d_size);
    std::vector<DataPoint<T>> splitsTree(splitsTreeSize);

    std::vector<bool> initialized(splitsTreeSize, false);
    auto start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel
    {
    #pragma omp single
    {
    buildTree(std::begin(dataPoints), d_size, 0, 1, 0, 0, initialized, splitsTree);
    }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Build tree took: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << '\n';
    std::vector<T> flatTree = parkdtree::utils::unpackRiskyArray(splitsTree, initialized);
    d_root = parkdtree::utils::convertToKnodes<T>(std::begin(flatTree), splitsTreeSize, d_dimension, 0, 1, 0);
}

template <typename T>
template <typename Iterator>
void KDTree<T>::buildTree(Iterator start, int size, 
                          int depth, int regionWidth, int regionStartIndex, 
                          int branchStartingIndex, std::vector<bool>& initialized,
                          std::vector<DataPoint<T>>& splitsTree)
{
    auto idx = regionStartIndex + branchStartingIndex;
    initialized[idx] = true;

    if (size <= 1)
    {
        splitsTree[idx] = start[0];
        return;
    }

    int dimension = parkdtree::utils::selectSplittingDimension(depth, d_dimension);
    int splitPointIdx = parkdtree::utils::sortAndSplit<T>(start, size, dimension);

    splitsTree[idx] = start[splitPointIdx];

    regionStartIndex += regionWidth;
    regionWidth *= 2;
    branchStartingIndex *= 2;
    depth += 1;

    int threadNum = omp_get_thread_num();
    int maxDepth = std::log2(omp_get_num_threads());
    int surplusProcesses = omp_get_num_threads() - static_cast<int>(std::pow(2.0, maxDepth));
    bool stopSpawningTasks =
        depth > maxDepth + 1 ||
        (depth == maxDepth + 1 && threadNum >= surplusProcesses);

    #pragma omp task default(shared) final(stopSpawningTasks)
    {
        buildTree(start + splitPointIdx + 1, size - splitPointIdx - 1, depth,
                  regionWidth, regionStartIndex, branchStartingIndex + 1,
                  initialized, splitsTree);
    }

    if (splitPointIdx > 0)
    {        
        buildTree(start, splitPointIdx, depth, regionWidth,
                  regionStartIndex, branchStartingIndex,
                  initialized, splitsTree);        
    }

    #pragma omp taskwait

}

template <typename T>
std::ostream& operator<<(std::ostream& os, const KDTree<T>& kdTree)
{
    return kdTree.print_tree(os, *kdTree.d_root, "", false);    
}

template <typename T>
std::ostream &KDTree<T>::print_node_values(std::ostream &os,
                                           const KNode<T> &node) const
{
    os << "(";

    for (int idx = 0; idx < node.get_dims(); ++idx) 
    {
        if (idx > 0)
            os << ",";

        if (node.get_data(idx) == std::numeric_limits<int>::min())
        {
            os << "n/a";
            break;
        } else
            os << node.get_data(idx);
    }
    return os << ")";
}

template <typename T>
std::ostream &KDTree<T>::print_tree(std::ostream &os, const KNode<T> &node,
                                    const std::string &prefix, bool isLeft) const
{
    os << prefix;
    os << (isLeft ? "├──" : "└──");

    // print the value of the node
    print_node_values(os, node);
    os << '\n';

    // enter the next tree level - left and right branch
    auto left = node.get_left();
    if (left)
        print_tree(os, *left, prefix + (isLeft ? "│   " : "    "), true);
    auto right = node.get_right();
    if (right)
        print_tree(os, *right, prefix + (isLeft ? "│   " : "    "), false);
    return os;
}
}