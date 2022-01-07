#pragma once

#include "knode.h"
#include "datapoint.h"
#include "fileutils.h"

template <typename T>
class KDTree
{
    KNode<T>* d_root;
    int d_dimension;
    int d_size;

public:
    KDTree(std::vector<T> const& data, int dimension);
    KNode<T>* getRoot();
    template <typename TS>
    friend std::ostream &operator<<(std::ostream &os, const KDTree<TS> &kdTree);
private:
    void generateKDTree(std::vector<T> const& data);
    template <typename Iterator>
    void buildTree(Iterator start, int size, 
                   int depth, int region_width, int region_start_index, 
                   int branch_starting_index, std::vector<bool> const& initialized,
                   std::vector<T> const& splitsTree);
    int selectSplittingDimension(int depth, int dims);
    int sortAndSplit(DataPoint<T> const& array, int size, int axis);
};