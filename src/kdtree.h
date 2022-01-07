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
    //template <typename TS>
    //friend std::ostream &operator<<(std::ostream &os, const KDTree<TS> &kdTree);
private:
    void generateKDTree(std::vector<T> const& data);
    template <typename Iterator>
    void buildTree(Iterator start, int size, 
                   int depth, int region_width, int region_start_index, 
                   int branch_starting_index, std::vector<bool>& initialized,
                   std::vector<DataPoint<T>>& splitsTree);
    int selectSplittingDimension(int depth, int dims);
    template <typename Iterator>
    int sortAndSplit(Iterator start, int size, int axis);
};

#include "kdtree.h"
#include "utils.h"
#include "treeprinter.h"

template <typename T>
KDTree<T>::KDTree(std::vector<T> const& data, int dimension)
:   
    d_dimension{dimension}
{
    d_size = data.size() / dimension;
    generateKDTree(data);
}

template <typename T>
void KDTree<T>::generateKDTree(std::vector<T> const& data)
{
    std::vector<DataPoint<T>> dataPoints(d_size);
    for (int dp = 0; dp < d_size; ++dp)
    {
        dataPoints[dp] = { data[dp * d_dimension + 0], 
                           data[dp * d_dimension + 1],
                           data[dp * d_dimension + 2] };
    }

    int splitsTreeSize = bigger_powersum_of_two(d_size);
    std::vector<DataPoint<T>> splitsTree(splitsTreeSize);

    std::vector<bool> initialized(splitsTreeSize, false);

    buildTree(std::begin(dataPoints), d_size, 0, 1, 0, 0, initialized, splitsTree);

    std::vector<T> flatTree = unpack_risky_array(splitsTree, initialized);
    d_root = convertToKnodes<T>(std::begin(flatTree), splitsTreeSize, d_dimension, 0, 1, 0);
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

    int dimension = selectSplittingDimension(depth, d_dimension);
    int splitPointIdx = sortAndSplit(start, size, dimension);

    splitsTree[idx] = start[splitPointIdx];

    regionStartIndex += regionWidth;
    regionWidth *= 2;
    branchStartingIndex *= 2;
    depth += 1;

    buildTree(start + splitPointIdx + 1, size - splitPointIdx - 1, depth,
               regionWidth, regionStartIndex, branchStartingIndex + 1,
               initialized, splitsTree);

    if (splitPointIdx > 0)
      buildTree(start, splitPointIdx, depth, regionWidth,
                 regionStartIndex, branchStartingIndex,
                 initialized, splitsTree);

}

template <typename T>
int KDTree<T>::selectSplittingDimension(int depth, int dims) 
{
  return depth % dims;
}

template <typename T>
template <typename Iterator>
int KDTree<T>::sortAndSplit(Iterator start, int size, int axis) 
{
  // the second part of median_idx is needed to unbalance the split towards the
  // left region (which is the one which may parallelize with the highest
  // probability).
  int median_idx = size / 2 - ((size + 1) % 2);

  auto comparer = [&](DataPoint<T> const& lhs, DataPoint<T> const& rhs){ return lhs[axis] < rhs[axis];};

  std::nth_element(start, start + median_idx, 
                   start + size, comparer);
  // if size is 2 we want to return the first element (the smallest one), since
  // it will be placed into the first empty spot in serial_split
  return median_idx;
}

/*template <typename TS>
std::ostream &operator<<(std::ostream &os, const KDTree<TS> &kdTree)
{
    return print_tree(os, kdTree.getRoot(), "", false);
    
}*/