#include "kdtree.h"
#include "utils.h"
#include "treeprinter.h"

template <typename T>
KDTree<T>::KDTree(std::vector<T> const& data, int dimension)
:   
    d_dimension{dimension}
{
    d_size = data.size() / dimension;
    d_root = generateKDTree(data);
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
    std::vector splitsTree(splitsTreeSize);

    std::vector<bool> initialized(splitsTreeSize, false);

    buildTree(dataPoints, d_size, 0, 1, 0, 0);

    std::vector<T> flatTree = unpack_risky_array(splitsTree, splitsTreeSize, d_dimension, initialized);
    d_root = convertToKnodes(flatTree, splitsTreeSize, d_dimension, 0, 1, 0);
}

template <typename T>
template <typename Iterator>
void KDTree<T>::buildTree(Iterator start, int size, 
                          int depth, int regionWidth, int regionStartIndex, 
                          int branchStartingIndex, std::vector<bool> const& initialized,
                          std::vector<T> const& splitsTree)
{
    auto idx = regionStartIndex + branchStartingIndex;
    initialized[idx] = true;

    if (size <= 1)
    {
        splitsTree[idx] = start[0];
        return;
    }

    int dimension = selectSplittingDimension(depth, d_dimension);
    int splitPointIdx = sortAndSplit(start, size, d_dimension);

    splitsTree[idx] = start[splitPointIdx];

    regionStartIndex += regionWidth;
    regionWidth *= 2;
    branchStartingIndex *= 2;
    depth += 1;

    build_tree(start + splitPointIdx + 1, size - splitPointIdx - 1, depth,
               regionWidth, regionStartIndex, branchStartingIndex + 1);

    if (splitPointIdx > 0)
      build_tree(start, splitPointIdx, depth, regionWidth,
                 regionStartIndex, branchStartingIndex);

}

template <typename T>
int KDTree<T>::selectSplittingDimension(int depth, int dims) 
{
  return depth % dims;
}

template <typename T>
int KDTree<T>::sortAndSplit(DataPoint<T> const& array, int size, int axis) 
{
  // the second part of median_idx is needed to unbalance the split towards the
  // left region (which is the one which may parallelize with the highest
  // probability).
  int median_idx = size / 2 - ((size + 1) % 2);

  auto comparer = [axis_ = axis](auto lhs, auto rhs){ return lhs[axis_] < rhs[axis_];};

  std::nth_element(std::begin(array), std::begin(array) + median_idx, 
                   std::end(array), comparer);
  // if size is 2 we want to return the first element (the smallest one), since
  // it will be placed into the first empty spot in serial_split
  return median_idx;
}

template <typename TS>
std::ostream &operator<<(std::ostream &os, const KDTree<TS> &kdTree)
{
    return print_tree(os, kdTree.getRoot(), "", false);
    
}