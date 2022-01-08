#pragma once
#include "datapoint.h"
#include "knode.h"

#include <algorithm>
#include <limits>

namespace parkdtree::utils
{

/**
 * @brief Return the smallest number N such that N >= n and N is a sum of powers
 *         of two.
 *
 * Return the smallest sum of powers of two greater than $n$.
 *
 * Example:
 *
 *     bigger_powersum_of_two(5) = 7 = 1 + 2 + 4
 *     bigger_powersum_of_two(3) = 3 = 1 + 2
 */

int biggerPowerSumOfTwo(int num)
{
    int base = 1;
    int sum = 0;

    while (sum < num) 
    {
        sum += base;
        base *= 2;
    }

    return sum;
}

int selectSplittingDimension(int depth, int dims)
{
  return depth % dims;
}

template <typename T, typename Iterator>
int sortAndSplit(Iterator start, int size, int axis) 
{
    // the second part of median_idx is needed to unbalance the split towards the
    // left region (which is the one which may parallelize with the highest
    // probability).
    int median_idx = size / 2 - ((size + 1) % 2);

    auto comparer = [&](parkdtree::DataPoint<T> const& lhs, parkdtree::DataPoint<T> const& rhs){ return lhs[axis] < rhs[axis];};

    std::nth_element(start, start + median_idx, 
                    start + size, comparer);
    // if size is 2 we want to return the first element (the smallest one), since
    // it will be placed into the first empty spot in serial_split
    return median_idx;
}

/**
 * @brief Transform the given array (which may contain uninitialized values)
 *          of data points in a 1D array such that `dims` contiguous items
 *          constitute a data point.
 *
 * Transform the given array of data points, that potentially contains several
 * uninitialized items, in a 1D array such that `dims` contiguous items
 * constitute a data point.
 *
 * `uninitialized` items are spotted using the boolean array `initialized`, and
 * are represented in the output with `dims` consecutive `EMPTY_PLACEHOLDER`.
 *
 * @param array 1D array of data points.
 * @param size  Number of data points in the array (i.e. `length(array)`).
 * @param dims  Number of components for each data point.
 * @param initialized A 1D boolean array (same size of `array`) whose i-th
 *                      element is `true` if and only if the i-th element of
 *                      `array` has been initialized.
 * @return data_type* A 1D array of size `size*dims`.
 */

template <typename T>
std::vector<T> unpackRiskyArray(std::vector<parkdtree::DataPoint<T>> const& array, std::vector<bool> const& initialized) 
{
    int numDataPts = array.size();
    int dim = array[0].size();

    // We must have at least some data
     if (numDataPts == 0 || dim == 0)
        throw std::invalid_argument{"Either zero length or zero dimension"}; 
 
    std::vector<T> unpacked(numDataPts * dim);
    for (int dp = 0; dp < numDataPts; ++dp) 
    {
        for (int idx = 0; idx < dim; ++idx)
        {
            unpacked[dp * dim + idx] = initialized[dp] ? array[dp][idx] : std::numeric_limits<int>::min();
        }
    }
  return unpacked;
}

/**
 * @brief Convert the given tree to a kind-of linked list structure. This
 * assumes that the given size is a powersum of two.
 *
 * @param tree The 1D array representation of the tree.
 * @param size Number of data points in `tree`
 * @param dims Number of components for each data point.
 * @param current_level_start Index of the first element of `tree` which
 *                             contains an element of the current node.
 * @param current_level_nodes Number of elements in this level of the tree (each
 *                             recursive call multiplies it by two).
 * @param start_offset Offset starting from `current_level_start` at which is
 *                      located the root node of the subtree represented by this
 *                      recursive call.
 */
template <typename T, typename Iterator>
parkdtree::KNode<T> *convertToKnodes(Iterator tree, int size, int dims,
                                    int current_level_start,
                                    int current_level_nodes, int start_offset) 
{
    int next_level_start = current_level_start + current_level_nodes * dims;
    int next_level_nodes = current_level_nodes * 2;
    int next_start_offset = start_offset * 2;

    if (next_level_start < size * dims) {
        auto left = convertToKnodes<T>(tree, size, dims, next_level_start,
                                    next_level_nodes, next_start_offset);
        auto right = convertToKnodes<T>(tree, size, dims, next_level_start,
                                    next_level_nodes, next_start_offset + 1);

        return new parkdtree::KNode<T>(&(*tree) + current_level_start +
                                        start_offset * dims,
                                        dims, left, right, current_level_start == 0);
    } else
        return new parkdtree::KNode<T>(&(*tree) + current_level_start +
                                        start_offset * dims,
                                        dims, nullptr, nullptr, false);
}
}
