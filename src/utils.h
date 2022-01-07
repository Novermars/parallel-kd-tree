#pragma once
#include "datapoint.h"
#include "knode.h"

#include <algorithm>
#include <limits>

/**
 * @def
 * @brief A placeholder used to fill holes in the 1D array which represents the
 *          three (this is needed for indexing reasons).
 */
#define EMPTY_PLACEHOLDER std::numeric_limits<int>::min()

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
int bigger_powersum_of_two(int n);

/**
 * @brief Return the biggest number N such that N <= n and N is a sum of powers
 *         of two.
 *
 * Return the biggest sum of powers of two smaller than $n$.
 *
 * Example:
 *
 *     smaler_powersum_of_two(5) = 3 = 1 + 2
 *     smaller_powersum_of_two(7) = 7 = 1 + 2 + 4
 */
int smaller_powersum_of_two(int n);

/**
 * @brief Transform the given array of data points in a 1D array such that
 *          `dims` contiguous items constitute a data point.
 *
 * @param array 1D array of data points.
 * @param size  Number of data points in the array (i.e. `length(array)`).
 * @param dims  Number of components for each data point.
 * @return data_type* A 1D array of size `size*dims`.
 */
template <typename T>
std::vector<T> unpack_array(std::vector<DataPoint<T>> const& array);

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
std::vector<T> unpack_risky_array(std::vector<DataPoint<T>> const& array, bool *initialized);

/**
 * @brief Rearrange `branch1`, `branch2` into a single array.
 *
 * Rearranges `branch1` and `branch2` into `dest` in such a way that:
 *
 *  1. We first take 1 node from branch1 and 1 node from branch2,
 *  2. Then 2 nodes from `branch1`, and 2 nodes from `branch2`;
 *  3. Then 4 nodes from `branch1`, and 4 nodes from `branch2`;
 *  4. And so on..
 *
 * @param dest    1D array in which we are going to store the content of
 *                  `branch1`, `branch2`.
 * @param branch1 The first branch, from which we are going to take nodes first.
 * @param branch1 The second branch, from which we are going to take nodes after
 *                  the first.
 * @param branches_size Size of `branch1` and `branch2` (number of data points,
 *                        **not** number of data points times the number of
 *                        dimensions).
 * @param dims    Number of dimensions for each data point.
 */
template <typename T>
void rearrange_branches(T *dest, T *branch1, T *branch2,
                        int branches_size, int dims);

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
template <typename T>
KNode<T> *convertToKnodes(std::vector<T> const& tree, int size, int dims,
                                    int current_level_start,
                                    int current_level_nodes, int start_offset);

