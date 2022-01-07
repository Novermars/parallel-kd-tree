#include "utils.h"
#include <exception>

int select_splitting_dimension(int depth, int dims) {
  return depth % dims;
}

int bigger_powersum_of_two(int n) {
  int base = 1;
  int N = 0;
  while (N < n) {
    N += base;
    base *= 2;
  }
  return N;
}

int smaller_powersum_of_two(int n) {
  int base = 1;
  int N = 0;
  while (N < n) {
    N += base;
    base *= 2;
  }
  return N - base / 2;
}

// transform the given DataPoint array in a 1D array such that `dims` contiguous
// items constitute a data point
template <typename T>
std::vector<T> unpack_array(std::vector<DataPoint<T>> const& array) 
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
            unpacked[dp * dim + idx] = array[dp][idx];
        }
    }
    return unpacked;
}

// unpack an array which may contain uninitialized items
// template <typename T>
// std::vector<T> unpack_risky_array(std::vector<DataPoint<T>> const& array, std::vector<bool> const& initialized) 
// {
//     int numDataPts = array.size();
//     int dim = array[0].size();

//     // We must have at least some data
//      if (numDataPts == 0 || dim == 0)
//         throw std::invalid_argument{"Either zero length or zero dimension"}; 
 
//     std::vector<T> unpacked(numDataPts * dim);
//     for (int dp = 0; dp < numDataPts; ++dp) 
//     {
//         for (int idx = 0; idx < dim; ++idx)
//         {
//             unpacked[dp * dim + idx] = initialized[dp] ? array[dp][idx] : std::numeric_limits<int>::min();
//         }
//     }
//   return unpacked;
// }

/*
  This function rearranges branch1 and branch2 into dest such that we first
  take 1 node from branch1 and 1 node from branch2, then 2 nodes from branch1
  and 2 nodes from branch2, then 4 nodes from branch1 and 4 nodes from
  branch2..

  Note that this function is dimensions-safe (i.e. copies all the dimensions).

  Remember to add a split point before this function call (if you need to).
*/
template <typename T>
void rearrange_branches(T *dest, T *branch1, T *branch2,
                        int branches_size, int dims) {
  int already_added = 0;
  // number of nodes in each branch (left and right)at the current level of
  // the tree
  int nodes = 1;
  while (already_added < 2 * branches_size) {
    // we put into the three what's inside the left subtree
    memcpy(dest + already_added * dims, branch1,
                nodes * dims * sizeof(T));
    branch1 += nodes * dims;

    // we put into the three what's inside the right subtree
    memcpy(dest + (nodes + already_added) * dims, branch2,
                nodes * dims * sizeof(T));
    branch2 += nodes * dims;

    // we just added left and right branch
    already_added += nodes * 2;
    // the next level will have twice the number of nodes of the current level
    nodes *= 2;
  }
}

/*
    Convert the given tree to a linked list structure. This assumes that
    the given size is a powersum of two.

    - tree contains the array representation of the tree
    - size is the number of elements in `tree`
    - dims is the number of components for each data point
    - current_level_start contains the index of the first element of tree which
        contains an element of the current node
    - current_level_nodes contains the number of elements in this level of the
        tree (each recursive call multiplies it by two)
    - start_offset contains the offset (starting from current_level_start) for
        the root node of the subtree represented by this recursive call.
*/
// template <typename T>
// KNode<T> *convertToKnodes(std::vector<T> const& tree, int size, int dims,
//                                     int current_level_start,
//                                     int current_level_nodes, int start_offset) {
//   int next_level_start = current_level_start + current_level_nodes * dims;
//   int next_level_nodes = current_level_nodes * 2;
//   int next_start_offset = start_offset * 2;

//   if (next_level_start < size * dims) {
//     auto left = convert_to_knodes(tree, size, dims, next_level_start,
//                                   next_level_nodes, next_start_offset);
//     auto right = convert_to_knodes(tree, size, dims, next_level_start,
//                                    next_level_nodes, next_start_offset + 1);

//     return new KNode<T>(tree + current_level_start +
//                                     start_offset * dims,
//                                 dims, left, right, current_level_start == 0);
//   } else
//     return new KNode<T>(tree + current_level_start +
//                                     start_offset * dims,
//                                 dims, nullptr, nullptr, false);
// }
