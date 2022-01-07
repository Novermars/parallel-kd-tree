#pragma once

#include "knode.h"
#include "utils.h"
#include <iostream>

/**
 * @brief Print the content of a node (i.e. the data point stored inside).
 */
template <typename T>
std::ostream &print_node_values(std::ostream &os,
                                const KNode<T> &node) 
{
    //std::cout << node.get_data(0) << '\n';
    os << "(";

    for (int i = 0; i < node.get_dims(); i++) {
        if (i > 0)
        os << ",";
        if (node.get_data(i) == EMPTY_PLACEHOLDER) {
        os << "n/a";
        break;
        } else
        {
            os << node.get_data(i);
        }
    }
    return os << ")";
}

/**
 * @brief Print a node of the tree and print recursively the left and right
 *          nodes (if available).
 */
template <typename T>
std::ostream &print_tree(std::ostream &os, const KNode<T> &node,
                         const std::string &prefix, bool isLeft) 
{
  os << prefix;
  os << (isLeft ? "├──" : "└──");

  // print the value of the node
  print_node_values(os, node);
  os << std::endl;

  // enter the next tree level - left and right branch
  auto left = node.get_left();
  if (left)
    print_tree(os, *left, prefix + (isLeft ? "│   " : "    "), true);
  auto right = node.get_right();
  if (right)
    print_tree(os, *right, prefix + (isLeft ? "│   " : "    "), false);
  return os;
}

/*
    Start the recursion from the root node.
*/
template <typename T>
std::ostream &operator<<(std::ostream &os, const KNode<T> &node) {
  return print_tree(os, node, "", false);
}

