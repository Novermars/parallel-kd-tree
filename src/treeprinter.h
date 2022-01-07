#pragma once

#include "knode.h"
#include <iostream>

/**
 * @brief Auxiliary function which pretty-prints a k-d tree.
 *
 * This auxiliary function prints the given k-d tree in a pretty way, i.e.
 * the tree develops horizontally on the screen to exploit the space available
 * in the best possible way.
 *
 * @param os   The output file descriptor.
 * @param node The root of the k-d tree to be printed.
 * @return std::ostream& The output file descriptor.
 */
template <typename T>
std::ostream &operator<<(std::ostream &os, const KNode<T> &node);
