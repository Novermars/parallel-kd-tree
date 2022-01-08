#pragma once
#include "knode.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>
#include <tuple>

#include "fileutils.h"

/*
    This function reads the file row by row, and for each row stores the
    numbers found in an std::vector. The accepted dimension of each data point
    is the dimension of the data point in the first row of the file (a check
    is performed for each row though).
 */
template <typename T>
std::tuple<std::vector<T>, int> readDataset(const std::string &filename) {
    std::ifstream file(filename);

    // local variable which holds the last known number of dimensions per data
    // point. used to check that all data points have the same number of
    // components
    int temp_dims = -1;

    constexpr char separator = ',';

    std::vector<T> lines_buffer;
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
        std::vector<T> row_buffer;

        int start = 0;
        for (unsigned idx = 1; idx < line.length(); ++idx) {
            // we found a component
            if (line[idx] == separator) {
            row_buffer.push_back(std::stod(line.substr(start, idx)));
            start = idx + 1;
            }
        }
        row_buffer.push_back(std::stod(line.substr(start, line.length())));

        // we check that all the components have the same number of dimensions
        if (temp_dims != -1 && temp_dims != static_cast<int>(row_buffer.size()))
            throw std::invalid_argument(
                "Invalid number of dimensions for data point number " +
                std::to_string(lines_buffer.size()));
        temp_dims = row_buffer.size();

        // we put everything into line_buffer
        for (int idx = 0; idx < temp_dims; ++idx)
            lines_buffer.push_back(row_buffer[idx]);
        }
    } else 
    {
        throw std::invalid_argument("File not found.");
    }

    return {lines_buffer, temp_dims};
}


template <typename T>
void write_file(const std::string &filename, KNode<T> *root, int dims)
{
    std::ofstream outdata;
    outdata.open(filename, std::fstream::out);
    if (!outdata) {
        throw std::invalid_argument("File not found.");
    }

    std::queue<KNode<T> *> to_visit;
    to_visit.push(root);

    while (to_visit.size() > 0) {
        KNode<T> *node = to_visit.front();
        to_visit.pop();

        for (int idx = 0; idx < dims; ++idx) {
        outdata << node->get_data(idx);
        if (idx < dims - 1)
            outdata << ",";
        }

        if (node->get_left() != nullptr)
            to_visit.push(node->get_left());
        if (node->get_right() != nullptr)
            to_visit.push(node->get_right());

        outdata << '\n';
    }
}
