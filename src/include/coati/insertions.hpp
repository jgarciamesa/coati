/*
# Copyright (c) 2021 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
*/

#ifndef INSERTIONS_HPP
#define INSERTIONS_HPP

#include <Eigen/Sparse>
#include <iostream>
#include <numeric>
#include <vector>

using namespace std;

typedef Eigen::SparseVector<int, Eigen::RowMajor> SparseVectorInt;

struct insertion_data_t {
    vector<string> sequences, names;
    SparseVectorInt insertions;
    insertion_data_t(const vector<string>& s, const vector<string>& n,
                     SparseVectorInt i)
        : sequences{s}, names{n}, insertions{i} {}
    insertion_data_t(const string& s, const string& n, SparseVectorInt i)
        : sequences{{s}}, names{{n}}, insertions{i} {}
    insertion_data_t() : sequences{}, names{}, insertions{SparseVectorInt()} {}
};

bool insertion_flags(const string& ref, const string& seq,
                     SparseVectorInt& insertions_vector);
bool merge_indels(vector<insertion_data_t>& ins_data,
                  insertion_data_t& merged_data);
void add_gap(vector<insertion_data_t>& ins_data, vector<int> seq_indexes,
             int pos);

#endif
