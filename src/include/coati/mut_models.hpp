/*
# Copyright (c) 2020 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
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

#ifndef MUT_MODELS_H
#define MUT_MODELS_H

#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/CXX11/Tensor>
#include <coati/utils.hpp>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <map>

typedef Eigen::Matrix<double, 64, 1>Vector64f;
typedef Eigen::Matrix<double, 4,  4>Matrix4f;
typedef Eigen::Matrix<double, 5,  1>Vector5d;

using namespace std;

void mg94_q(Matrix64f& Q);
void mg94_p(Matrix64f& P);
void marginalize_p(Eigen::Tensor<double, 3>& p, Matrix64f& P);
void nts_ntv(uint8_t c1, uint8_t c2, int& nts, int& ntv);
double k(uint8_t c1, uint8_t c2, int model=0);
void ecm_p(Matrix64f& P);
vector<string> gotoh_marginal(vector<string> sequences, float& w, Matrix64f& P_m);
vector<string> gotoh(vector<string> sequences, float& w, Matrix64f& P);
double transition(string codon, int position, char nucleotide, Eigen::Tensor<double, 3>& p);
double transition(string codon, string codon2, Matrix64f& p);
vector<string> backtracking(Eigen::MatrixXd Bd, Eigen::MatrixXd Bp, Eigen::MatrixXd Bq, string seqa, string seqb);
float alignment_score(vector<string> alignment, Matrix64f& P);

#endif
