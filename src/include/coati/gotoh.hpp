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

#ifndef GOTOH_HPP
#define GOTOH_HPP

#include <coati/mutation_coati.hpp>

int mg94_marginal(vector<string> sequences, alignment_t& aln, Matrix64f& P);
int gotoh_noframeshifts(vector<string> sequences, alignment_t& aln,
                        Matrix64f& P_m);
double transition(string codon, int position, char nucleotide,
                  const Eigen::Tensor<double, 3>& p);
int backtracking(Eigen::MatrixXd Bd, Eigen::MatrixXd Bp, Eigen::MatrixXd Bq,
                 string seqa, string seqb, alignment_t& aln);
int backtracking_noframeshifts(Eigen::MatrixXd Bd, Eigen::MatrixXd Bp,
                               Eigen::MatrixXd Bq, string seqa, string seqb,
                               alignment_t& aln);

#endif
