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

#include <doctest/doctest.h>

#include <coati/gotoh.hpp>

/* Dynamic Programming implementation of Marginal MG94 model*/
int mg94_marginal(vector<string> sequences, alignment_t& aln, Matrix64f& P_m) {
    // P matrix for marginal Muse and Gaut codon model
    Eigen::Tensor<double, 3> p(64, 3, 4);

    mg94_marginal_p(p, P_m);

    string seq_a = sequences[0];
    string seq_b = sequences[1];
    int m = sequences[0].length();
    int n = sequences[1].length();

    // ensure that length of first sequence (reference) is multiple of 3
    if(m % 3 != 0) {
        cout << "Reference coding sequence length must be a multiple of 3 ("
             << m << "). Exiting!" << endl;
        exit(EXIT_FAILURE);
    }

    // DP matrices for match/mismatch (D), insertion (P), and deletion (Q)
    Eigen::MatrixXf D = Eigen::MatrixXf::Constant(
        m + 1, n + 1, std::numeric_limits<float>::max());
    Eigen::MatrixXf P = Eigen::MatrixXf::Constant(
        m + 1, n + 1, std::numeric_limits<float>::max());
    Eigen::MatrixXf Q = Eigen::MatrixXf::Constant(
        m + 1, n + 1, std::numeric_limits<float>::max());

    // backtracking info matrices for match/mismatch (Bd), insert (Bp), and
    // deletion (Bq)
    Eigen::MatrixXd Bd = Eigen::MatrixXd::Constant(m + 1, n + 1, -1);
    Eigen::MatrixXd Bp = Eigen::MatrixXd::Constant(m + 1, n + 1, -1);
    Eigen::MatrixXd Bq = Eigen::MatrixXd::Constant(m + 1, n + 1, -1);

    double insertion = log(0.001);
    double deletion = log(0.001);
    double insertion_ext = log(1.0 - (1.0 / 6.0));
    double deletion_ext = log(1.0 - (1.0 / 6.0));
    double no_insertion = log(1.0 - 0.001);
    double no_deletion = log(1.0 - 0.001);
    double no_insertion_ext = log(1.0 / 6.0);
    double no_deletion_ext = log(1.0 / 6.0);

    Vector5d nuc_freqs;
    nuc_freqs << log(0.308), log(0.185), log(0.199), log(0.308), log(0.25);

    // DP and backtracking matrices initialization

    // fill first values on D that are independent
    D(0, 0) = 0.0;  // 0.0;
    Bd(0, 0) = 0;
    D(0, 1) = -insertion - nuc_freqs[nt4_table[seq_b[0]]] - no_insertion_ext;
    P(0, 1) = -insertion - nuc_freqs[nt4_table[seq_b[0]]] - no_insertion_ext;
    Bd(0, 1) = 1;
    Bp(0, 1) = 2;
    D(1, 0) = -no_insertion - deletion - no_deletion_ext;
    Q(1, 0) = -no_insertion - deletion - no_deletion_ext;
    Bd(1, 0) = 2;
    Bq(1, 0) = 2;

    // fill first row of D
    if(n + 1 >= 2) {
        for(int j = 2; j < n + 1; j++) {
            D(0, j) = D(0, j - 1) - insertion_ext -
                      nuc_freqs[nt4_table[seq_b[j - 1]]];
            P(0, j) = P(0, j - 1) - insertion_ext -
                      nuc_freqs[nt4_table[seq_b[j - 1]]];
            Bd(0, j) = 1;
            Bp(0, j) = 1;
        }
    }

    // fill first column of D
    if(m + 1 >= 2) {
        for(int i = 2; i < m + 1; i++) {
            D(i, 0) = D(i - 1, 0) - deletion_ext;
            Q(i, 0) = Q(i - 1, 0) - deletion_ext;
            Bd(i, 0) = 2;
            Bq(i, 0) = 1;
        }
    }

    string codon;
    double p1, p2, q1, q2, d;

    for(int i = 1; i < m + 1; i++) {
        codon = seq_a.substr((((i - 1) / 3) * 3), 3);  // current codon
        for(int j = 1; j < n + 1; j++) {
            // insertion
            p1 = P(i, j - 1) - insertion_ext -
                 nuc_freqs[nt4_table[seq_b[j - 1]]];
            p2 = Bd(i, j - 1) == 0
                     ? D(i, j - 1) - insertion -
                           nuc_freqs[nt4_table[seq_b[j - 1]]] - no_insertion_ext
                 : Bd(i, j - 1) == 1 ? D(i, j - 1) - insertion_ext -
                                           nuc_freqs[nt4_table[seq_b[j - 1]]]
                                     : numeric_limits<double>::max();
            P(i, j) = min(p1, p2);
            Bp(i, j) =
                p1 < p2
                    ? 1
                    : 2;  // 1 is insertion extension, 2 is insertion opening

            // deletion
            q1 = Q(i - 1, j) - deletion_ext;
            q2 = Bd(i - 1, j) == 0
                     ? D(i - 1, j) - no_insertion - deletion - no_deletion_ext
                 : Bd(i - 1, j) == 1 ? D(i - 1, j) - no_deletion_ext - deletion
                                     : D(i - 1, j) - deletion_ext;
            Q(i, j) = min(q1, q2);
            Bq(i, j) =
                q1 < q2 ? 1
                        : 2;  // 1 is deletion extension, 2 is deletion opening

            // match/mismatch
            if(Bd(i - 1, j - 1) == 0) {
                d = D(i - 1, j - 1) - no_insertion - no_deletion -
                    log(transition(codon, (i) % 3, seq_b[j - 1], p));
            } else if(Bd(i - 1, j - 1) == 1) {
                d = D(i - 1, j - 1) - no_deletion -
                    log(transition(codon, (i) % 3, seq_b[j - 1], p));
            } else {
                d = D(i - 1, j - 1) -
                    log(transition(codon, (i) % 3, seq_b[j - 1], p));
            }

            // D[i,j] = highest weight between insertion, deletion, and
            // match/mismatch
            //	in this case, lowest (-log(weight)) value
            if(d < P(i, j)) {
                if(d < Q(i, j)) {
                    D(i, j) = d;
                    Bd(i, j) = 0;
                } else {
                    D(i, j) = Q(i, j);
                    Bd(i, j) = 2;
                }
            } else {
                if(P(i, j) < Q(i, j)) {
                    D(i, j) = P(i, j);
                    Bd(i, j) = 1;
                } else {
                    D(i, j) = Q(i, j);
                    Bd(i, j) = 2;
                }
            }
        }
    }

    aln.weight = D(m, n);  // weight

    // backtracking to obtain alignment
    return backtracking(Bd, Bp, Bq, seq_a, seq_b, aln);
}

/* Dynamic Programming with no frameshifts*/
int gotoh_noframeshifts(vector<string> sequences, alignment_t& aln,
                        Matrix64f& P_m) {
    // P matrix for marginal Muse and Gaut codon model
    Eigen::Tensor<double, 3> p(64, 3, 4);

    mg94_marginal_p(p, P_m);

    string seq_a = sequences[0];
    string seq_b = sequences[1];
    int m = sequences[0].length();
    int n = sequences[1].length();

    // ensure that length of first sequence (reference) is multiple of 3
    if((m % 3 != 0) || (n % 3 != 0)) {
        cout << "The length of both sequences must be a multiple of 3. Exiting!"
             << endl;
        exit(EXIT_FAILURE);
    }

    // DP matrices for match/mismatch (D), insertion (P), and deletion (Q)
    Eigen::MatrixXf D = Eigen::MatrixXf::Constant(
        m + 1, n + 1, std::numeric_limits<float>::max());
    Eigen::MatrixXf P = Eigen::MatrixXf::Constant(
        m + 1, n + 1, std::numeric_limits<float>::max());
    Eigen::MatrixXf Q = Eigen::MatrixXf::Constant(
        m + 1, n + 1, std::numeric_limits<float>::max());

    // backtracking info matrices for match/mismatch (Bd), insert (Bp), and
    // deletion (Bq)
    Eigen::MatrixXd Bd = Eigen::MatrixXd::Constant(m + 1, n + 1, -1);
    Eigen::MatrixXd Bp = Eigen::MatrixXd::Constant(m + 1, n + 1, -1);
    Eigen::MatrixXd Bq = Eigen::MatrixXd::Constant(m + 1, n + 1, -1);

    double insertion = 0.001;
    double deletion = 0.001;
    double insertion_ext = 1.0 - (1.0 / 6.0);
    double deletion_ext = 1.0 - (1.0 / 6.0);

    Vector5d nuc_freqs;
    nuc_freqs << 0.308, 0.185, 0.199, 0.308, 0.25;

    // DP and backtracking matrices initialization

    // fill first values on D that are independent
    D(0, 0) = 0.0;
    Bd(0, 0) = 0;
    D(0, 3) = P(0, 3) = -log(insertion) - log(nuc_freqs[nt4_table[seq_b[0]]]) -
                        log(nuc_freqs[nt4_table[seq_b[1]]]) -
                        log(nuc_freqs[nt4_table[seq_b[2]]]) -
                        log(1.0 - insertion_ext) - 2 * log(insertion_ext);
    D(3, 0) = Q(3, 0) = -log(1.0 - insertion) - log(deletion) -
                        2 * log(deletion_ext) - log(1.0 - deletion_ext);
    Bd(0, 3) = 1;
    Bd(3, 0) = Bp(0, 3) = Bq(3, 0) = 2;

    // fill first row of D
    if(n + 1 >= 6) {
        for(int j = 6; j < n + 1; j += 3) {
            D(0, j) = P(0, j) = D(0, j - 3) - 3 * log(insertion_ext) -
                                log(nuc_freqs[nt4_table[seq_b[j - 3]]]) -
                                log(nuc_freqs[nt4_table[seq_b[j - 2]]]) -
                                log(nuc_freqs[nt4_table[seq_b[j - 1]]]);
            Bd(0, j) = 1;
            Bp(0, j) = 1;
        }
    }

    // fill first column of D
    if(m + 1 >= 6) {
        for(int i = 6; i < m + 1; i += 3) {
            D(i, 0) = D(i - 3, 0) - 3 * log(deletion_ext);
            Q(i, 0) = Q(i - 3, 0) - 3 * log(deletion_ext);
            Bd(i, 0) = 2;
            Bq(i, 0) = 1;
        }
    }

    string codon;
    double p1, p2, q1, q2, d, temp;

    // Cells with only match/mismatch (1,1) & (2,2)
    codon = seq_a.substr(0, 3);
    for(int i = 1; i < 3; i++) {
        D(i, i) = D(i - 1, i - 1) - log(1.0 - insertion) - log(1.0 - deletion) -
                  log(transition(codon, i % 3, seq_b[i - 1], p));
        Bd(i, i) = 0;
    }

    // Second and third row/column (match/mismatch && insertion || deletion)
    for(int i = 1; i < 3; i++) {
        for(int j = i + 3; j < n + 1; j += 3) {  // rows
            codon = seq_a.substr(0, 3);
            // match/mismatch
            if(Bd(i - 1, j - 1) == 0) {
                d = D(i - 1, j - 1) - log(1.0 - insertion) -
                    log(1.0 - deletion) -
                    log(transition(codon, (i) % 3, seq_b[j - 1], p));
            } else if(Bd(i - 1, j - 1) == 1) {
                d = D(i - 1, j - 1) - log(1.0 - deletion) -
                    log(transition(codon, (i) % 3, seq_b[j - 1], p));
            } else {
                d = D(i - 1, j - 1) -
                    log(transition(codon, (i) % 3, seq_b[j - 1], p));
            }
            // insertion
            p1 = P(i, j - 3) - 3 * log(insertion_ext) -
                 log(nuc_freqs[nt4_table[seq_b[j - 3]]]) -
                 log(nuc_freqs[nt4_table[seq_b[j - 2]]]) -
                 log(nuc_freqs[nt4_table[seq_b[j - 1]]]);
            p2 = Bd(i, j - 3) == 0
                     ? D(i, j - 3) - log(insertion) - 2 * log(insertion_ext) -
                           log(nuc_freqs[nt4_table[seq_b[j - 3]]]) -
                           log(nuc_freqs[nt4_table[seq_b[j - 2]]]) -
                           log(nuc_freqs[nt4_table[seq_b[j - 1]]]) -
                           log(1.0 - insertion_ext)
                 : Bd(i, j - 1) == 1
                     ? D(i, j - 1) - 3 * log(insertion_ext) -
                           log(nuc_freqs[nt4_table[seq_b[j - 3]]]) -
                           log(nuc_freqs[nt4_table[seq_b[j - 2]]]) -
                           log(nuc_freqs[nt4_table[seq_b[j - 1]]])
                     : numeric_limits<double>::max();
            P(i, j) = min(p1, p2);
            Bp(i, j) =
                p1 < p2
                    ? 1
                    : 2;  // 1 is insertion extension, 2 is insertion opening
            D(i, j) = P(i, j) < d ? P(i, j) : d;
            Bd(i, j) = P(i, j) < d ? 1 : 0;
        }

        for(int j = i + 3; j < m + 1; j += 3) {            // columns
            codon = seq_a.substr((((j - 1) / 3) * 3), 3);  // current codon
            // match/mismatch
            if(Bd(j - 1, i - 1) == 0) {
                d = D(j - 1, i - 1) - log(1.0 - insertion) -
                    log(1.0 - deletion) -
                    log(transition(codon, (j) % 3, seq_b[i - 1], p));
            } else if(Bd(j - 1, i - 1) == 1) {
                d = D(j - 1, i - 1) - log(1.0 - deletion) -
                    log(transition(codon, (j) % 3, seq_b[i - 1], p));
            } else {
                d = D(j - 1, i - 1) -
                    log(transition(codon, (j) % 3, seq_b[i - 1], p));
            }
            // deletion
            q1 = Q(j - 3, i) - 3 * log(deletion_ext);
            q2 = Bd(j - 3, i) == 0
                     ? D(j - 3, i) - log(1.0 - insertion) - log(deletion) -
                           log(1.0 - deletion_ext) - 2 * log(deletion_ext)
                 : Bd(j - 3, i) == 1 ? D(j - 3, i) - log(1.0 - deletion_ext) -
                                           log(deletion) - 2 * log(deletion_ext)
                                     : D(j - 3, i) - 3 * log(deletion_ext);
            Q(j, i) = min(q1, q2);
            Bq(j, i) =
                q1 < q2 ? 1
                        : 2;  // 1 is deletion extension, 2 is deletion opening
            D(j, i) = Q(j, i) < d ? Q(j, i) : d;
            Bd(j, i) = Q(j, i) < d ? 1 : 2;
        }
    }

    // Cells considering all 3 events (insertion, deletion, match/mismatch)
    for(int i = 3; i < m + 1; i++) {
        codon = seq_a.substr((((i - 1) / 3) * 3), 3);  // current codon
        for(int j = 3; j < n + 1; j += 3) {
            temp = j;
            j += i % 3;
            if(j > n) continue;
            // match/mismatch
            if(Bd(i - 1, j - 1) == 0) {
                d = D(i - 1, j - 1) - log(1.0 - insertion) -
                    log(1.0 - deletion) -
                    log(transition(codon, (i) % 3, seq_b[j - 1], p));
            } else if(Bd(i - 1, j - 1) == 1) {
                d = D(i - 1, j - 1) - log(1.0 - deletion) -
                    log(transition(codon, (i) % 3, seq_b[j - 1], p));
            } else {
                d = D(i - 1, j - 1) -
                    log(transition(codon, (i) % 3, seq_b[j - 1], p));
            }
            // insertion
            p1 = P(i, j - 3) - 3 * log(insertion_ext) -
                 log(nuc_freqs[nt4_table[seq_b[j - 3]]]) -
                 log(nuc_freqs[nt4_table[seq_b[j - 2]]]) -
                 log(nuc_freqs[nt4_table[seq_b[j - 1]]]);
            p2 = Bd(i, j - 3) == 0
                     ? D(i, j - 3) - log(insertion) - 2 * log(insertion_ext) -
                           log(nuc_freqs[nt4_table[seq_b[j - 3]]]) -
                           log(nuc_freqs[nt4_table[seq_b[j - 2]]]) -
                           log(nuc_freqs[nt4_table[seq_b[j - 1]]]) -
                           log(1.0 - insertion_ext)
                 : Bd(i, j - 1) == 1
                     ? D(i, j - 1) - 3 * log(insertion_ext) -
                           log(nuc_freqs[nt4_table[seq_b[j - 3]]]) -
                           log(nuc_freqs[nt4_table[seq_b[j - 2]]]) -
                           log(nuc_freqs[nt4_table[seq_b[j - 1]]])
                     : numeric_limits<double>::max();
            P(i, j) = min(p1, p2);
            Bp(i, j) =
                p1 < p2
                    ? 1
                    : 2;  // 1 is insertion extension, 2 is insertion opening
            // deletion
            q1 = Q(i - 3, j) - 3 * log(deletion_ext);
            q2 = Bd(i - 3, j) == 0
                     ? D(i - 3, j) - log(1.0 - insertion) - log(deletion) -
                           log(1.0 - deletion_ext) - 2 * log(deletion_ext)
                 : Bd(i - 3, j) == 1 ? D(i - 3, j) - log(1.0 - deletion_ext) -
                                           log(deletion) - 2 * log(deletion_ext)
                                     : D(i - 3, j) - 3 * log(deletion_ext);
            Q(i, j) = min(q1, q2);
            Bq(i, j) =
                q1 < q2 ? 1
                        : 2;  // 1 is deletion extension, 2 is deletion opening
            // D[i,j] = highest weight between insertion, deletion, and
            // match/mismatch
            //	in this case, lowest (-log(weight)) value
            if(d < P(i, j)) {
                if(d < Q(i, j)) {
                    D(i, j) = d;
                    Bd(i, j) = 0;
                } else {
                    D(i, j) = Q(i, j);
                    Bd(i, j) = 2;
                }
            } else {
                if(P(i, j) < Q(i, j)) {
                    D(i, j) = P(i, j);
                    Bd(i, j) = 1;
                } else {
                    D(i, j) = Q(i, j);
                    Bd(i, j) = 2;
                }
            }
            j = temp;
        }
    }

    aln.weight = D(m, n);  // weight

    // backtracking to obtain alignment
    return backtracking_noframeshifts(Bd, Bp, Bq, seq_a, seq_b, aln);
}
/* Return value from marginal MG94 model p matrix for a given transition */
double transition(string codon, int position, char nuc,
                  const Eigen::Tensor<double, 3>& p) {
    position = position == 0 ? 2 : --position;

    if(nuc != 'N') {
        return p(cod_int(codon), position, nt4_table[nuc]);
    } else {
        double val = 0.0;
        for(int i = 0; i < 4; i++) {
            val += p(cod_int(codon), position, i);
        }
        return val / 4.0;
    }
}

/* Recover alignment given backtracking matrices for DP alignment */
int backtracking(Eigen::MatrixXd Bd, Eigen::MatrixXd Bp, Eigen::MatrixXd Bq,
                 string seqa, string seqb, alignment_t& aln) {
    int i = seqa.length();
    int j = seqb.length();

    // vector<string> alignment;
    aln.f.seq_data.push_back(string());
    aln.f.seq_data.push_back(string());

    while((i != 0) || (j != 0)) {
        // match/mismatch
        if(Bd(i, j) == 0) {
            aln.f.seq_data[0].insert(0, 1, seqa[i - 1]);
            aln.f.seq_data[1].insert(0, 1, seqb[j - 1]);
            i--;
            j--;
            // insertion
        } else if(Bd(i, j) == 1) {
            while(Bp(i, j) == 1) {
                aln.f.seq_data[0].insert(0, 1, '-');
                aln.f.seq_data[1].insert(0, 1, seqb[j - 1]);
                j--;
            }
            aln.f.seq_data[0].insert(0, 1, '-');
            aln.f.seq_data[1].insert(0, 1, seqb[j - 1]);
            j--;
            // deletion
        } else {
            while(Bq(i, j) == 1) {
                aln.f.seq_data[0].insert(0, 1, seqa[i - 1]);
                aln.f.seq_data[1].insert(0, 1, '-');
                i--;
            }
            aln.f.seq_data[0].insert(0, 1, seqa[i - 1]);
            aln.f.seq_data[1].insert(0, 1, '-');
            i--;
        }
    }

    return 0;
}

/* Recover alignment given backtracking matrices for DP alignment */
int backtracking_noframeshifts(Eigen::MatrixXd Bd, Eigen::MatrixXd Bp,
                               Eigen::MatrixXd Bq, string seqa, string seqb,
                               alignment_t& aln) {
    int i = seqa.length();
    int j = seqb.length();

    // vector<string> alignment;
    aln.f.seq_data.push_back(string());
    aln.f.seq_data.push_back(string());

    while((i != 0) || (j != 0)) {
        // match/mismatch
        if(Bd(i, j) == 0) {
            aln.f.seq_data[0].insert(0, 1, seqa[i - 1]);
            aln.f.seq_data[1].insert(0, 1, seqb[j - 1]);
            i--;
            j--;
            // insertion
        } else if(Bd(i, j) == 1) {
            while(Bp(i, j) == 1) {
                for(int h = 0; h < 3; h++) {
                    aln.f.seq_data[0].insert(0, 1, '-');
                    aln.f.seq_data[1].insert(0, 1, seqb[j - 1]);
                    j--;
                }
            }
            for(int h = 0; h < 3; h++) {
                aln.f.seq_data[0].insert(0, 1, '-');
                aln.f.seq_data[1].insert(0, 1, seqb[j - 1]);
                j--;
            }
            // deletion
        } else {
            while(Bq(i, j) == 1) {
                for(int h = 0; h < 3; h++) {
                    aln.f.seq_data[0].insert(0, 1, seqa[i - 1]);
                    aln.f.seq_data[1].insert(0, 1, '-');
                    i--;
                }
            }
            for(int h = 0; h < 3; h++) {
                aln.f.seq_data[0].insert(0, 1, seqa[i - 1]);
                aln.f.seq_data[1].insert(0, 1, '-');
                i--;
            }
        }
    }

    return 0;
}
