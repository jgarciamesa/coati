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

#include <doctest/doctest.h>
#include <coati/utils.hpp>

#define PRINT_SIZE 100

/* Read fasta format file */
int read_fasta(string file, vector<string>& seq_names, vector<string>& sequences) {
	ifstream input(file);
	if(!input.good()) {
		cerr << "Error opening '" << file << "'." << endl;
		return EXIT_FAILURE;
	}

	string line, name, content;
	while(getline(input, line).good() ) {
		if(line[0] == ';') continue;
		if(line[0] == '>') { // Identifier marker
			if(!name.empty()) {
				sequences.push_back(content);
				name.clear();
			}
			// Add name of sequence
			name = line.substr(1);
			seq_names.push_back(name);
			content.clear();
		} else if(line.empty()) {
			continue;	// omit empty lines
		} else if(!name.empty()) {
			// Remove spaces
			line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
			content += line;
		}
	}
	if(!name.empty()) { // Add last sequence
		sequences.push_back(content);
	}

	return 0;
}

/* Write alignment in Fasta format */
int write_fasta(vector<string> alignment, string output, vector<string> seq_names) {
	ofstream outfile;
	outfile.open(output);
	if(!outfile) {
		cerr << "Opening output file failed.\n";
		return EXIT_FAILURE;
	}

	outfile << ">" << seq_names[0] << endl << alignment[0] << endl;
	outfile << ">" << seq_names[1] << endl << alignment[1] << endl;
	outfile.close();

	return EXIT_SUCCESS;
}

/* Write alignment in PHYLIP format */
int write_phylip(vector<string> alignment, string output, vector<string> seq_names) {
	ofstream outfile;
	outfile.open(output);
	if(!outfile) {
		cerr << "Opening output file failed.\n";
		return EXIT_FAILURE;
	}

	// write aligned sequences to file
	outfile << seq_names.size() << " " << alignment[0].length() << endl;
	int i = PRINT_SIZE-4-max(seq_names[0].length(),seq_names[1].length());
	outfile << seq_names[0] << "\t" << alignment[0].substr(0,i) << endl;
	outfile << seq_names[1] << "\t" << alignment[1].substr(0,i) << endl << endl;
	for(; i<alignment[0].length(); i+=PRINT_SIZE) {
		outfile << alignment[0].substr(i,PRINT_SIZE) << endl;
		outfile << alignment[1].substr(i,PRINT_SIZE) << endl;
		outfile << endl;
	}

	return EXIT_SUCCESS;
}

// TEST_CASE("[utils.cc] write_phylip") {
//
// }

/* Hamming distance between two codons */
int cod_distance(uint8_t cod1, uint8_t cod2) {
	int distance = 0;

	distance += (((cod1 & 48) >> 4) == ((cod2 & 48) >> 4) ? 0:1);
	distance += (((cod1 & 12) >> 2) == ((cod2 & 12) >> 2) ? 0:1);
	distance += ((cod1 & 3) == (cod2 & 3) ? 0:1);

	return distance;
}

/* Cast codon to position in codon list AAA->1, AAAC->2 ... TTT->63 */
int cod_int(string codon) {
	return ((uint8_t) nt4_table[codon[0]]<<4)+((uint8_t) nt4_table[codon[1]]<<2)
		+((uint8_t) nt4_table[codon[2]]);
}


/* Read substitution rate matrix from a CSV file */
int parse_matrix_csv(string file, Matrix64f& P, double& br_len) {
	ifstream input(file);
	if(!input.good()) {
		cerr << "Error opening '" << file << "'." << endl;
		return EXIT_FAILURE;
	}

	string line;
	if(input.good()) {
		// Read branch length
		getline(input,line);
		br_len = stod(line);
	}

	vector<string> vec;
	int count = 0;

	while(getline(input, line)) {
		boost::algorithm::split(vec,line,boost::is_any_of(","));
		P(cod_int(vec[0]), cod_int(vec[1])) = stod(vec[2]);
		count++;
	}

	if(count != 4096){
		cout << "Error reading substitution rate CSV file. Exiting!" << endl;
		return EXIT_FAILURE;
	}

	return 0;
}
