# COATi
Codon-Aware Multiple Sequence Alignments

## Installation

### Dependencies

* CMake 3.12.0 or higher
* Boost 1.47 or higher
* Eigen 3.3

### Compiling

```
	cd coati/build
	cmake -DCMAKE_BUILD_TYPE=Release ..
	make -j4
```


## Pairwise Alignment: coati-alignpair

```
Usage:	coati-alignpair file.fasta [options]

Allowed options:
  -h [ --help ]                 Display this message
  -f [ --fasta ] arg            fasta file path
  -m [ --model ] arg (=m-coati) substitution model: m-coati (default), ecm, m-ecm
  -w [ --weight ] arg           Write alignment score to file
  -o [ --output ] arg           Alignment output file
  -s [ --score ]                Calculate alignment score given marginal model (m-coati by default)
  -r [ --rate ] arg             Substitution rate matrix (CSV)
 ```
