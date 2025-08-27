# GBLA LC "Gröbner Bases by Linear Algebra and Linear Codes"

## Overview

This project is a C++ application designed to study linear codes, particularly focusing on the **set of leader codewords** and their application to the **Code Equivalence Problem**. It is based on the research article:

[**"Computing an Invariant of a Linear Code"**](https://link.springer.com/chapter/10.1007/978-3-030-43120-4_17)
Lecture Notes in Computer Science, Springer, November 2019.  
By Mijail Borges-Quintana, Miguel Ángel Borges-Trenard, Edgar Martínez-Moro, and Gustavo Torres-Guerrero.

The application implements algorithms for computing leader codewords incrementally and provides tools for analyzing and testing code equivalence under permutations.

## Try it Online with GitPod

Run this project immediately in your browser without any local setup:

[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/gussttaav/GBLA_LC)


## Features

### Code Generation and Configuration
- **Generate a New Linear Code:** Input generator matrices for code creation.
- **Select Word Metric:** Configure the metric (e.g., Hamming metric) for calculations.
- **Select Monomial Order:** Set monomial order (e.g., `degrevlex`) for Gröbner representations.
- **Select Word Output Format:** Choose representation format for codewords.

### Code Analysis
- **Gröbner Representation:** Generate the Gröbner basis of the code.
- **Coset Leaders:** Compute coset leaders and perform their analysis.
- **Codewords and Leader Codewords:** Retrieve and categorize codewords and leader codewords.

### Permutation and Equivalence Testing
- **Permutation:** Apply permutations to codes.
- **Equivalence Testing:** Check if two codes are equivalent under permutation.
- **Automorphism Group:** Compute the automorphism group to analyze code symmetry.

## Installation

### Prerequisites
- A C++ compiler (e.g., `g++`).
- GNU `make` utility.

### Build Instructions
The project is configured with a `Makefile` for both **debug** and **release** builds.

1. Clone the repository or download the source code.
2. Open a terminal in the project source directory.
3. Run the following commands:
   - For a debug build:
     ```bash
     make debug
     ```
   - For a release build:
     ```bash
     make release
     ```

### Cleaning Builds
To clean the build artifacts:
```bash
make clean
```
This will remove both debug and release binaries.

## Usage

After building the application, the executable can be found in:

* `bin/Debug/GBLA` (for debug builds)
* `bin/Release/GBLA` (for release builds)

Run the application:
```bash
./bin/Debug/GBLA
```
or
```bash
./bin/Release/GBLA
```

You can store the definitions of the linear codes you want to work with in a file named codes.txt, which must be located in the same folder as the GBLA executable.

The finite field over which the code is constructed is specified in a line that begins with the following format: 
--> GFOrd[, primitive] 
Where Ord is the order of the field, and if it is not prime, a primitive element must be specified. (See the [codes.txt](codes.txt) file as an example).

Below the line that specifies the finite field, all codes that are specified will be constructed over this field.
Each code begins on a new line with the name that will identify it in the program when loading it, followed in the same line by the length (n) and dimension (k) of the specified code, formatted as n=num1, k=num2. In the immediately following line, the control matrix must appear, as it is a necessary requirement for creating the code. If you also want to specify the generator matrix, it should appear on the line immediately below the control matrix. The matrix should be specified in a single line.

## Research Context

This project is based on the theoretical results presented in the referenced article, focusing on:

1. Efficient Leader Codeword Computation:
An incremental algorithm to compute a compact set of leader codewords, essential for code analysis.

2. Code Equivalence Problem:
Uses leader codewords to derive an invariant for testing code equivalence, leveraging connections to the Graph Isomorphism Problem.

For detailed theory, refer to the article:
[Computing an Invariant of a Linear Code](https://www.springerprofessional.de/computing-an-invariant-of-a-linear-code/17808632).

## License
This application is open-source and distributed under the [MIT License](LICENSE). Contributions are welcome!

## Acknowledgments
The project is a practical implementation of research by Mijail Borges-Quintana et al. Special thanks to the authors for their contributions to coding theory.