# Using LAPACK in C++

## Download and Install LAPACK

* The source code can be downloaded from [here](https://netlib.org/lapack/), e.g. `lapack-3.11.0.tar.gz`. I put it into `~/progs`.

* Unzip the file using

  ```bash
  cd ~/progs
  tar -zxvf lapack-3.11.0.tar.gz
  cd lapack-3.11.0
  ```

* Before building LAPACK, `gfortran` has to be installed using `sudo apt install gfortran`.

* For simplicity, I use `cmake` to install LAPACK in WSL with the following commands. Notice that by default static libraries (`*.a`) are built, so I explicitly specify the option `-DBUILD_SHARED_LIBS=on` to build dynamic libraries.

  ```bash
  mkdir build
  cd build
  cmake -DBUILD_SHARED_LIBS=on -DCMAKE_INSTALL_LIBDIR=~/progs/lib/lapack ..
  cmake --build . -j --target install
  ```

* After the last command, we will see that the following files are installed in `~/progs/lib/lapack`:

  ```bash
  libblas.so
  libblas.so.3
  libblas.so.3.11.0
  liblapack.so
  liblapack.so.3
  liblapack.so.3.11.0
  ```

* Up to now, the compiler does not know where to find them yet. I add the following lines in `~/.bashrc`:

  ```bash
  export LIBRARY_PATH=$LIBRARY_PATH:/home/flcon/progs/lib/lapack
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBRARY_PATH
  ```

  `LIBRARY_PATH` includes paths to search for libraries when compiling using `g++` and `LD_LIBRARY_PATH` includes paths to search for libraries when running executables.

## A concrete example using `DGESV` to solve $Ax=B$

### Problem

Now consider the following two matrices:
$$
A=\left[\begin{array}{ccc}
1 & 3 & 0 \\
2 & 4 & -4 \\
-1 & 9 & 8 \\
\end{array}\right], B=\left[\begin{array}{cc}
1 & 9 \\
2 & 3 \\
9 & 1 \\
\end{array}\right]
$$
and we want to solve for $x$ such that $Ax=B$.

### Code

The code is as follows:

```c++
// Saved as mylapacktest.cpp
// Compile and Run:
//      g++ mylapacktest.cpp -llapack -lblas -o mylapacktest
//      ./mylapacktest
#include <iostream>
#include <iomanip>

using namespace std;

// Declare function
extern "C" void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);

// Print matrix
template <typename T>
void print_matrix(T* mat, int nrows, int ncols) {
    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            cout << fixed << setprecision(3) << mat[i+nrows*j] << "\t";
        }
        cout << endl;
    }
}

int main() {
    int N = 3, NRHS = 2;
    int LDA = N, LDB = N;
    // Matrices
    double A[] = {
        1, 2, -1,
        3, 4, 9,
        0, -4, 8,
    };
    double B[] = {
        1, 2, 9,
        9, 3, 1,
    };
    int IPIV[N];
    int INFO;
    // Print before calculation
    cout << "Pre-calculation" << endl;
    cout << "A=" << endl;
    print_matrix(A, N, N);
    cout << "B=" << endl;
    print_matrix(B, N, NRHS);
    // Calculation
    dgesv_(&N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO);
    // Print after calculation
    cout << endl << "Post-calculation" << endl;
    cout << "INFO=" << INFO << endl;
    cout << "A=" << endl;
    print_matrix(A, N, N);
    cout << "B=" << endl;
    print_matrix(B, N, NRHS);
    cout << "IPIV=" << endl;
    print_matrix(IPIV, N, 1);
    return 0;
}
```

The output is:

```
Pre-calculation
A=
1.000   3.000   0.000
2.000   4.000   -4.000
-1.000  9.000   8.000
B=
1.000   9.000
2.000   3.000
9.000   1.000

Post-calculation
INFO=0
A=
2.000   4.000   -4.000
-0.500  11.000  6.000
0.500   0.091   1.455
B=
-2.750  16.500
1.250   -2.500
-0.625  5.000
IPIV=
2
3
3
```

### Analysis

* First, notice that matrices in C/C++ is stored in a row-major way, but those in Fortran is stored in a column-major way (see more details [here](https://en.wikipedia.org/wiki/Row-_and_column-major_order)). As a result, when initializing `A` and `B`, we have to use the "transpose" of the original matrices.

* Second, the signature of `DGESV` can be found [here](https://netlib.org/lapack/explore-html/d7/d3b/group__double_g_esolve_ga5ee879032a8365897c3ba91e3dc8d512.html#ga5ee879032a8365897c3ba91e3dc8d512). An underscore is sometimes added, depending on the compiler and OS. `DGESV` is a subroutine in Fortran, so outputs are saved in arguments in place. Even when passing `int` or `double`, we have to pass by pointer.

* In the output, `INFO=0` indicates successful calculation. `A` is overwritten by its LU-decomposition $A=PLU$ where
  $$
  L=\left[\begin{array}{ccc}
  1 & 0 & 0 \\
  -0.5 & 1 & 0 \\
  0.5 & 0.091 & 1 \\
  \end{array}\right], U=\left[\begin{array}{ccc}
  2 & 4 &-4 \\
  0 & 11 & 6 \\
  0 & 0 & 1.455 \\
  \end{array}\right], P=\left[\begin{array}{ccc}
  0 & 0 & 1\\
  1 & 0 & 0 \\
  0 & 1 & 0\\
  \end{array}\right]
  $$
  `L` and `U` are stored in `A` after calling the function where the diagonal elements of `L` are 1 by default. `P` is obtained from `IPIV`. Starting from the first row, row $i$ is swapped with row `IPIV[i]`. As a result, in this case, we first swap row 1 with row 2 (`123->213`), then swap row 2 with row 3 (`213->231`). In other words, to obtain `A` from `LU`, we move the 3rd row to the top, leading to the matrix `P` above.

### References

* https://scicomp.stackexchange.com/a/26406/37375
* https://netlib.org/lapack/explore-html/index.html
* https://kaba.hilvi.org/homepage/blog/lapack_blas.htm



