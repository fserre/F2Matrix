# F2Matrix
F2Matrix is a linear algebra library for small matrices (of maximal size $8\times 8$) in $\mathbb F2$, i.e. the Galois field with two elements. It supports matrix addition (and therefore subtraction), multiplication, transposition, inversion and rank computation. 

This library aims to provide high performance, an relies heavily on bit manipulation instructions. Particularly, it requires a BMI2 capable processor (e.g. Intel Haswell or newer).

### Use
Include the file `matrix.hpp` in your .cpp file, and turn on the flag allowing intrinsics on your compiler.

##### Creating and displaying a matrix
Matrices are represented using the structure `Matrix`. Template parameters `m`and `n` correspond repectively to the number of rows and columns. The value of a Matrix is stored as a 64 bits integer, in row major order, where the least significant bit is the upper right value of the matrix. Rows are always stored as 8 bits; matrices with less than 8 columns ignore the most significant bits of each row.

The method `print` can be used to display them:

```
Matrix<3,3> p(0x010101LL); //creates a 3x3 matrix m with all bits on the right column set to 1
p.print(); //diplays this matrix
```

The method `printTex` outputs on an output stream the corresponding TeX text (for better visualisation).

##### Addition and multiplication
Addition and multiplication of matrices can be performed using the corresponding C++ operator:

```
Matrix<3,3> q(0x010204LL); //creates a 3x3 matrix m with all bits on the inverse diagonal set to 1
(p + q).print(); //Addition of two matrices
Matrix<3,5> r(0x150A1FLL); //creates a 3x5 matrix
(q * r).print(); //Multiplication of two matrices
```

##### Rank
The rank of a matrix can be computed using the method `rk`:
```
std::cout << q.rk() << std::endl; //Rank of the matrix q
```

For improved performance in the case of small matrices, ranks can be precomputed by calling the function `preCompRanks` once, and then uncommenting the line

```
#include "ranks.hpp"
```

##### Matrix inversion and transposition
These operations can be performed using the (non destructive) methods `transpose` and `inverse`. In the case where the matrix is singular, the method `inverse` returns a zero filled matrix.

##### Submatrix
A submatrix can be extracted using the method `getSubMatrix<i,j,h,w>()`, where `i` and `j` are the coordinates of the upper right corner of the submatrix with respect to the upper right corner of the main matrix, and `h` and `w` the size of the submatrix.

