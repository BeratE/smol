# Simple Matrix Operation Library

SMOL, as the name suggests, is a small C library for relatively simple matrix operations on small matrix structures. 
Most notes in this document are not relevant for the user, but merely serve as a mental note and guidance for myself.
The use and development of this library is mainly for educational purposes. There are several, more sophisticated libraries for intensive linear algebra computations, like BLAS and LAPACK.

## The Data Structure
SMOL operates on what is essentially a pointer to a continous block of (one-dimensional) memory, which is interpreted as a matrix in row-major order (compatible to how multi-dimensional arrays are stored).
The interpretation requires usually two further attributes; the number of rows and the number of columns. 
An entry i,j in the memory array of m rows and n columns is therefore indexed by [i*n+j].

*SMOL_Matrix* is a struct, which encapsulates the array and the size attribute into a single type. This type can be passed around, copied and assigned like any other first-class type.
Have in mind that the internal matrix array is stored as a pointer. This leads to only the pointer being copied on assignment, rather than the memory block the pointer points to. 
This effect has to be taken into account, if a deep copy is desired.

Depending on the size of the rows and columns of the matrix structure, a *SMOL_Matrix* can either be interpreted as a zero value (null matrix), a scalar, 
a proper matrix or a row- or column vector. The type of a matrix can be inferred by either checking its size properties or calling *SMOL_TypeOf(&matrix)*,
which will return an *SMOL_TYPE* enumerated integer value.

It is easy to convert a two-dimensional array into a matrix struct. Given the array *arr[3][2]*, you can do *(SMOL_Matrix){.fields=&arr, .nRows=3, .nCols=2}*
to construct a matrix struct pointing to the given array. If you want a deep copy of the array, you would do a *SMOL_CopyMatrix* operation on the constructed matrix struct.

## The Operations
Most of the operations follow a losely structured naming convention to convey information to the user about memory allocation and call-by-reference argument alteration.
Functions with the suffix *Matrix* can be generally thought of constructing a new completely matrix. These functions will always allocate memory, which has to be freed at some point.
So beware that you do not do something like this: *SMOL_Multiply(&A, &SMOL_RandomMatrix());*, as you will have no way to retrieve the pointers to the memory allocated by the call to *SMOL_RandomMatrix()*.

Most operations are declared in imperative form, such as Add, Transpose or Multiply. 
This means, that the operations alters the left-hand side argument - the subject of operation - which is always the first argument passed by reference.
The operation A = A * B translates to *SMOL_Multiply(&A, &B);*. The result of the operation is stored in the left-hand operator.

## Todo
* Variadic Multiplication
* Quaternions

## Licence
This project is licenced under the GNU GPLv3.

