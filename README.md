# Simple Matrix Operation Library

SMOL, as the name suggests, is a small C library for relatively simple matrix operations on rather small matrix objects. 
Most notes in this document are not relevant for the user, but merely serve as a mental note and guidance for myself.
The use and development of this library is mainly for educational purposes. There are several, more sophisticated libraries for intensive linear algebra computations, like BLAS and LAPACK.

## The Data Structure
SMOL operates on what is essentially a pointer to a contigous block of (one-dimensional) memory, which is interpreted as a matrix in row-major order.
The interpretation requires usually two further attributes - if no context is given -, which are the number of rows and the number of columns. 
An entry i,j in the memory array is therefore indexed by [i*columns + j].

The *SMOL_Matrix* is a struct, which encapsulates the array and the size attribute into a single type. This type can be passed around, copied and assigned like any other first-class type.
Have in mind that the internal matrix array is stored as a pointer. This leads to only the pointer being copied on assignment, rather than the memory block the pointer points to. 
This effect has to be taken into account, if a deep copy is desired.

Depending on the size of the rows and columns of the matrix structure, a *SMOL_Matrix* can either be interpreted as a zero value (null matrix), a scalar, 
a proper matrix or a row- or column vector. The type of a matrix can be inferred by either checking its size properties or calling *SMOL_TypeOf(&matrix)*,
which will return an *SMOL_TYPE* enumerated integer value.

It is easy to convert a two-dimensional array into a matrix struct. Given the array *arr[3][2]*, you can do *(SMOL_Matrix){.fields=&arr, .nRows=3, .nCols=2}*
to construct a matrix struct pointing to the given array. If you want a deep copy of the array, you would do a *SMOL_CopyMat* operation on the constructed matrix struct.

## The Operations
Most of the operations follow a losely structured naming convention to convey information to the user about memory allocation and call-by-reference argument alteration.
Functions with the suffix *Matrix* can be generally thought of constructing a new matrix from the given (non-matrix) parameters. These functions will always allocate memory,
so beware that you do not do something like this: *SMOL_MultiplyMat(&SMOL_RandomMatrix(), &SMOL_RandomMatrix());*. As you will have no way to retrieve the pointers to the memory
allocated by the calls to *SMOL_RandomMatrix()*.
The functions with the *Mat* suffix perform in a smilar way, with the exception, that these functions are imperative operations on a given matrix, returning a new matrix struct.
Functions that return a new element, may return the *SMOL_NULLMATRIX* if an error occured. So it is wise to check for a null pointer in the matrix struct, or by checking if
(SMOL_TypeOf(mat) != SMOL_TYPE_NULL_MATRIX), which is in my opinion more readable.
In general, there is no real difference in the way the Matrix and Mat functions perform. It may just be my obsessive desire to needlessly classify everything.
The functions in imperative form without any suffix - like *SMOL_Fill* - will perform an operation on the matrix passed as call-by-reference directly, therefore altering the contents of
the original matrix.
Furthermore, special functions for the use with vector types are reserved. These follow the convention of having a *Vec* suffix following the function name, e.g. *SMOL_LengthVec(&vec)*;

