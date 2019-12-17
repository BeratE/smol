#include "matrix.h"
#include <string.h>
#include <time.h>
#include <math.h>

static long _randseed = 0; // RNG Seed, 0 if unitialized. 

void SMOL_PrintError(enum SMOL_STATUS status)
/* Prints the status message in human readable form. */
{
    fprint("#!SMOL: ");
    switch(status) {
    case(SMOL_STATUS_OK):
	fprint("The operation performed successfully. \n");
	break;
    case(SMOL_STATUS_INVALID_TYPE):
	fprint("The type of the parameter is not accepted. \n");
	break;
    case(SMOL_STATUS_INCOMPATIBLE_SIZES):
	fprint("The sizes of the arguments is not compatible. \n");
	break;
    case(SMOL_STATUS_ARRAY_OUT_OF_BOUNDS):
	fprint("The arguments are out of the bounds of the matrix. \n");
	break;
    default:
	fprint("Something something..\n");
    };
}

void SMOL_Free(SMOL_Matrix* mat)
/* Free memory allocated for the given matrix. */
{
    if (mat->fields != NULL) {
	free(mat->fields);
	mat->fields = NULL;
	mat->nRows = 0;
	mat->nCols = 0;
    }
}

void SMOL_FreeV(int count, ...)
/* Free memory allocated by variable number of matrices. */
{
    va_list args;
    va_start(args, count);
    while (count--)
	SMOL_Free(va_arg(args, SMOL_Matrix*));
    va_end(args);
}

int SMOL_TypeOf(const SMOL_Matrix *mat)
/* Return the SMOL_TYPE of the given matrix. */
{
    if (mat->fields == NULL)
	return SMOL_TYPE_NULL;
    if (mat->nCols > 1 && mat->nCols > 1)
	return SMOL_TYPE_MATRIX;
    if (mat->nCols == 1 && mat->nRows > 1)
	return SMOL_TYPE_COLUMN_VECTOR;
    if (mat->nRows == 1 && mat->nCols > 1)
	return SMOL_TYPE_ROW_VECTOR;
    if (mat->nCols == 1 && mat->nRows == 1)
	return SMOL_TYPE_SCALAR;

    return 0;
}

int SMOL_AllocMatrix(SMOL_Matrix *lhs, size_t nRows, size_t nCols)
/* Returns a newly allocated Matrix of given size.
 * ! Allocated memory has to be freed at some point.*/
{
    lhs->nRows = nRows;
    lhs->nCols = nCols;
    lhs->fields = malloc(sizeof(double)*nRows*nCols);
    memset(lhs->fields, 0.0, sizeof(double)*nRows*nCols);

    return SMOL_STATUS_OK;
}

int SMOL_RandomMatrix(SMOL_Matrix *lhs, size_t nRows, size_t nCols)
/* Creates a matrix of given size and fills it with random values. */
{
    if (_randseed == 0) { // Init seed
	_randseed = time(NULL)/2;
	srand(_randseed);
    }
    
    SMOL_AllocMatrix(lhs, nRows, nCols);
    for (size_t r = 0; r < nRows; r++) {
	for (size_t c = 0; c < nCols; c++) {
	    lhs->fields[r*nCols+c] = rand();
	}
    }
    
    return SMOL_STATUS_OK;
}

int SMOL_CopyMatrix(SMOL_Matrix *lhs, const SMOL_Matrix *mat)
/* Returns a deep copy of given matrix. */
{
    SMOL_AllocMatrix(lhs, mat->nRows, mat->nCols);
    memcpy(lhs->fields, mat->fields, sizeof(double)*lhs->nRows*lhs->nCols);
    
    return SMOL_STATUS_OK;
}

int SMOL_EyeMatrix(SMOL_Matrix *lhs, size_t size)
/* Assigns given matrix of to the identity. */
{
    SMOL_AllocMatrix(lhs, size, size);
    for (unsigned int i = 0; i < size; i++)
	lhs->fields[i*size+i] = 1.0;
    
    return SMOL_STATUS_OK;
}

int SMOL_PerspectiveMatrix(SMOL_Matrix *lhs, double fov, double ratio, double near, double far)
/* Construct the perspective matrix from the given parameters. */
{
    SMOL_AllocMatrix(lhs, 4, 4);
    const double t = near * tan((fov/2)*M_PI/180);
    const double r = t * ratio;
    
    SMOL_SetField(lhs, 0, 0, near/r);
    SMOL_SetField(lhs, 1, 1, near/t);
    SMOL_SetField(lhs, 2, 2, (near-far)/(far-near));
    SMOL_SetField(lhs, 2, 3, (-2*near*far)/(far-near));
    SMOL_SetField(lhs, 3, 2, -1.0);
    
    return SMOL_STATUS_OK;
}

int SMOL_CameraMatrix(SMOL_Matrix *lhs, const double* vec3eye, const double* vec3front, const double* vec3up)
/* Construct a camera matrix from the eye position and the given front & up vectors. */
{
    SMOL_Matrix front, right, up;
    SMOL_CopyMatrix(&front, &(SMOL_Matrix){.fields=(double*)vec3front, .nRows=3, .nCols=1});
    SMOL_CopyMatrix(&right, &front);
    SMOL_VectorCross(&right, &(SMOL_Matrix){.fields=(double*)vec3up, .nRows=3, .nCols=1});
    SMOL_CopyMatrix(&up, &right);
    SMOL_VectorCross(&up, &front);

    SMOL_AllocMatrix(lhs, 4, 4);
    
    SMOL_Scale(&front, -1.0);    
    SMOL_VectorNormalize(&up);
    SMOL_VectorNormalize(&front);
    SMOL_VectorNormalize(&right);

    SMOL_SetColumn(lhs, 0, &right);
    SMOL_SetColumn(lhs, 1, &up);
    SMOL_SetColumn(lhs, 2, &front);

    for(int i = 0; i < 3; i++)
	SMOL_SetField(lhs, i, 3, -vec3eye[i]);

    SMOL_FreeV(3, &front, &right, &up);

    return SMOL_STATUS_OK;
}

int SMOL_Fill(SMOL_Matrix *lhs, double value)
/* Fills the given matrix with the given value. */
{
    memset(lhs->fields, value, sizeof(double)*lhs->nRows*lhs->nCols);
    return SMOL_STATUS_OK;
}

int SMOL_Add(SMOL_Matrix* lhs, const SMOL_Matrix* rhs)
/* Add the rhs to the lhs; A = A + B. */
{
    if (lhs->nRows != rhs->nRows || lhs->nCols != rhs->nCols)
	return SMOL_STATUS_INCOMPATIBLE_SIZES;

    for(size_t r = 0; r < lhs->nRows; r++) {
	for (size_t c = 0; c < lhs->nCols; c++) 
	    lhs->fields[r*lhs->nCols+c] += rhs->fields[r*lhs->nCols+c];
    }

    return SMOL_STATUS_OK;
}

int SMOL_AddV(SMOL_Matrix *lhs, const SMOL_Matrix *rhs, ...)
/* Variadic addition. */
{
    va_list args;
    va_start(args, count);
    
    int status = 0;
    while (status == SMOL_STATUS_OK && count--)
	status = SMOL_Add(lhs, va_arg(args, SMOL_Matrix*));
    
    va_end(args);
    return status;
}

int SMOL_SubtractMat(SMOL_Matrix* lhs, const SMOL_Matrix* rhs)
/* Subtact the rhs from the lhs; A = A - B. */
{
    if (lhs->nRows != rhs->nRows || lhs->nCols != rhs->nCols)
	return SMOL_STATUS_INCOMPATIBLE_SIZES;

    for(size_t r = 0; r < lhs->nRows; r++) {
	for (size_t c = 0; c < lhs->nCols; c++) 
	    lhs->fields[r*lhs->nCols+c] -= rhs->fields[r*lhs->nCols+c];
    }

    return SMOL_STATUS_OK;
}

int SMOL_SubtractV(SMOL_Matrix *lhs, const SMOL_Matrix *rhs, ...)
/* Variadic subtraction. */
{
    va_list args;
    va_start(args, count);
    
    int status = 0;
    while (status == SMOL_STATUS_OK && count--)
	status = SMOL_Add(lhs, va_arg(args, SMOL_Matrix*));
    
    va_end(args);
    return status;
}

int SMOL_Multiply(SMOL_Matrix *lhs, const SMOL_Matrix *rhs)
/* Multiply the rhs of size NxK to the lhs of size MxN resulting in a MxK; A = A * B. */
{
    if (lhs->nCols != rhs->nRows)
	return SMOL_STATUS_INCOMPATIBLE_SIZES;

    SMOL_Matrix copy;
    SMOL_CopyMatrix(&copy, lhs);
    
    SMOL_Free(lhs);
    SMOL_AllocMatrix(lhs, copy.nRows, rhs->nCols);
    
    for (unsigned int rA = 0; rA < lhs->nRows; rA++) {
	for (unsigned int cB = 0; cB < lhs->nCols; cB++) {
	    for (unsigned int i = 0; i < copy.nCols; i++)
		lhs->fields[rA*rhs->nCols+cB] += copy.fields[rA*copy.nCols+i]*rhs->fields[i*rhs->nCols+cB];
	}
    }

    SMOL_Free(&copy);
    return SMOL_STATUS_OK;
}

int SMOL_MultiplyV(SMOL_Matrix *lhs, const SMOL_Matrix *rhs, ...)
/* Variadic multiplication. */
{
    va_list args;
    va_start(args, count);
    
    int status = 0;
    while (status == SMOL_STATUS_OK && count--)
	status = SMOL_Multiply(lhs, va_arg(args, SMOL_Matrix*));
    
    va_end(args);
    return status;
}

int SMOL_Scale(SMOL_Matrix *lhs, double scalar)
/* Multiply given Matrix with a scalar value. */
{
    for (unsigned int r = 0; r < lhs->nRows; r++) {
	for (unsigned int c = 0; c < lhs->nCols; c++)
	    lhs->fields[r*lhs->nCols+c] *= scalar;
    }

    return SMOL_STATUS_OK;
}

int SMOL_Transpose(SMOL_Matrix *lhs) 
/* Return the transposed matrix of the given matrix; A = A^T. */
{
    SMOL_Matrix copy;
    SMOL_CopyMatrix(&copy, lhs);
    
    SMOL_Free(lhs);
    lhs->nCols = copy.nRows;
    lhs->nRows = copy.nCols;
    
    for (unsigned int r = 0; r < lhs->nRows; r++) {
	for (unsigned int c = 0; c < lhs->nCols; c++)
	    lhs->fields[r*lhs->nCols+c] = copy.fields[c*copy.nCols+r];
    }

    SMOL_Free(&copy);
    return SMOL_STATUS_OK;
}

int SMOL_VectorNormalize(SMOL_Matrix *lhs)
/* Normalize the given vector. */
{
    if (SMOL_TypeOf(lhs) < SMOL_TYPE_VECTOR)
	return SMOL_STATUS_INVALID_TYPE;
    
    double l;
    SMOL_VectorLength(&l, lhs);
    for (size_t i = 0; i < lhs->nCols*lhs->nRows; i++)
	lhs->fields[i] /= l;

    return SMOL_STATUS_OK;
}

int SMOL_VectorCross(SMOL_Matrix *lhs, const SMOL_Matrix *rhs)
/* Assigns the cross product of the lhs and rhs vector the the lhs; A = A x B. */
{
    if (SMOL_TypeOf(lhs) < SMOL_TYPE_VECTOR ||
	SMOL_TypeOf(rhs) < SMOL_TYPE_VECTOR)
	return SMOL_STATUS_INVALID_TYPE;
    
    if (lhs->nCols != 1 || rhs->nCols != 1 || lhs->nRows != 3 || rhs->nRows != 3)
	return SMOL_STATUS_INCOMPATIBLE_SIZES;

    SMOL_Matrix copy;
    SMOL_CopyMatrix(&copy, lhs);
    lhs->fields[0] = copy.fields[1]*rhs->fields[2] - copy.fields[2]*rhs->fields[1];
    lhs->fields[1] = copy.fields[2]*rhs->fields[0] - copy.fields[0]*rhs->fields[2];
    lhs->fields[2] = copy.fields[0]*rhs->fields[1] - copy.fields[1]*rhs->fields[0];
    
    return SMOL_STATUS_OK;
}

int SMOL_VectorLength(double *lhs, const SMOL_Matrix *vec)
/* Return the length of the given vector. */
{
    if (SMOL_TypeOf(vec) < SMOL_TYPE_VECTOR)
	return SMOL_STATUS_INVALID_TYPE;
    
    *lhs = SMOL_VectorDot(lhs, vec, vec);
    *lhs = sqrt(*lhs);
    return SMOL_STATUS_OK;
}

int SMOL_VectorLengthSquare(double* lhs, const SMOL_Matrix *vec)
/* Return the length of the vector squared. */
{
    if (SMOL_TypeOf(vec) < SMOL_TYPE_VECTOR)
	return SMOL_STATUS_INVALID_TYPE;
    
    *lhs = SMOL_VectorDot(lhs, vec, vec);
    return SMOL_STATUS_OK;
}

int SMOL_VectorDot(double *lhs, const SMOL_Matrix *vecA, const SMOL_Matrix *vecB)
/* Return the dot product between vector A and B. */
{
    if (SMOL_TypeOf(vecA) < SMOL_TYPE_VECTOR ||
	SMOL_TypeOf(vecB) < SMOL_TYPE_VECTOR)
	return SMOL_STATUS_INVALID_TYPE;

    if (vecA->nCols != vecB->nCols || vecA->nRows != vecB->nRows)
	return SMOL_STATUS_INCOMPATIBLE_SIZES;
    
    *lhs = 0.0;
    for (size_t i = 0; i < vecA->nCols*vecA->nRows; i++)
	*lhs += vecA->fields[i] * vecB->fields[i];

    return SMOL_STATUS_OK;
}

int SMOL_SetField(SMOL_Matrix* lhs, size_t row, size_t col, double value)
/* Sets the value of given matrix at position [row, col]. */
{
    if (row > lhs->nRows || col > lhs->nCols)
	return SMOL_STATUS_ARRAY_OUT_OF_BOUNDS;
    
    lhs->fields[row*lhs->nCols + col] = value;

    return SMOL_STATUS_OK;
}

int SMOL_SetRow(SMOL_Matrix *lhs, size_t row, const SMOL_Matrix *vec)
/* Set the given row of the matrix to the given vector. */
{
    if (SMOL_TypeOf(vec) < SMOL_TYPE_VECTOR)
	return SMOL_STATUS_INCOMPATIBLE_SIZES;
    if (row > lhs->nRows || vec->nCols*vec->nRows > lhs->nCols)
	return SMOL_STATUS_ARRAY_OUT_OF_BOUNDS;
    
    memcpy(&lhs->fields[row*lhs->nCols], vec->fields, sizeof(double)*vec->nCols*vec->nRows);

    return SMOL_STATUS_OK;
}

int SMOL_SetColumn(SMOL_Matrix *lhs, size_t col, const SMOL_Matrix *vec)
/* Set the given column of the matrix to the given vector. */
{
    if (SMOL_TypeOf(vec) < SMOL_TYPE_VECTOR)
	return SMOL_STATUS_INCOMPATIBLE_SIZES;
    if (col > lhs->nCols || vec->nCols*vec->nRows > lhs->nRows)
	return SMOL_STATUS_ARRAY_OUT_OF_BOUNDS;
    
    for(size_t i = 0; i < vec->nCols*vec->nRows; i++)
	lhs->fields[i*lhs->nCols+col] = vec->fields[i];

    return SMOL_STATUS_OK;
}

int SMOL_GetField(double* lhs, const SMOL_Matrix* mat, size_t row, size_t col)
/* Returns the value of the matrix entry at [row, col]. */
{
    if (row > mat->nRows || col > mat->nCols)
	return SMOL_STATUS_ARRAY_OUT_OF_BOUNDS;
    
    *lhs = mat->fields[row*mat->nCols + col];

    return SMOL_STATUS_OK;
}

int SMOL_GetRow(SMOL_Matrix *lhs, const SMOL_Matrix *mat, size_t row)
/* Get the row vector of the given matrix and row. */
{
    if (row > mat->nRows)
	return SMOL_STATUS_ARRAY_OUT_OF_BOUNDS;
    
    SMOL_AllocMatrix(lhs, 1, mat->nCols);
    memcpy(lhs->fields, &mat->fields[row*mat->nCols], sizeof(double)*mat->nCols);

    return SMOL_STATUS_OK;
}

int SMOL_GetColumn(SMOL_Matrix *lhs, const SMOL_Matrix *mat, size_t col)
/* Get the column vector of the given matrix and column. */
{
    if (col > mat->nCols)
	return SMOL_STATUS_ARRAY_OUT_OF_BOUNDS;

    SMOL_AllocMatrix(lhs, mat->nRows, 1);
    for(size_t i = 0; i < mat->nRows; i++)
	lhs->fields[i] = mat->fields[i*mat->nCols+col];
    
    return SMOL_STATUS_OK;
}
