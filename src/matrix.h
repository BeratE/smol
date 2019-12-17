#ifndef SMOL_MATRIX_H
#define SMOL_MATRIX_H
/* Simple Matrix Operation Library (smol). Small collection of handy matrix operations. */

#include <stdarg.h>
#include <stdlib.h>

typedef struct SMOL_Matrix {
    size_t nRows;
    size_t nCols;
    double *fields;
} SMOL_Matrix;

/* Enums */
enum SMOL_TYPE {SMOL_TYPE_NULL,
		SMOL_TYPE_MATRIX,
		SMOL_TYPE_VECTOR,
		SMOL_TYPE_ROW_VECTOR,
		SMOL_TYPE_COLUMN_VECTOR,
		SMOL_TYPE_SCALAR,};

enum SMOL_STATUS {SMOL_STATUS_OK,
		  SMOL_STATUS_INVALID_TYPE,
		  SMOL_STATUS_INCOMPATIBLE_SIZES,
		  SMOL_STATUS_ARRAY_OUT_OF_BOUNDS};

/* Macros */ 
#define SMOL_NULLMATRIX (SMOL_Matrix){.nRows = 0, .nCols = 0, .fields = NULL}

void SMOL_PrintError(enum SMOL_STATUS status);
void SMOL_Free(SMOL_Matrix* mat);
void SMOL_FreeV(int count, ...);
int SMOL_TypeOf(const SMOL_Matrix *mat);

int SMOL_AllocMatrix(SMOL_Matrix* lhs, size_t nRows, size_t nCols);
int SMOL_RandomMatrix(SMOL_Matrix *lhs, size_t nRows, size_t nCols);
int SMOL_CopyMatrix(SMOL_Matrix *lhs, const SMOL_Matrix *rhs);
int SMOL_EyeMatrix(SMOL_Matrix *lhs, size_t size);
int SMOL_PerspectiveMatrix(SMOL_Matrix* lhs, double fov, double ratio, double near, double far);
int SMOL_CameraMatrix(SMOL_Matrix* lhs, const double* vec3eye, const double* vec3front, const double* vec3up);

int SMOL_Fill(SMOL_Matrix *lhs, double value);
int SMOL_Add(SMOL_Matrix *lhs, const SMOL_Matrix *rhs);
int SMOL_Subtract(SMOL_Matrix *lhs, const SMOL_Matrix *rhs);
int SMOL_Multiply(SMOL_Matrix *lhs, const SMOL_Matrix *rhs);
int SMOL_Scale(SMOL_Matrix *lhs, double scalar);
//int SMOL_Rotate(SMOL_Matrix *lhs, const double* vec3axis, double angle);
int SMOL_Transpose(SMOL_Matrix *lhs);

int SMOL_VectorNormalize(SMOL_Matrix *lhs);
int SMOL_VectorCross(SMOL_Matrix *lhs, const SMOL_Matrix *rhs);
int SMOL_VectorLength(double* lhs, const SMOL_Matrix *vec);
int SMOL_VectorLentghSquare(double *lhs, const SMOL_Matrix *vec);
int SMOL_VectorDot(double* lhs, const SMOL_Matrix *vecA, const SMOL_Matrix *vecB);

int SMOL_SetField(SMOL_Matrix* lhs, size_t row, size_t col, double value);
int SMOL_SetRow(SMOL_Matrix *lhs, size_t row, const SMOL_Matrix *vec);
int SMOL_SetColumn(SMOL_Matrix *lhs, size_t col, const SMOL_Matrix *vec);
int SMOL_GetField(double *lhs, const SMOL_Matrix* mat, size_t row, size_t col);
int SMOL_GetRow(SMOL_Matrix *lhs, const SMOL_Matrix *mat, size_t row);
int SMOL_GetColumn(SMOL_Matrix *lhs, const SMOL_Matrix *mat, size_t col);

#endif // SMOL_MATRIX_H 
