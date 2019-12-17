#include "minunit.h"
#include "matrix.h"

void setup ()
{
}

void teardown ()
{
}

MU_TEST(dotproduct) {
    SMOL_Matrix v1t, v2;
    SMOL_AllocMatrix(&v1t, 1, 2);
    SMOL_AllocMatrix(&v2, 2, 1);
    
    SMOL_SetField(&v1t, 0, 0, 0.5);
    SMOL_SetField(&v1t, 0, 1, 0.5);
    
    SMOL_SetField(&v2, 0, 0, -0.5);
    SMOL_SetField(&v2, 1, 0, 0.5);
    
    SMOL_Matrix r1;
    SMOL_CopyMatrix(&r1, &v1t);
    SMOL_Multiply(&r1, &v2);
    mu_assert_int_eq(1, r1.nRows);
    mu_assert_int_eq(1, r1.nCols);
    mu_assert_double_eq(r1.fields[0], 0.0);

    SMOL_FreeV(3, &v1t, &v2, &r1);
}

MU_TEST(addition) {
    SMOL_Matrix m1, m2, m3;
    SMOL_EyeMatrix(&m1, 4);
    SMOL_AllocMatrix(&m2, 4, 4);
    SMOL_AllocMatrix(&m3, 2, 4);

    for (int r = 0; r < 4; r++) {
	for (int c = 0; c < 4; c++)
	    SMOL_SetField(&m2, r, c, r+c);
    }
    // Sucessfull addition
    SMOL_Matrix r1;
    SMOL_CopyMatrix(&r1, &m1);
    SMOL_Add(&r1, &m2);
    for (int r = 0; r < 4; r++) {
	for (int c = 0; c < 4; c++) {
	    double ex = r + c;
	    if (r == c) 
		ex++;
	    double val = 0;
	    SMOL_GetField(&val, &r1, r, c);
	    mu_assert_double_eq(ex, val);
	}
    }

    SMOL_Matrix r2;
    SMOL_CopyMatrix(&r2, &m1);
    SMOL_Add(&r2, &m3);
    mu_check(r2.fields != NULL);

    SMOL_FreeV(5, &m1, &m2, &m3, &r1, &r2);
}

MU_TEST(multiplication) {
    double a1[] = { 2.0, -2.0, 1.0,
		    5.0, -7.0, 3.0};
    double a2[] = { 0.0, 5.0,
		    3.0, 7.0,
		    7.0, 4.0};
    SMOL_Matrix m1 = (SMOL_Matrix){.nRows=2,.nCols=3,.fields=a1};
    SMOL_Matrix m2 = (SMOL_Matrix){.nRows=3,.nCols=2,.fields=a2};

    SMOL_Matrix r1;
    SMOL_CopyMatrix(&r1, &m1);
    SMOL_Multiply(&r1, &m2);

    mu_assert_int_eq(2, r1.nRows);
    mu_assert_int_eq(2, r1.nCols);
    double f1, f2, f3, f4;
    SMOL_GetField(&f1, &r1, 0, 0);
    SMOL_GetField(&f2, &r1, 0, 1);
    SMOL_GetField(&f3, &r1, 1, 0);
    SMOL_GetField(&f4, &r1, 1, 1);
    mu_assert_double_eq(1.0,   f1);
    mu_assert_double_eq(0.0,   f2);
    mu_assert_double_eq(0.0,   f3);
    mu_assert_double_eq(-12.0, f4);

    SMOL_FreeV(1, &r1);
}

MU_TEST(crossproduct) {
    double a1[] = { 4.0, 5.0, 7.0};
    double a2[] = { 3.0, 2.0, 0.5};

    SMOL_Matrix m1 = (SMOL_Matrix){.nRows=3,.nCols=1,.fields=a1};
    SMOL_Matrix m2 = (SMOL_Matrix){.nRows=3,.nCols=1,.fields=a2};
    SMOL_Matrix m3;
    SMOL_CopyMatrix(&m3, &m1);
    SMOL_VectorCross(&m3, &m2);
    
    mu_assert(m3.fields != NULL, "");
    double f1, f2, f3;
    SMOL_GetField(&f1, &m3, 0, 0);
    SMOL_GetField(&f2, &m3, 1, 0);
    SMOL_GetField(&f3, &m3, 2, 0);
    mu_assert_double_eq(-11.5, f1);
    mu_assert_double_eq( 19,   f2);
    mu_assert_double_eq(-7,    f3);
    
    SMOL_Free(&m3);
}

MU_TEST_SUITE(matrix_suite) {
    MU_SUITE_CONFIGURE(&setup, &teardown);
    MU_RUN_TEST(addition);
    MU_RUN_TEST(multiplication);
    MU_RUN_TEST(dotproduct);
    MU_RUN_TEST(crossproduct);
}

int main ()
{
    MU_RUN_SUITE(matrix_suite);
    MU_REPORT();
    return 0;
}
