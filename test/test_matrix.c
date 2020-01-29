#include "minunit.h"
#include "matrix.h"

SMOL_Matrix eye;
SMOL_Matrix random;

void setup ()
{
}

void teardown ()
{
}

MU_TEST(matrix_multiplication)
{
    printf("multiplying");

    // Identity
    SMOL_Matrix result;
    SMOL_Multiply(&result, &random, &eye);
    for(size_t r = 0; r < random.nRows; r++) {
	for(size_t c = 0; c < random.nCols; c++) {
	    mu_assert_double_eq(random.fields[r*random.nCols+c], result.fields[r*result.nCols+c]);
	}
    }

    SMOL_Free(&result);

    // Dot product
    SMOL_Matrix vec1 = (SMOL_Matrix){.nCols=1, .nRows=3, .fields=(double[]){0.0, 0.0, 1.0}};
    SMOL_Matrix vec2 = (SMOL_Matrix){.nCols=3, .nRows=1, .fields=(double[]){1.0, 0.0, 0.0}};
    SMOL_Matrix d;
    SMOL_Multiply(&d, &vec2, &vec1);
    mu_assert_int_eq(1, d.nRows);
    mu_assert_int_eq(1, d.nCols);
    mu_assert_double_eq(0.0, d.fields[0]);

    SMOL_Free(&d);
    
    // Variadic multiplication
    SMOL_MultiplyV(&result, 4, &eye, &eye, &eye, &eye);

    for(size_t r = 0; r < result.nRows; r++) {
	for(size_t c = 0; c < result.nCols; c++) {
	    if(r == c)
		mu_assert_double_eq(1.0, result.fields[r*result.nCols+c]);
	    else
		mu_assert_double_eq(0.0, result.fields[r*result.nCols+c]);
	}
    }
    
    SMOL_Free(&result);
}

MU_TEST(matrix_addition)
{
    printf("adding");
    SMOL_Matrix res;
    SMOL_CopyMatrix(&res, &random);
    SMOL_AddV(&res, 4, &eye, &eye, &eye, &eye);

    for(size_t r = 0; r < res.nRows; r++) {
	for(size_t c = 0; c < res.nCols; c++) {
	    if(r == c)
		mu_assert_double_eq(random.fields[r*random.nCols+c]+4.0, res.fields[r*res.nCols+c]);
	}
    }

    SMOL_Free(&res);
}

MU_TEST(matrix_subtraction)
{
    printf("subtracting");
    SMOL_Matrix res;
    SMOL_CopyMatrix(&res, &random);
    SMOL_SubtractV(&res, 4, &eye, &eye, &eye, &eye);

    for(size_t r = 0; r < res.nRows; r++) {
	for(size_t c = 0; c < res.nCols; c++) {
	    if(r == c)
		mu_assert_double_eq(random.fields[r*random.nCols+c]-4.0, res.fields[r*res.nCols+c]);
	}
    }

    SMOL_Free(&res);
}

MU_TEST(matrix_invert)
{
    printf("inverting");
    SMOL_Matrix inv, res;
    SMOL_CopyMatrix(&inv, &random);
    SMOL_Invert(&inv);

    SMOL_Multiply(&res, &random, &inv);

    for(size_t r = 0; r < res.nRows; r++) {
	for(size_t c = 0; c < res.nCols; c++) {
	    if(r == c)
		mu_assert_double_eq(1.0, res.fields[r*res.nCols+c]);
	    else
		mu_assert_double_eq(0.0, res.fields[r*res.nCols+c]);
	}
    }
    
    SMOL_FreeV(2, &inv, &res);
}

MU_TEST_SUITE(matrix_suite)
{
    SMOL_EyeMatrix(&eye, 4);
    SMOL_RandomMatrix(&random, 4, 4, 10);
    
    MU_SUITE_CONFIGURE(&setup, &teardown);

    printf("Running matrix test suite..\n");
    
    printf("Todays random matrix:\n");
    SMOL_PrintMatrix(&random);

    MU_RUN_TEST(matrix_addition);
    MU_RUN_TEST(matrix_subtraction);
    MU_RUN_TEST(matrix_multiplication);
    MU_RUN_TEST(matrix_invert);

    SMOL_FreeV(2, &eye, &random);
}

int main ()
{
    MU_RUN_SUITE(matrix_suite);
    MU_REPORT();
    return 0;
}
