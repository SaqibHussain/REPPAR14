/***************************************************
File:			MatrixMult.c
Description:		MPI Matrix Multiplication
Author:			Saqib Hussain
Organisation:		School of Engineering, Cranfield University
Email:			s.m.hussain@cranfield.ac.uk
Copyright:		Copyright Saqib Hussain 2014
Note:			The following code is a modification of that found 
			online @ http://preeyakorn.rmutl.ac.th/mpi/matrixmultiply_blanks.c
*****************************************************/

#define _CRT_SECURE_NO_DEPRECATE

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define TAG_MATRIX_PARTITION	0x4560

typedef struct
{  unsigned int   m, n;  // Rows, cols
double        *data;  // Data, ordered by row, then by col
double       **rows;  // Pointers to rows in data
} TMatrix;

TMatrix createMatrix     (const unsigned int rows, const unsigned int cols);
void    freeMatrix       (TMatrix *matrix);
int     validMatrix      (TMatrix matrix);
TMatrix initMatrix       (void);
TMatrix matrixMultiply   (TMatrix A, TMatrix B);
void    doMatrixMultiply (TMatrix A, TMatrix B, TMatrix C);
void    printMatrix      (char name[128],TMatrix A);
int readMatrix(char *filename, TMatrix *A);

int main (int argc, char *argv[])
{
	int processor_rank  = 0;
	int processor_count = 1;
	MPI_Status   status;
	TMatrix      A,B,C,D;
	unsigned int m, n, i, j;
	double       time0, time1;
	FILE *file;
	A = initMatrix(); B = initMatrix(); C = initMatrix(); D = initMatrix();
	MPI_Init(&argc, &argv); 
	MPI_Comm_size (MPI_COMM_WORLD, &processor_count);
	MPI_Comm_rank (MPI_COMM_WORLD, &processor_rank );
	if (processor_rank == 0)
	{  
		file = fopen(argv[3], "w");
		if (file == NULL) { printf("Error"); exit(1); }

		readMatrix(argv[1], &A);
		readMatrix(argv[2], &B);

		time0 = MPI_Wtime();
		n = A.n;
		m = n / processor_count;		
		C = createMatrix(n,n);

		// Broadcast (send) size of matrix
		MPI_Bcast((void *)&n, 1, MPI_INT, 0, MPI_COMM_WORLD); 

		// Broadcast (send) all of B matrix
		MPI_Bcast((void *)B.data, n*n,MPI_DOUBLE, 0, MPI_COMM_WORLD);

		// Send each process it's own part of A
		for (i = 1; i < processor_count; i++)
			MPI_Send((void *)A.rows[i*m], m*n, MPI_DOUBLE, i, TAG_MATRIX_PARTITION, MPI_COMM_WORLD);

		// Multiply own part of matrix A with B into already existing matrix C
		A.m = m;
		doMatrixMultiply(A,B,C);

		A.m = n;

		// Receive part of C matrix from each process
		for (i = 1; i < processor_count; i++)
			MPI_Recv((void *)C.rows[m*i], m*n, MPI_DOUBLE,i, TAG_MATRIX_PARTITION,MPI_COMM_WORLD, &status);

		// Record finish time
		time1 = MPI_Wtime();
		//printMatrix("A",A);
		//printMatrix("B",B);
		//printMatrix("C",C);
		// Print time statistics
		printf ("Total time using [%2d] processors : [%f] seconds\n", processor_count, time1 - time0);
		fprintf(file, "%f Seconds\n\n", time1 - time0);
		//for (i = 0; i < C.m; i++)
		//{  
		//	for (j = 0; j < C.n; j++) 
		//		fprintf (file, "%7.3f ", C.rows[i][j]);
		//	fprintf (file, "\n");
		//}
		fclose(file);
	}
	else
	{
		// Broadcast (receive) size of matrix
		MPI_Bcast((void *)&n, 1, MPI_INT, 0, MPI_COMM_WORLD); 

		// Allocate memory for matrices
		m = n / processor_count;
		A = createMatrix(m, n);
		B = createMatrix(n ,n);

		// Broadcast (receive) B matrix
		MPI_Bcast((void *)B.data, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);


		// Receive its part of A matrix
		MPI_Recv((void *)A.data, m*n, MPI_DOUBLE, 0, TAG_MATRIX_PARTITION,MPI_COMM_WORLD, &status);

		// Multiply local matrices
		C = matrixMultiply(A,B);

		// Send result back to 0
		MPI_Send((void *)C.data, m*n, MPI_DOUBLE, 0, TAG_MATRIX_PARTITION, MPI_COMM_WORLD);
	}

	// Free matrix data
	freeMatrix(&A); freeMatrix(&B); freeMatrix(&C);

	// Wait for everyone to stop   
	MPI_Barrier(MPI_COMM_WORLD);

	// Always use MPI_Finalize as the last instruction of the program
	MPI_Finalize();
	return 0;
}

TMatrix createMatrix(const unsigned int rows, const unsigned int cols)
{  
	TMatrix           matrix;
	unsigned long int m, n;
	unsigned int      i,j;

	m = rows; n = cols;
	matrix.m    = rows;
	matrix.n    = cols;
	matrix.data = (double *) malloc(sizeof(double) * m * n);
	matrix.rows = (double **) malloc(sizeof(double *) * m);

	if (validMatrix(matrix))
	{  
		matrix.m = rows; 
		matrix.n = cols;
		for (i = 0; i < rows; i++)
		{  
			matrix.rows[i] = matrix.data + (i * cols);
		}

	}
	else
	{  
		freeMatrix(&matrix);
	}
	return matrix;
}

void freeMatrix (TMatrix *matrix)
{  
	if (matrix == NULL) return;
	if (matrix -> data) { free(matrix -> data); matrix -> data = NULL; }
	if (matrix -> rows) { free(matrix -> rows); matrix -> rows = NULL; }
	matrix -> m = 0;
	matrix -> n = 0;
}

int validMatrix (TMatrix matrix)
{  
	if ((matrix.data == NULL) || (matrix.rows == NULL) ||
		(matrix.m == 0) || (matrix.n == 0))
		return 0;
	else return 1;
}

TMatrix initMatrix()
{  
	TMatrix matrix;
	matrix.m = 0;
	matrix.n = 0;
	matrix.data = NULL;
	matrix.rows = NULL;
	return matrix;
}

TMatrix matrixMultiply(TMatrix A, TMatrix B)
{  
	TMatrix C;
	C = initMatrix();
	if (validMatrix(A) && validMatrix(B) && (A.n == B.m))
	{  
		C = createMatrix(A.m, B.n);
		if (validMatrix(C))
		{  
			doMatrixMultiply(A, B, C);
		}
	}
	return C;
}

void doMatrixMultiply(TMatrix A, TMatrix B, TMatrix C)
{  
	unsigned int i, j, k;
	double sum;
	for (i = 0; i < A.m; i++) // Rows
	{  
		for (j = 0; j < B.n; j++) // Cols
		{  
			sum = 0;
			for (k = 0; k < A.n; k++) 
				sum += A.rows[i][k] * B.rows[k][j];
			C.rows[i][j] = sum;
		}
	}
}

void printMatrix(char name[128], TMatrix A)
{  
	unsigned int i, j;
	printf("%s:\n", name);
	if (validMatrix(A))
	{  
		for (i = 0; i < A.m; i++)
		{  
			for (j = 0; j < A.n; j++) 
				printf ("%7.3f ", A.rows[i][j]);
			printf ("\n");
		}
	}
}

int readMatrix(char *filename, TMatrix *A)
{  FILE *fp;
unsigned int m, n, i, j;
float d;
int result = 0;

if ((fp = fopen (filename, "r")) == NULL) return 0;

do
{
	if (fscanf (fp, "%d%d",  &m, &n) != 2) break;
	if ((m == 0) || (n == 0)) break;
	*A = createMatrix(m,n);
	if (!validMatrix(*A)) break;

	for (i = 0; i < m; i ++)
	{  for (j = 0; j < n; j ++)
	{  if (fscanf (fp, "%f", &d) != 1) break;
	A -> rows[i][j] = d;
	}
	if (j != n) break;
	}
	if (i != m) break;

	result = 1;
} while (0);

fclose (fp);

return result;
}
