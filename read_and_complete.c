#include <stdio.h>
#include <time.h>
#include "maxdet_completion.h"
#define COUNT 1817
//COUT.max = 1817
int main()
{
	int dim = 25;
	clock_t begin,diff;
	int msec=0;
	int k = 0;
	double time_used=0;
	double temp;
	int empty_length = 0;
	int* empty_idx = NULL;
	double sum_zero = 0;
	int success = 0;
	gsl_matrix *M = gsl_matrix_alloc(dim, dim);
	FILE *f;
	char *matrix_dir = "OUT.new";
	f = fopen(matrix_dir, "r");
	if (!f) {
		fprintf(stderr, "File does not exist.\n");
		return 1;
	}
        for(int matrix_num = 0;matrix_num<COUNT;matrix_num++){
        	printf("Matrix Num: %d\n",matrix_num);
        	for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				if (fscanf(f, " %le", &(temp)) != 1) {
					fprintf(stderr, "Invalid data in %s.\n", matrix_dir);
					fclose(f);
					return EXIT_FAILURE;
				}
				if(i == j){
					temp = temp+1e-5;
					}
				gsl_matrix_set(M,i,j,temp);
			}
		}
		empty_length = 0;
		for (int i = 0; i < dim; i++) {
		for (int j = i; j < dim; j++) {
				if(gsl_isnan(gsl_matrix_get(M,i,j))){
						empty_length = empty_length + 1;
					}
			}
		}
		k = 0;
		empty_idx = malloc(empty_length*2*sizeof(int));
		for (int i = 0; i < dim; i++) {
		for (int j = i; j < dim; j++) {
		if(gsl_isnan(gsl_matrix_get(M,i,j))){
				empty_idx[k] = i;
				empty_idx[k+1] = j;
				k = k+2;
					}
				
			}
		}
		
		begin = clock();
		matrix_completion(M);
		diff = clock() - begin;
		
		msec = msec + diff*1000/CLOCKS_PER_SEC;
		time_used = time_used + ((double) diff)/CLOCKS_PER_SEC;
		
		sum_zero = 0;
		gsl_set_error_handler_off();
		if(gsl_linalg_cholesky_decomp1(M) == 0){
			gsl_linalg_cholesky_invert(M);
			for(int i = 0;i<empty_length;i++){
				sum_zero = sum_zero + fabs(2*gsl_matrix_get(M,empty_idx[2*i],empty_idx[2*i+1]));
			}
			printf("The sum of %d fillin in Minv matrix is:%f.\n", empty_length,sum_zero);
			success = success +1;
		}
		free(empty_idx);
        }
	printf("Time spent for %d completion is %d seconds %d milliseconds.\n",success,msec/1000,msec%1000);
	gsl_matrix_free(M);
	fclose(f);
        
	return 0;
}









/*
empty_idx = get_empty_idx(M, &empty_length);
		for (int i = 0; i < empty_length; i++) {
			gsl_matrix_set(M,empty_idx[2*i],empty_idx[2*i+1],0.9);
			gsl_matrix_set(M,empty_idx[2*i+1],empty_idx[2*i],0.9);
		}
		
		make_psd(M,workspace,empty_length,empty_idx);
		
		gsl_matrix_memcpy(workspace, M);
			if (gsl_linalg_cholesky_decomp1(workspace) == 0) {
				optim_newton(M, workspace, empty_length, empty_idx,100);
			} 


printf("flag\n");


printf("The M matrix is:\n");
		for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
				printf("%12.5f ", gsl_matrix_get(M, i, j));
			}
			printf("\n");
		}
		printf("\n");

*/

/*
	printf("The M matrix is:\n");
		for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
				printf("%12.5f ", gsl_matrix_get(M, i, j));
			}
			printf("\n");
		}
		printf("\n");
		
	printf("The M_remained matrix is:\n");
		for (int i = 0; i < dim_remained; i++) {
		for (int j = 0; j < dim_remained; j++) {
				printf("%12.5f ", gsl_matrix_get(M_remained, i, j));
			}
			printf("\n");
		}
		printf("\n");
		
	printf("The M_recovered matrix is:\n");
		for (int i = 0; i < dim_remained; i++) {
		for (int j = 0; j < dim_remained; j++) {
				printf("%12.5f ", gsl_matrix_get(M_remained, i, j) + gsl_matrix_get(M_reducer, i, j));
			}
			printf("\n");
		}
		printf("\n");
		
		printf("The adj is:\n");
		for (int i = 0; i < dim_remained; i++) {
		for (int j = 0; j < dim_remained; j++) {
				printf("%d ", adj[i][j]);
			}
			printf("\n");
		}
		printf("\n");
		
		
			printf("The set is:\n");
		for (int i = 0; i < dim_remained; i++) {
		for (int j = 0; j < dim_remained; j++) {
				printf("%d ", set[i][j]);
			}
			printf("\n");
		}
		printf("\n");
		
	printf("The size is:\n");
		for (int i = 0; i < dim_remained; i++) {
		printf("%d ", size[i]);
			
		}
		printf("\n");
		
		printf("The alpha is:\n");
		for (int i = 0; i < dim_remained; i++) {
		printf("%d ", alpha[i]);
		}
		printf("\n");
		
		printf("The alpha_inv is:\n");
		for (int i = 0; i < dim_remained; i++) {
		printf("%d ", alpha_inv[i]);
		}
		printf("\n");
		
		
		
		

	*/






