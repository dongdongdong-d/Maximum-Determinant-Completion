#include <math.h>
#include <string.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
void completion_one(gsl_matrix * M_target, int col, int row, int *nei, int num_nei)
{
	if (num_nei == 0) {
		gsl_matrix_set(M_target, col, row, 0);
		gsl_matrix_set(M_target, row, col, 0);
	} else {
		int i, j;
		double temp;
		gsl_matrix *C = gsl_matrix_alloc(num_nei, num_nei);
		gsl_vector *B = gsl_vector_alloc(num_nei);
		gsl_vector *D = gsl_vector_alloc(num_nei);
		for (i = 0; i < num_nei; i++) {
			gsl_vector_set(B, i, gsl_matrix_get(M_target, col, nei[i]));
			gsl_vector_set(D, i, gsl_matrix_get(M_target, row, nei[i]));

			for (j = i; j < num_nei; j++) {
				gsl_matrix_set(C, i, j, gsl_matrix_get(M_target, nei[j], nei[i]));
				gsl_matrix_set(C, j, i, gsl_matrix_get(M_target, nei[j], nei[i]));
			}
		}


		gsl_linalg_HH_svx(C, B);
		gsl_blas_ddot(B, D, &temp);
		gsl_matrix_set(M_target, col, row, temp);
		gsl_matrix_set(M_target, row, col, temp);


		gsl_matrix_free(C);
		gsl_vector_free(B);
		gsl_vector_free(D);
	}
}

void schur_complement(gsl_matrix * M, gsl_matrix * M_remained, gsl_matrix * M_reducer, int *reduced_idx, int *full_idx, int dim_remained,
		      int dim_full)
{
	int col, row, i, j;
	double temp;
	gsl_matrix *C = gsl_matrix_alloc(dim_full, dim_full);
	gsl_matrix *B = gsl_matrix_alloc(dim_remained, dim_full);
	gsl_matrix *BT = gsl_matrix_alloc(dim_full, dim_remained);
	for (i = 0; i < dim_remained; i++) {
		for (j = 0; j < dim_full; j++) {
			row = reduced_idx[i];
			col = full_idx[j];
			temp = gsl_matrix_get(M, row, col);
			gsl_matrix_set(B, i, j, temp);
			gsl_matrix_set(BT, j, i, temp);
		}
	}

	for (i = 0; i < dim_full; i++) {
		for (j = 0; j < dim_full; j++) {
			row = full_idx[i];
			col = full_idx[j];
			temp = gsl_matrix_get(M, row, col);
			gsl_matrix_set(C, i, j, temp);
		}
	}

	gsl_linalg_cholesky_decomp1(C);
	gsl_blas_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1, C, BT);
	gsl_blas_dtrsm(CblasLeft, CblasLower, CblasTrans, CblasNonUnit, 1, C, BT);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, B, BT, 0, M_reducer);

	for (i = 0; i < dim_remained; i++) {
		for (j = i; j < dim_remained; j++) {
			row = reduced_idx[i];
			col = reduced_idx[j];
			temp = gsl_matrix_get(M, row, col) - gsl_matrix_get(M_reducer, i, j);
			gsl_matrix_set(M_remained, i, j, temp);
			gsl_matrix_set(M_remained, j, i, temp);
		}
	}

	gsl_matrix_free(C);
	gsl_matrix_free(B);
	gsl_matrix_free(BT);
}

void MCS_seach_chordal_check(int *alpha, int *alpha_inv, int **adj, int dim_remained, int *num_fill_in, int *fill_in)
{
	int **set = malloc(sizeof(int *) * dim_remained);
	int *size = malloc(sizeof(int) * dim_remained);
	int *f = malloc(sizeof(int) * dim_remained);
	int *idx = malloc(sizeof(int) * dim_remained);
	int v, w, x, current_order, current_weight_level = 0, is_level_active;
	for (v = 0; v < dim_remained; v++) {
		set[v] = malloc(sizeof(int) * dim_remained);
		size[v] = 0;
		f[v] = 0;
		idx[v] = 0;
		for (w = 0; w < dim_remained; w++) {
			set[v][w] = 0;
		}
		set[0][v] = 1;
	}

	for (current_order = dim_remained - 1; current_order >= 0; current_order--) {
		for (v = 0; v < dim_remained; v++) {
			if (set[current_weight_level][v]) {
				break;
			}
		}
		set[current_weight_level][v] = 0;
		alpha[v] = current_order;
		alpha_inv[current_order] = v;
		size[v] = -1;
		for (w = 0; w < dim_remained; w++) {
			if (adj[v][w] && size[w] >= 0) {
				set[size[w]][w] = 0;
				size[w] = size[w] + 1;
				set[size[w]][w] = 1;
			}
		}
		current_weight_level = current_weight_level + 1;
		while (current_weight_level > 0) {
			is_level_active = 0;
			for (w = 0; w < dim_remained; w++) {
				if (set[current_weight_level][w]) {
					is_level_active = 1;
					break;
				}
			}
			if (is_level_active) {
				break;
			} else {
				current_weight_level = current_weight_level - 1;
			}
		}
	}


	*num_fill_in = 0;
	for (current_order = 0; current_order < dim_remained; current_order++) {
		w = alpha_inv[current_order];
		f[w] = w;
		idx[w] = current_order;
		for (v = 0; v < dim_remained; v++) {
			if (adj[v][w] == 1 && alpha[v] < current_order) {
				x = v;
				while (idx[x] < current_order) {
					idx[x] = current_order;
					if (adj[x][w] == 0) {
						fill_in[*num_fill_in] = x;
						fill_in[*num_fill_in + 1] = w;
						*num_fill_in = *num_fill_in + 2;
					}
					x = f[x];
				}
				if (f[x] == x) {
					f[x] = w;
				}
			}
		}
	}

	*num_fill_in = (*num_fill_in + 1) / 2;


	free(f);
	free(idx);
	for (v = 0; v < dim_remained; v++) {
		free(set[v]);
	}
	free(set);
	free(size);
}

//This method seems not good.
void BK(int *R, int *P, int *X, int **adj, int length, int **R_res, int *count)
{
	int R_sum = 0, P_sum = 0, X_sum = 0;
	int i, v, u, u_nbs, u_pivot = 0, u_pivot_nbs;
	int *R_new = NULL;
	int *P_new = NULL;
	int *X_new = NULL;
	for (int i = 0; i < length; i++) {
		R_sum += R[i];
		P_sum += P[i];
		X_sum += X[i];
	}

	if (P_sum == 0 && X_sum == 0) {
		R_res[*count] = R;
		*count = *count + 1;
		//if (P) free(P);
		//if (X) free(X);
	} else {
		u_pivot_nbs = -1;
		for (u = 0; u < length; u++) {
			if (P[u] == 1 || X[u] == 1) {
				u_nbs = 0;
				for (i = 0; i < length; i++) {
					u_nbs = adj[u][i] + u_nbs;
				}
				if (u_nbs > u_pivot_nbs) {
					u_pivot = u;
					u_pivot_nbs = u_nbs;
				}
			}
		}
		for (v = 0; v < length; v++) {
			if (P[v] == 1 && adj[u_pivot][v] == 0) {
				R_new = (int *) calloc(3 * length, sizeof(int));
				P_new = R_new + length;
				X_new = R_new + 2 * length;
				memcpy(R_new, R, length*sizeof(int));
				for (i = 0; i < length; i++) {
					//R_new[i] = R[i];
					P_new[i] = (P[i] == 1 && adj[v][i] == 1);
					X_new[i] = (X[i] == 1 && adj[v][i] == 1);
				}
				R_new[v] = 1;
				BK(R_new, P_new, X_new, adj, length, R_res, count);
				P[v] = 0;
				X[v] = 1;
				//free(R_new);
			}
		}
		free(R); R = NULL;
	}
}
void guess_fill_in(gsl_matrix * M_target, int col, int row, int *nei, int **adj)
{
	int i, j;
	int num_nei = 0;
	int dim = M_target->size1;
	for (i = 0; i < dim; i++) {
		if (adj[row][i] && adj[col][i]) {
			nei[num_nei] = i;
			num_nei = num_nei + 1;
		}
	}

	int **adj_sub = malloc(sizeof(int *) * num_nei);
	int **R_res = malloc(sizeof(int **) * 100);	       // It should be changed to a linked list.
	int R_count = 0;
	int *R = malloc(sizeof(int) * num_nei);
	int *P = malloc(sizeof(int) * num_nei);
	int *X = malloc(sizeof(int) * num_nei);
	for (i = 0; i < num_nei; i++) {
		R[i] = 0;
		P[i] = 1;
		X[i] = 0;
		adj_sub[i] = malloc(sizeof(int) * num_nei);
		for (j = 0; j < num_nei; j++) {
			adj_sub[i][j] = adj[nei[i]][nei[j]];
			if (j == i) {
				adj_sub[i][j] = 0;
			}
		}
	}
	BK(R, P, X, adj_sub, num_nei, R_res, &R_count);
	// ##### linear algebra part starts
	int dim_clique, current_clique, dim_sub, signum;
	int *current_clique_idx;
	gsl_matrix *sub;
	// gsl_matrix* sub_inv;
	gsl_permutation *p;
	gsl_vector *e;
	double temp, b, a, delta, right = 1000, left = -1000, right_new, left_new, inv_11, inv_22, inv_12;
	for (current_clique = 0; current_clique < R_count; current_clique++) {
		dim_clique = 0;
		for (i = 0; i < num_nei; i++) {
			if (R_res[current_clique][i] == 1) {
				dim_clique = dim_clique + 1;
			}
		}
		dim_sub = dim_clique + 2;
		current_clique_idx = malloc(dim_sub * sizeof(int));
		j = 0;
		for (i = 0; i < num_nei; i++) {
			if (R_res[current_clique][i] == 1) {
				current_clique_idx[j] = nei[i];
				j++;
			}
		}
		current_clique_idx[dim_clique] = row;
		current_clique_idx[dim_clique + 1] = col;
		sub = gsl_matrix_alloc(dim_sub, dim_sub);
		e = gsl_vector_alloc(dim_sub);
		p = gsl_permutation_alloc(dim_sub);
		for (i = 0; i < dim_sub; i++) {
			for (j = i; j < dim_sub; j++) {
				temp = gsl_matrix_get(M_target, current_clique_idx[i], current_clique_idx[j]);
				gsl_matrix_set(sub, i, j, temp);
				gsl_matrix_set(sub, j, i, temp);
			}
		}
		gsl_matrix_set(sub, dim_clique, dim_clique + 1, 0);
		gsl_matrix_set(sub, dim_clique + 1, dim_clique, 0);

		gsl_linalg_LU_decomp(sub, p, &signum);
		gsl_vector_set_basis(e, dim_clique);
		gsl_linalg_LU_svx(sub, p, e);
		inv_11 = gsl_vector_get(e, dim_clique);
		inv_12 = gsl_vector_get(e, dim_clique + 1);
		gsl_vector_set_basis(e, dim_clique + 1);
		gsl_linalg_LU_svx(sub, p, e);
		inv_22 = gsl_vector_get(e, dim_clique + 1);
		// gsl_linalg_LU_invert(sub,p,sub_inv);
		b = 2 * inv_12;
		a = inv_12 * inv_12 - inv_11 * inv_22;
		delta = 4 * inv_11 * inv_22;

		if (delta > 0) {
			if (a < 0) {
				left_new = (-b + sqrt(delta)) / a / 2;
				right_new = (-b - sqrt(delta)) / a / 2;
				if (left_new > left) {
					left = left_new;
				}
				if (right_new < right) {
					right = right_new;
				}
			} else {
				left_new = (-b - sqrt(delta)) / a / 2;
				right_new = (-b + sqrt(delta)) / a / 2;
				if (left_new > left) {
					left = left_new;
				}
				if (right_new < right) {
					right = right_new;
				}
			}

		}

		free(current_clique_idx);
		gsl_matrix_free(sub);
		gsl_vector_free(e);
		gsl_permutation_free(p);
	}
	temp = 0.5 * (left + right);
	gsl_matrix_set(M_target, col, row, temp);
	gsl_matrix_set(M_target, row, col, temp);

	// ################

	for (i = 0; i < R_count; i++) {
		free(R_res[i]);
	}
	free(R_res);
	for (i = 0; i < num_nei; i++) {
		free(adj_sub[i]);
	}
	free(adj_sub);
}

void guess_fill_in2(gsl_matrix * M_target, int col, int row, int *nei, int *nnbs, int **adj)
{
	int i, j;
	int num_nei = 0;
	int dim = M_target->size1;
	for (i = 0; i < dim; i++) {
		if (adj[row][i] && adj[col][i]) {
			nei[num_nei] = i;
			num_nei = num_nei + 1;
		}
	}

	int minimum_nbs = num_nei;
	int minimum_nbs_idx = num_nei;

	for (i = 0; i < num_nei; i++) {
		nnbs[i] = 0;
		for (j = 0; j < num_nei; j++) {
			nnbs[i] = nnbs[i] + adj[nei[i]][nei[j]];
		}
		if (minimum_nbs > nnbs[i]) {
			minimum_nbs = nnbs[i];
			minimum_nbs_idx = i;
		}
	}

	while (minimum_nbs < num_nei) {
		for (j = minimum_nbs_idx; j < num_nei; j++) {
			nei[j] = nei[j + 1];
		}
		num_nei = num_nei - 1;
		minimum_nbs = num_nei;
		for (i = 0; i < num_nei; i++) {
			nnbs[i] = 0;
			for (j = 0; j < num_nei; j++) {
				nnbs[i] = nnbs[i] + adj[nei[i]][nei[j]];
			}
			if (minimum_nbs > nnbs[i]) {
				minimum_nbs = nnbs[i];
				minimum_nbs_idx = i;
			}
		}
	}

	completion_one(M_target, row, col, nei, num_nei);
}


void optim_newton(gsl_matrix * M, gsl_matrix * workspace, int empty_length, int *empty_idx, int iter_max, double tol)
{
	gsl_vector *df = gsl_vector_alloc(empty_length);
	gsl_matrix *ddf = gsl_matrix_alloc(empty_length, empty_length);
	double x_temp;
	double df_norm_now = tol + 1;
	int i, j, coli, rowi, colj, rowj, iter = 0;
	while (df_norm_now > tol) {
		gsl_linalg_cholesky_invert(workspace);
		for (i = 0; i < empty_length; i++) {
			rowi = empty_idx[2 * i];
			coli = empty_idx[2 * i + 1];
			x_temp = 2 * gsl_matrix_get(workspace, rowi, coli);
			gsl_vector_set(df, i, x_temp);
			for (j = i; j < empty_length; j++) {
				rowj = empty_idx[2 * j];
				colj = empty_idx[2 * j + 1];
				x_temp =
				    2 * (gsl_matrix_get(workspace, rowi, rowj) * gsl_matrix_get(workspace, coli, colj) +
					 gsl_matrix_get(workspace, rowi, colj) * gsl_matrix_get(workspace, coli, rowj));
				gsl_matrix_set(ddf, i, j, x_temp);
				gsl_matrix_set(ddf, j, i, x_temp);
			}
		}
		df_norm_now = fabs(gsl_blas_dnrm2(df));
		// printf("abs(df) = %f\n",df_norm_now);
		gsl_linalg_cholesky_decomp1(ddf);
		gsl_linalg_cholesky_svx(ddf, df);
		for (i = 0; i < empty_length; i++) {
			rowi = empty_idx[2 * i];
			coli = empty_idx[2 * i + 1];
			x_temp = gsl_matrix_get(M, rowi, coli) + gsl_vector_get(df, i);
			gsl_matrix_set(M, rowi, coli, x_temp);
			gsl_matrix_set(M, coli, rowi, x_temp);
		}
		iter = iter + 1;
		if (iter > iter_max) {
			break;
		}
		gsl_matrix_memcpy(workspace, M);
		gsl_linalg_cholesky_decomp1(workspace);
	}
	gsl_vector_free(df);
	gsl_matrix_free(ddf);
}



void make_psd(gsl_matrix * M, gsl_matrix * workspace, int empty_length, int *empty_idx, int iter_max)
{
	int dim = M->size1, iter_now = 0;
	double shift, shift_up_cum, shift_down_cum = 0;
	gsl_eigen_symm_workspace *eigen_workspace = gsl_eigen_symm_alloc(dim);
	gsl_vector *eigen_values = gsl_vector_alloc(dim);
	double *diag_origin = malloc(dim * sizeof(double));
	gsl_matrix_memcpy(workspace, M);
	gsl_eigen_symm(workspace, eigen_values, eigen_workspace);
	gsl_sort_vector(eigen_values);
	shift_up_cum = -gsl_vector_get(eigen_values, 0) + 1e-5;
	shift = shift_up_cum;
	for (int i = 0; i < dim; i++) {
		diag_origin[i] = gsl_matrix_get(M, i, i);
	}
	while (shift_down_cum < shift_up_cum) {
		for (int i = 0; i < dim; i++) {
			gsl_matrix_set(M, i, i, gsl_matrix_get(M, i, i) + shift);
		}

		gsl_matrix_memcpy(workspace, M);
		gsl_linalg_cholesky_decomp1(workspace);
		optim_newton(M, workspace, empty_length, empty_idx, 100, 1e-5);
		gsl_matrix_memcpy(workspace, M);
		gsl_eigen_symm(workspace, eigen_values, eigen_workspace);
		gsl_sort_vector(eigen_values);
		shift_down_cum = shift_down_cum + gsl_vector_get(eigen_values, 0) * 0.99;
		shift = -gsl_vector_get(eigen_values, 0) * 0.99;
		iter_now = iter_now + 1;
		if (iter_now > iter_max) {
			break;
		}
	}

	for (int i = 0; i < dim; i++) {
		gsl_matrix_set(M, i, i, diag_origin[i]);
	}
	gsl_matrix_memcpy(workspace, M);
	gsl_eigen_symm(workspace, eigen_values, eigen_workspace);
	gsl_sort_vector(eigen_values);

	free(diag_origin);
	gsl_vector_free(eigen_values);
	gsl_eigen_symm_free(eigen_workspace);

}

void chordal_completion(gsl_matrix * M_remained, int dim_remained, int *alpha_inv, int *nei, int **adj)
{
	int col, row, i, j, k, ii;
	for (i = dim_remained - 1; i >= 0; i--) {
		col = alpha_inv[i];
		for (j = dim_remained - 1; j >= 0; j--) {
			row = alpha_inv[j];
			if (adj[col][row] == 0) {
				k = 0;
				for (ii = 0; ii < dim_remained; ii++) {
					if (adj[row][ii] && adj[col][ii]) {
						nei[k] = ii;
						k = k + 1;
					}
				}
				completion_one(M_remained, row, col, nei, k);
				adj[col][row] = 1;
				adj[row][col] = 1;
			}
		}
	}
}

void matrix_completion(gsl_matrix * M)
{
	gsl_set_error_handler_off();
	int dim = M->size1, dim_remained = 0, dim_full = 0, i = 0, j = 0, k = 0;
	int *reduced_idx = malloc(dim * sizeof(int));
	int *full_idx = malloc(dim * sizeof(int));
	for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++) {
			if (gsl_isnan(gsl_matrix_get(M, i, j))) {
				reduced_idx[dim_remained] = i;
				dim_remained = dim_remained + 1;
				break;
			} else if (j == dim - 1) {
				full_idx[dim_full] = i;
				dim_full = dim_full + 1;
			}
		}

	}

	if (dim_remained > 0) {
		if (dim_remained == 2) {
			completion_one(M, reduced_idx[0], reduced_idx[1], full_idx, dim_full);
		} else {
			int empty_length = 0, col = 0, row = 0;
			double temp;
			int **adj = malloc(sizeof(int *) * dim_remained);
			int *empty_idx = NULL;
			int *nei = malloc(sizeof(int) * dim_remained);
			gsl_matrix *M_reducer = gsl_matrix_alloc(dim_remained, dim_remained);
			gsl_matrix *M_remained = gsl_matrix_alloc(dim_remained, dim_remained);
			schur_complement(M, M_remained, M_reducer, reduced_idx, full_idx, dim_remained, dim_full);
			for (i = 0; i < dim_remained; i++) {
				adj[i] = malloc(sizeof(int) * dim_remained);
				for (j = 0; j < dim_remained; j++) {
					if (gsl_isnan(gsl_matrix_get(M_remained, i, j))) {
						adj[i][j] = 0;
						empty_length = empty_length + 1;
					} else {
						adj[i][j] = 1;
					}
				}
			}

			empty_idx = malloc(sizeof(int) * empty_length);
			empty_length = empty_length / 2;
			for (i = 0; i < dim_remained; i++) {
				for (j = i; j < dim_remained; j++) {
					if (gsl_isnan(gsl_matrix_get(M_remained, i, j))) {
						empty_idx[2 * k] = i;
						empty_idx[2 * k + 1] = j;
						k++;
					}
				}
			}

			if (dim_remained == 3) {
				for (i = 0; i < empty_length; i++) {
					row = empty_idx[2 * i + 1];
					col = empty_idx[2 * i];
					// Here, k is the number of shared neibours
					k = 0;
					for (j = 0; j < dim_remained; j++) {
						if (adj[row][j] && adj[col][j]) {
							nei[k] = i;
							k = k + 1;
						}
					}
					completion_one(M_remained, row, col, nei, k);
				}
			} else {
				int *alpha = malloc(sizeof(int) * dim_remained);
				int *alpha_inv = malloc(sizeof(int) * dim_remained);
				int *fill_in = malloc(sizeof(int) * empty_length * 2);
				int num_fill_in;
				MCS_seach_chordal_check(alpha, alpha_inv, adj, dim_remained, &num_fill_in, fill_in);

				if (num_fill_in > 0) {
					gsl_matrix *workspace_chol = gsl_matrix_alloc(dim_remained, dim_remained);
					for (i = 0; i < num_fill_in; i++) {
						row = fill_in[2 * i];
						col = fill_in[2 * i + 1];
						// ############
						guess_fill_in(M_remained, col, row, nei, adj);
						// ############
						adj[col][row] = 1;
						adj[row][col] = 1;
					}
					chordal_completion(M_remained, dim_remained, alpha_inv, nei, adj);
					gsl_matrix_memcpy(workspace_chol, M_remained);
					if (gsl_linalg_cholesky_decomp1(workspace_chol) == 0) {
						optim_newton(M_remained, workspace_chol, empty_length, empty_idx, 100, 1e-5);
					} else {
						make_psd(M_remained, workspace_chol, empty_length, empty_idx, 100);
						gsl_matrix_memcpy(workspace_chol, M_remained);
						if (gsl_linalg_cholesky_decomp1(workspace_chol) == 0) {
							optim_newton(M_remained, workspace_chol, empty_length, empty_idx, 100, 1e-5);
						} else {
							exit(1);
						}
					}
					gsl_matrix_free(workspace_chol);
				} else {
					chordal_completion(M_remained, dim_remained, alpha_inv, nei, adj);
				}

				free(fill_in);
				free(alpha);
				free(alpha_inv);
			}

			for (i = 0; i < empty_length; i++) {
				row = empty_idx[2 * i + 1];
				col = empty_idx[2 * i];
				temp = gsl_matrix_get(M_remained, col, row) + gsl_matrix_get(M_reducer, col, row);
				gsl_matrix_set(M, reduced_idx[row], reduced_idx[col], temp);
				gsl_matrix_set(M, reduced_idx[col], reduced_idx[row], temp);
			}
			for (i = 0; i < dim_remained; i++) {
				free(adj[i]);
			}
			free(nei);
			free(adj);
			free(empty_idx);
			gsl_matrix_free(M_reducer);
			gsl_matrix_free(M_remained);
		}
	}
	free(reduced_idx);
	free(full_idx);

}
