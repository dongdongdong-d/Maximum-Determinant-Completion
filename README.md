This repository is a program doing matrix completion. The completed matrix should be a positive definite matrix with maximized determinants. The input is a matrix with specified diagonal entries and principal minors, composed of specified entries, are positive.

The general steps of the algorithm are the following:
1. A Schur complement is applied to the incomplete matrix to reduce the size of the matrix by eliminating those full columns and rows.
2. Build an undirected graph based on the pattern of the missing values. The graph is built in the way:
    2.1 If the matrix is of dimension N*N, N nodes are created.
    2.2 If the entry (i,j) is specified, node i and node j are connected by an edge. Otherwise, they are not connected.
3. Determine if the graph is chordal by using the maximum cardinality search[1], which will also provide a chordal completion order (not a minimum completion). A chordal completion order makes a matrix chordal.
4. If the graph is chordal, use the method shown in [2] to complete the matrix.
5. If the graph is not chordal, do a positive definite completion using the following heuristic way first:
    5.1 Use the chordal completion order obtained in step 3 to make the matrix chordal. If we want to fill in entry (i,j), we consider the largest principal submatrices with an only entry (i,j) unspecified. The determinant of those principal submatrices should be larger than zero. Thus, we can have a set of inequalities to bound the entry (i,j).
    5.2 The largest principal submatrices are found on the graph. The task is equivalent to finding the largest cliques consisting of shared neighbor nodes of node i and node j. The largest cliques are found by Bron-Kerbosch algorithm with pivoting, where the pivot is the most connected node.[3]
    5.3 Then, apply step 4 to get a matrix completion. If the matrix completion is positive definite, jump to step 7. Otherwise, proceed to step 6.
6. The matrix is now completed but not positive definite from step 5. We apply the barrier method shown in [4] to get a definite positive completion. The general idea is:
    6.1 Add xI to the matrix, where I is an identity matrix, and x makes the matrix positive definite.
    6.2 Maximize the determinant of the matrix, with those unspecified entries as function arguments.
    6.3 Reduce the matrix by yI, where y keeps the matrix positive definite.
    6.4 Repeat 6.2 and 6.3 until the cumulative reduction cancels out x in 6.1. Then a positive definite completion is obtained.
7. The matrix is now completed to be a positive definite matrix by either 5 or 6. We can now maximize the determinant of the matrix by the Newton Raphson method. The gradient and Hessian can be computed analytically. Note step 7 and step 6.2 use the same function.

Some Notes:

1. maxdet_completion.c is the source code of the library. 
    ---The major function in this library is matrix_completion.
2. read_and_complete.c reads partial matrices in "OUT.new" and completes them. 
3. The program depends on gsl. We can compile it with intel mkl.


Reference:
[1]Tarjan, Robert E., and Mihalis Yannakakis. "Simple linear-time algorithms to test chordality of graphs, test acyclicity of hypergraphs, and selectively reduce acyclic hypergraphs." SIAM Journal on computing 13.3 (1984): 566-579.
[2]Grone, Robert, et al. "Positive definite completions of partial Hermitian matrices." Linear algebra and its applications 58 (1984): 109-124.
[3]Tomita, Etsuji, Akira Tanaka, and Haruhisa Takahashi. "The worst-case time complexity for generating all maximal cliques and computational experiments." Theoretical computer science 363.1 (2006): 28-42.
[4]Glunt, W., et al. "Positive definite completions and determinant maximization." Linear algebra and its applications 288 (1999): 1-10.
