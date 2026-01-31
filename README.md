the codes is matrix multipliction in parallel with openmp,
/* matmul_omp_input.c
 * Compile:
 *   gcc -O3 -fopenmp -march=native -o matmul_omp_input matmul_omp_input.c
 * Run:
 *   ./matmul_omp_input
 * The program will prompt for m p n (A is m x p, B is p x n), number of threads,
 * and algorithm choice (1=naive, 2=blocked). If blocked, it asks block size.
 */
