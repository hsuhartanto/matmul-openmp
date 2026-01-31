/* matmul_omp_input.c
 * Compile:
 *   gcc -O3 -fopenmp -march=native -o matmul_omp_input matmul_omp_input.c
 * Run:
 *   ./matmul_omp_input
 * The program will prompt for m p n (A is m x p, B is p x n), number of threads,
 * and algorithm choice (1=naive, 2=blocked). If blocked, it asks block size.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>

static double rand_d(void) { return (double)rand() / (double)RAND_MAX; }

static double *alloc_mat(int r, int c) {
    double *m = malloc(sizeof(double) * (size_t)r * (size_t)c);
    if (!m) { perror("malloc"); exit(EXIT_FAILURE); }
    return m;
}

/* Naive OpenMP matmul: C = A(m x p) * B(p x n) */
static void matmul_naive_omp(const double *A, const double *B, double *C,
                             int m, int p, int n) {
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            C[i*(size_t)n + j] = 0.0;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < m; ++i) {
        for (int k = 0; k < p; ++k) {
            double aik = A[i*(size_t)p + k];
            for (int j = 0; j < n; ++j) {
                C[i*(size_t)n + j] += aik * B[k*(size_t)n + j];
            }
        }
    }
}

/* Blocked (tiling) OpenMP matmul: C = A * B */
static void matmul_blocked_omp(const double *A, const double *B, double *C,
                               int m, int p, int n, int Bsize) {
    #pragma omp parallel for schedule(static)
    for (size_t idx = 0; idx < (size_t)m * n; ++idx) C[idx] = 0.0;

    for (int ii = 0; ii < m; ii += Bsize) {
        int iimax = (ii + Bsize < m) ? ii + Bsize : m;
        for (int kk = 0; kk < p; kk += Bsize) {
            int kkmax = (kk + Bsize < p) ? kk + Bsize : p;
            for (int jj = 0; jj < n; jj += Bsize) {
                int jjmax = (jj + Bsize < n) ? jj + Bsize : n;
                #pragma omp parallel for collapse(2) schedule(static)
                for (int i = ii; i < iimax; ++i) {
                    for (int k = kk; k < kkmax; ++k) {
                        double aik = A[i*(size_t)p + k];
                        for (int j = jj; j < jjmax; ++j) {
                            C[i*(size_t)n + j] += aik * B[k*(size_t)n + j];
                        }
                    }
                }
            }
        }
    }
}

int main(void) {
    int m, p, n, threads;
    int alg_choice = 1;
    int Bsize = 64;

    printf("Enter matrix dimensions (m p n) for A(m x p) * B(p x n): ");
    if (scanf("%d %d %d", &m, &p, &n) != 3) {
        fprintf(stderr, "Invalid input. Expected three integers.\n");
        return EXIT_FAILURE;
    }
    if (m <= 0 || p <= 0 || n <= 0) {
        fprintf(stderr, "Matrix dimensions must be positive.\n");
        return EXIT_FAILURE;
    }

    printf("Enter number of threads to use: ");
    if (scanf("%d", &threads) != 1 || threads <= 0) {
        fprintf(stderr, "Invalid thread count; using 1 thread.\n");
        threads = 1;
    }

    printf("Select algorithm: 1=naive, 2=blocked: ");
    if (scanf("%d", &alg_choice) != 1) alg_choice = 1;
    if (alg_choice == 2) {
        printf("Enter block size (e.g., 32,64): ");
        if (scanf("%d", &Bsize) != 1 || Bsize <= 0) {
            fprintf(stderr, "Invalid block size; using 64.\n");
            Bsize = 64;
        }
    }

    omp_set_num_threads(threads);
    printf("Computing A(%d x %d) * B(%d x %d) with %d thread(s), alg=%d\n",
           m, p, p, n, threads, alg_choice);

    srand((unsigned)time(NULL));
    double *A = alloc_mat(m, p);
    double *B = alloc_mat(p, n);
    double *C = alloc_mat(m, n);

    for (size_t i = 0; i < (size_t)m * p; ++i) A[i] = rand_d();
    for (size_t i = 0; i < (size_t)p * n; ++i) B[i] = rand_d();

    double t0 = omp_get_wtime();
    if (alg_choice == 1) {
        matmul_naive_omp(A, B, C, m, p, n);
    } else {
        matmul_blocked_omp(A, B, C, m, p, n, Bsize);
    }
    double t1 = omp_get_wtime();

    double checksum = 0.0;
    for (size_t i = 0; i < (size_t)m * n; ++i) checksum += C[i];

    printf("Time elapsed: %.6f s\n", t1 - t0);
    printf("Checksum (sum of C elements): %.12e\n", checksum);

    free(A); free(B); free(C);
    return 0;
}

