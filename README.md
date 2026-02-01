the codes is matrix multipliction in parallel with openmp,
/* matmul_omp_input.c
 * Compile:
 *   gcc -O3 -fopenmp -march=native -o matmul_omp_input matmul_omp_input.c
 * Run:
 *   ./matmul_omp_input
 * The program will prompt for m p n (A is m x p, B is p x n), number of threads,
 * and algorithm choice (1=naive, 2=blocked). If blocked, it asks block size.
 */

examples of some testings:

root@user11-gpu02-cuda10:/var/nfs/hs/matmul-openmp# ls
README.md  matmul-openmp.c  matmul_omp_input
root@user11-gpu02-cuda10:/var/nfs/hs/matmul-openmp# ./matmul_omp_input
Enter matrix dimensions (m p n) for A(m x p) * B(p x n): 2048 2048 2048
Enter number of threads to use: 2
Select algorithm: 1=naive, 2=blocked: 1
Computing A(2048 x 2048) * B(2048 x 2048) with 2 thread(s), alg=1
Time elapsed: 3.448814 s
Checksum (sum of C elements): 2.147198825109e+09

root@user11-gpu02-cuda10:/var/nfs/hs/matmul-openmp# ./matmul_omp_input
Enter matrix dimensions (m p n) for A(m x p) * B(p x n): 2048 2048 2048
Enter number of threads to use: 4
Select algorithm: 1=naive, 2=blocked: 1
Computing A(2048 x 2048) * B(2048 x 2048) with 4 thread(s), alg=1
Time elapsed: 2.728874 s
Checksum (sum of C elements): 2.148178858097e+09

root@user11-gpu02-cuda10:/var/nfs/hs/matmul-openmp# ./matmul_omp_input
Enter matrix dimensions (m p n) for A(m x p) * B(p x n): 2048 2048 2048
Enter number of threads to use: 8
Select algorithm: 1=naive, 2=blocked: 1
Computing A(2048 x 2048) * B(2048 x 2048) with 8 thread(s), alg=1
Time elapsed: 2.462738 s
Checksum (sum of C elements): 2.147265618046e+09

