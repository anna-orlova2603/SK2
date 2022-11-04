#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double generate(x, y)
{
    return ((double)rand()/(double)(RAND_MAX)) * (y - x) + x;
}

int main(int argc, char *argv[])
{ 
    int size, rank;
    long N = 0, i = 0;
    double integral = 1.0/364.0, epsilon = strtod(argv[1], NULL);
    double a1 = 0.0, a2 = 1.0, b1 = 0.0, b2 = 1.0, c1 = 0.0, c2 = 1.0;
    int par = 1, root = 0;
    double x, y, z;
    double sum = 0.0, f_sum, func, int_appr = 0.0, error = 0.0, sum_2 = 0;
    double time_start, time_end, ttime, timee;
    int seed, M = 0; // seeds = {1, 13, 20, 614, 1300, 2022, 809, 5575, 5555}

    sscanf(argv[1], "%lf", &epsilon);
    sscanf(argv[2], "%d", &seed);
    sscanf(argv[3], "%d", &M);

    double buff[M];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand(seed + rank); 

    time_start = MPI_Wtime();
 
    while (1)
    {
        N++;
        
        for(int j = 0; j < M; j++)
        { 
            x = generate(a1, a2);
            y = generate(b1, b2);
            z = generate(c1, c2);

            if (a1 <= x <= a2 && b1 <= y <= b2 && c1 <= z <= c2 && x * y >= z && y <= x)  
            {
                func = x*pow(y, 2)*pow(z, 3);
                i++;
            }
            else 
            {
                func = 0.0; 
            }

            buff[j] = func;
        }

        MPI_Reduce(&buff, &f_sum, M, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == root)
        {
            sum = sum + f_sum/((double)size);
            int_appr = par * sum/ ((double)N); 
            error = fabs(integral - int_appr);
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&error, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

        if (error < epsilon)
            {
                break;
            }
    } 

    time_end = MPI_Wtime();
    timee = time_end - time_start;
    MPI_Reduce(&timee, &ttime, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);

    MPI_Finalize();

    if (rank == root)
    {
        printf("Integral approximate = %f\n", int_appr);
        printf("Integral numerical = %f\n", integral);
        printf("Error:  %e   %.10f\n", error, error);
        printf("Epsilon:  %.10f\n", epsilon);
        printf("Seed:  %d\n", seed);
        printf("M size (buff):  %d\n", M);
        printf("Process number (all points) N*size:  %ld\n", M*N*size);
        printf("Max time:  %.8f\n", ttime);
        printf("Size: %d\n", size);
    }
    return 0;
}