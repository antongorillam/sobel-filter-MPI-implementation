#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <string.h>
#include <vector>
/* Use MPI */
#include <mpi.h>


using namespace std;

#define BUFFERSIZE 16384
#define MIN(a, b) ((a) < (b) ? (a) : (b))

/**
 * Reads a CSV file from the given path and returns its contents as a 2D vector of integers.
 *
 * @param path The path to the CSV file.
 * @return A 2D vector of integers representing the contents of the CSV file.
 */
vector<vector<int>> getMatrix(string path) {
    FILE* fp;
    fp = fopen(path.c_str(), "r");
    if (fp == NULL) {
        printf("Error opening file %s\n", path.c_str());
        exit(1);
    }
    vector<vector<int>> data;
    char buffer[BUFFERSIZE];
    while (fgets(buffer, BUFFERSIZE, fp) != NULL) {
        char* field = strtok(buffer, ",");
        vector<int> row;
        while (field != NULL) {
            int intField = strtod(field, NULL);
            row.push_back(intField);
            field = strtok(NULL, ",");
        }
        data.push_back(row);
    }
    fclose(fp);
    return data;
}

/**
 * Prints the contents of a 2D vector of integers to the console.
 *
 * @param mat The 2D vector of integers to print.
 * @return 1 if the function executed successfully.
 */
int printMatrix(vector<vector<int>> mat) {
    for (int i = 0; i < mat.size(); i++)
    { // iterate over the rows of the 2D vector
        for (int j = 0; j < mat[i].size(); j++) // iterate over the columns of each row 
        {
            std::cout << mat[i][j] << " "; // output the element at row i, column j
        }
        std::cout << std::endl; // output a newline character after each row
    }
    return 1;
}

int main(int argc, char* argv[])
{
    /* local variable */
    /* Initialize MPI */

    int p, P, tag, n;
    tag = 42;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);

    printf("p: %d\n", p);

    vector<vector<int>> mat = getMatrix("images/mat_files/8049.jpg.csv");
    printf("heigh: %lu, width: %lu\n", mat.size(), mat[0].size());
    // printMatrix(mat);

    MPI_Finalize(); exit(0); return 0;
    // if (N < P) {
    //     fprintf(stdout, "Too few discretization points...\n");
    //     exit(1);
    // }

    // unew = (double *) malloc((int) L*sizeof(double));
    // u = (double *) calloc(u_size, sizeof(double));

    // /* Jacobi iteration */
    // for (int step = 0; step < K; step++) {
    //     if (p==0 && step % 1000 == 0) {
    //         printf("At step: %d/%d\n", step, K);
    //     }
    //     /* RB communication of overlap */

    //     if (p == 0) {
    //         MPI_Send(&u[u_size-2], 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD);
    //         MPI_Recv(&u[u_size-1], 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status);
    //     } else if (p == P-1 && p % 2 == 0) { // last processor and even
    //         MPI_Send(&u[1], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);
    //         MPI_Recv(&u[0], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status);
    //     } else if (p == P-1 && p % 2 != 0) { // last processor and odd
    //         MPI_Recv(&u[0], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status);
    //         MPI_Send(&u[1], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);                    
    //     } else if (p % 2 == 0) { // Even: red
    //         MPI_Send(&u[u_size-2], 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD);
    //         MPI_Recv(&u[u_size-1], 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status);
    //         MPI_Send(&u[1], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);
    //         MPI_Recv(&u[0], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status);
    //     } else { // Odd: black
    //         MPI_Recv(&u[0], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status);
    //         MPI_Send(&u[1], 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);
    //         MPI_Recv(&u[u_size-1], 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status);
    //         MPI_Send(&u[u_size-2], 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD);
    //     }

    //     // /* local iteration step */
    //     for (int i = 0; i < L; i++) {
    //         n = p*L+MIN(p,R)+i;
    //         unew[i] = (u[i]+u[i+2]-h*h*f((double) n * h))/(2.0-h*h*r((double) n * h));

    //     }
    //     for (int i = 0; i < L; i++) {
    //         u[i+1] = unew[i]; 
    //     }
    // }

    // FILE *fp;
    // if (p==0){ // Master process
    //     fp = fopen("hm2_test_res.csv", "w");
    // 	for (int i = 0; i < L; i++) {
    // 	    fprintf(fp, "%f, ", unew[i]);
    //     }
    //     fclose(fp);		
    //     MPI_Send("hi", 2, MPI_CHAR, 1, tag, MPI_COMM_WORLD);
    // } else {

    //     char message[2]; // Nonesense message
    //     MPI_Recv(message, 2, MPI_CHAR, p-1, tag, MPI_COMM_WORLD, &status);
    //     fp = fopen("hm2_test_res.csv", "a");
    //     for (int i = 0; i < L; i++) {
    // 	    fprintf(fp, "%f, ", unew[i]);
    //     }
    //     // fprintf(fp, "\n");
    //     if (p!=P-1) { // If it's not the last processor, keep sending the signals forward
    //         MPI_Send(message, 2, MPI_CHAR, p+1, tag, MPI_COMM_WORLD);
    //     }
    //     fclose(fp);
    // }


    // free(u);
    // free(unew);

    // printf("Process %d finished\n", p);
    // /* That's it */
    MPI_Finalize();
    exit(0);
}