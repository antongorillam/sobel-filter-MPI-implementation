#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <string.h>
#include <vector>
#include <cmath>
#include <time.h>
#include <chrono>

/* Use MPI */
#include <mpi.h>


using namespace std;

#define BUFFERSIZE 16384
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define FILTER_HEIGHT 3
#define FILTER_WIDTH 3

const vector<vector<double>> Gx = {
        {-1, 0, 1},
        {-2, 0, 2},
        {-1, 0, 1}
    };
    
const vector<vector<double>> Gy = {
        {-1, -2, -1},
        { 0,  0,  0},
        { 1,  2,  1}
    };


/**
 * Reads a CSV file from the given path and returns its contents as a 2D vector of integers.
 *
 * @param path The path to the CSV file.
 * @return A 2D vector of integers representing the contents of the CSV file.
 */
vector<vector<double>> getMatrix(string path) {
    FILE* fp;
    fp = fopen(path.c_str(), "r");
    if (fp == NULL) {
        printf("Error opening file %s\n", path.c_str());
        exit(1);
    }
    vector<vector<double>> data;
    char buffer[BUFFERSIZE];
    while (fgets(buffer, BUFFERSIZE, fp) != NULL) {
        char* field = strtok(buffer, ",");
        vector<double> row;
        while (field != NULL) {
            double intField = strtod(field, NULL);
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
int printMatrix(vector<vector<double>> mat) {
    for (int i = 0; i < mat.size(); i++)
    { // iterate over the rows of the 2D vector
        for (int j = 0; j < mat[i].size(); j++) { // iterate over the columns of each row 
            cout << mat[i][j] << " "; // output the element at row i, column j
        }
        cout << endl; // output a newline character after each row
    }
    return 1;
}

double sobelHelper(vector<vector<double>> filter, vector<vector<double>> matrix, int start_row, int start_col) {
    double sum = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            sum += filter[i][j] * matrix[i + start_row][j + start_col];
        }
    }
    return sum;
    
}

vector<vector<double>> sobel(vector<vector<double>> mat) {
  
    int rows = mat.size(); int cols = mat[0].size(); 
    vector<vector<double>> mag(rows, vector<double> (cols, 0));
    printf("heigh: %lu, width: %lu\n", mag.size(), mag[0].size());

    for (int i = 0; i < rows - 2; i++) {
        for (int j = 0; j < cols - 2; j++) {
            double s1 = sobelHelper(Gx, mat, i, j);
            double s2 = sobelHelper(Gy, mat, i, j);
            mag[i + 1][j + 1] = sqrt(pow(s1, 2) + pow(s2, 2));
        }   
    }
    return mag;
}

/**
 * Reshapes a 1D vector into a 2D vector with the specified number of rows and columns.
 *
 * @param flat_vector The 1D vector to reshape.
 * @param rows The number of rows in the reshaped 2D vector.
 * @param cols The number of columns in the reshaped 2D vector.
 * @return The reshaped 2D vector.
 * @throws std::runtime_error if the number of elements in the 1D vector is not equal to rows * cols.
 */
vector<vector<double>> reshape(vector<double> flat_vector, int rows, int cols) {

    if (rows * cols != flat_vector.size()) {
        printf("Invalid reshape!");
        MPI_Finalize(); exit(0); return 0;
    }
    
    vector<vector<double>> vector_2d(rows, vector<double>(cols));
    int k = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            vector_2d[i][j] = flat_vector[k];
            k++;
        }
    }
    return vector_2d;
}

int main(int argc, char* argv[])
{
    /* local variable */
    /* Initialize MPI */
    clock_t startTime;
    int p, totalProcs, tag, n;
    tag = 42;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &totalProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
    
    int matrix_size[2];
    vector<double> flatten_list;
    // Master processor reads v
    if (p==0) {
        clock_t startTime = clock(); // Master processor keeps the time

        vector<vector<double>> mat = getMatrix("images/mat_files/100075.jpg.csv");

        int N = (int) mat.size(); int M = (int) mat[0].size();
        matrix_size[0] = N;
        matrix_size[1] = M;

        // Flatten the vector to a array object.
        flatten_list.resize(M * N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                flatten_list[i * M + j] = mat[i][j];
            }
        }
    }
    MPI_Bcast(&matrix_size, 2, MPI_INT, 0, MPI_COMM_WORLD);
    int N = matrix_size[0]; int M = matrix_size[1];

    vector<int> sendcnts(totalProcs);
    vector<int> displs(totalProcs);
    displs[0] = 0; 
    int W = M;
    int H = (int) N/totalProcs;
    int res = N - totalProcs * H;


    for (int i = 0; i < totalProcs - 1; i++) {
        sendcnts[i] = W * (i < res ? H + 1 : H);
        displs[i+1] = displs[i] + sendcnts[i];

    }
    sendcnts[totalProcs - 1] = W * H;

    if (p==0) {
        // cout << "p:" << p << " hej" << endl;
        for (int i = 0; i < sendcnts.size(); i++) {
            cout << sendcnts[i] << ", ";
        }
        cout << endl;
        for (int i = 0; i < sendcnts.size(); i++) {
            cout << displs[i] << ", ";
        }            
        cout << endl;
    }
    vector<double> local_flattened_list(sendcnts[p]);


    // Do the scatter stuff now
    MPI_Scatterv(flatten_list.data(), sendcnts.data(), displs.data(), MPI_DOUBLE, 
                local_flattened_list.data(), sendcnts[p], MPI_DOUBLE, 0, MPI_COMM_WORLD );

    free(flatten_list);
    H = H + (p < res ? 1 : 0);
    cout << "p: " << p << ", H: " << H << ", W: " << W << " flattened vec: " <<  local_flattened_list.size() << endl;
    // cout << "p" << p << ": " << local_flattened_list[0] << ", " << local_flattened_list[1] << endl;
    vector<vector<double>> local_vector = reshape(local_flattened_list, H, W);
    if (p==0){
        // printMatrix(local_vector);
        cout << "local_vector.size: " << local_vector.size() << ", ";
        cout << "local_vector[0].size: " << local_vector[0].size() << endl; 
    }
    



    // cout << "p: " << p << ", p_i: " << p_i << ", p_j: " << p_j <<  endl;

    // cout << "p: " << p << ", p_i: " << p_i << ", p_j: " << p_j <<  endl;


    // cout << "H_res: " << H_res << endl;    

    // double flatten_list[(int) matrix_size[0] * matrix_size[1]];
    

    // Send the size of the matrix to all the v

    // MPI_Scatter(v.data(), local_n, MPI_INT, local_v.data(), local_n, MPI_INT, 0, MPI_COMM_WORLD);
    


    if (p==0) {
        clock_t secElapsed = (clock() - startTime) / 1000000;
        // ofstream outfile("matrix.csv");
        // for (int i = 0; i < mag.size(); i++) {
        //     for (int j = 0; j < mag[i].size(); j++) {
        //         outfile << mag[i][j];
        //         if (j != mag[i].size() - 1) {
        //             outfile << ",";
        //         }
        //     }
        //     outfile << std::endl;
        // }
        // outfile.close();
        cout << "Time elapsed: " << secElapsed << " seconds" << endl;
    }
    



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