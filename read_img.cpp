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

#define BUFFERSIZE 2797568
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define FILTER_HEIGHT 3
#define FILTER_WIDTH 3

#define ODD(n) (n % 2 == 0 ? false : true)
#define EVEN(n) (n % 2 == 0 ? true : false)


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

int printVector(vector<double> vec) {
    for (int i = 0; i < vec.size(); i++)
    { // iterate over the rows of the 2D vector
        cout << vec[i] << " "; // output the element at row i, column j
    }
        cout << endl; // output a newline character after each row
    return 1;
}

double sobelTemp(vector<vector<double>> matrix) {
    if (matrix.size() != 3 || matrix[0].size() != 3) {
        printf("Error matrix has to be 3X3.");
        MPI_Finalize(); exit(0); return 0;
    }
    double sx, sy, s;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            sx += Gx[i][j] * matrix[i][j];
            sy += Gy[i][j] * matrix[i][j];
        }
    }
    s = sqrt(pow(sx, 2) + pow(sy, 2));
    return s;
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
    vector<vector<double>> mag(rows, vector<double> (cols - 2, 0));
    printf("heigh: %lu, width: %lu\n", mag.size(), mag[0].size());

    for (int i = 0; i < rows - 2; i++) {
        for (int j = 0; j < cols - 2; j++) {
            double s1 = sobelHelper(Gx, mat, i, j);
            double s2 = sobelHelper(Gy, mat, i, j);
            mag[i + 1][j] = sqrt(pow(s1, 2) + pow(s2, 2)); // j+1?
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
        MPI_Finalize(); exit(0); return vector<vector<double>>();
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
    clock_t startTime;
    int p, P, tag, n;
    tag = 42;
    MPI_Status status;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
        
    int matrix_size[2];
    vector<double> flatten_list;

    // Master processor reads image
    if (p==0) {
        clock_t startTime = clock(); // Master processor keeps the time

        vector<vector<double>> mat = getMatrix("images/mat_files/index.csv");

        int N = (int) mat.size(); 
        int M = (int) mat[0].size();
        matrix_size[0] = N;
        matrix_size[1] = M;
        cout << "N: " << N <<", M: " << M << endl;
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

    vector<int> sendcnts(P);
    vector<int> displs(P);
    displs[0] = 0; 
    int W = M;
    int H = (int) N/P;
    int res = N - P * H;


    for (int i = 0; i < P - 1; i++) {
        sendcnts[i] = W * (i < res ? H + 1 : H);
        displs[i+1] = displs[i] + sendcnts[i];

    }
    sendcnts[P - 1] = W * H;

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
                local_flattened_list.data(), sendcnts[p], MPI_DOUBLE, 0, MPI_COMM_WORLD);


    flatten_list.clear();

    H = H + (p < res ? 1 : 0);
    cout << "p: " << p << ", H: " << H << ", W: " << W << " flattened vec: " <<  local_flattened_list.size() << endl;
    // cout << "p" << p << ": " << local_flattened_list[0] << ", " << local_flattened_list[1] << endl;
    vector<vector<double>> local_mat = reshape(local_flattened_list, H, W);
    if (p==0){
        // printMatrix(local_vector);
        // cout << "local_vector.size: " << local_vector.size() << ", ";
        // cout << "local_vector[0].size: " << local_vector[0].size() << endl; 
    }

    int W_mag = W - 2;
    int H_mag = H + (p==0 || p==P-1 ? 1 : 0);
    vector<double> mag_top(W_mag); 
    vector<double> mag_down(W_mag); 
    vector<vector<double>> mat_temp(3 , vector<double> (3, 0)); // Create a 3x3 matrix with zeros  
    int j, k;
    double s;
    // Sending downwards and receiving from upwards
    for (int i = 0; i < W-2; i++) {         
        // First iteration, send an receive 3 elements                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
        
        // initialize the original temp_matrix
        if (p!=0) {
            for (j = 1; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    mat_temp[j][k] = local_mat[j-1][k+i];
                }
            }
        }
        // } else if (i!=0 && p!=0) {
        //     for (j = 0; j < 3; j ++) {
        //         for (k = 0; k < 2; k++) {
        //             mat_temp[j][k] = mat_temp[j][k+1];
        //         }
        //     }
        // }

        if (i==0 && EVEN(p)) {
            if (p != P-1)
                MPI_Send(local_mat[H-1].data(), 3, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD); // sending down
            
            if (p != 0)
                MPI_Recv(mat_temp[0].data(), 3, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status); // receiving from upper

        } else if (i==0 && ODD(p)) {
            if (p != 0)
                MPI_Recv(mat_temp[0].data(), 3, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status); // receive from upper

            if (p != P-1)
                MPI_Send(local_mat[H-1].data(), 3, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD); // sending down
        
        } else if (i!=0 && EVEN(p)) {
            // Shift the 1st row one step to the left
            for (k = 0; k < 2; k++)
                mat_temp[0][k] = mat_temp[0][k+1];
            
            if (p != P-1)
                MPI_Send(local_mat[H-1].data() + 2 + i, 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD); // sending down

            if (p != 0)
                MPI_Recv(mat_temp[0].data() + 2, 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status); // receiving from upper
            
        } else if (p!=0 && ODD(p)) {
            // Shift the 1st row one step to the left
            for (k = 0; k < 2; k++)
                mat_temp[0][k] = mat_temp[0][k+1];
                
            if (p != 0)
                MPI_Recv(mat_temp[0].data() + 2, 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status); // receive from upper
 
            if (p != P-1)
                MPI_Send(local_mat[H-1].data() + 2 + i, 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD); // sending down
        }

        mag_top[i] = sobelTemp(mat_temp);

    }


    // Sending upwards and receiving from downwards
    for (int i = 0; i < W-2; i++) {         
        // First iteration, send an receive 3 elements                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
        
        // initialize the original temp_matrix
        if (p!=P-1) {
            for (j = 0; j < 2; j++) {
                for (k = 0; k < 3; k++) {
                    mat_temp[j][k] = local_mat[H-2+j][k+i];
                }
            }
        }

        if (i==0 && EVEN(p)) {
            if (p != 0)
                MPI_Send(local_mat[0].data(), 3, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD); // sending up
            
            if (p != P-1)
                MPI_Recv(mat_temp[2].data(), 3, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status); // receiving from down

        } else if (i==0 && ODD(p)) {
            if (p != P-1)
                MPI_Recv(mat_temp[2].data(), 3, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status); // receive from down

            if (p != 0)
                MPI_Send(local_mat[0].data(), 3, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD); // sending up
        
        } else if (i!=0 && EVEN(p)) {
            // Shift the last row one step to the left
            for (k = 0; k < 2; k++)
                mat_temp[2][k] = mat_temp[2][k+1];
            
            if (p != 0)
                MPI_Send(local_mat[0].data() + 2 + i, 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD); // sending up

            if (p != P-1)
                MPI_Recv(mat_temp[2].data() + 2, 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status); // receiving from down
            
        } else if (p!=0 && ODD(p)) {
            // Shift the last row one step to the left
            for (k = 0; k < 2; k++)
                mat_temp[2][k] = mat_temp[2][k+1];
                
            if (p != P-1)
                MPI_Recv(mat_temp[2].data() + 2, 1, MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status); // receive from down
 
            if (p != 0)
                MPI_Send(local_mat[0].data() + 2 + i, 1, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD); // sending up
        }
        
        mag_down[i] = sobelTemp(mat_temp);
    }

    // Compute the mag for the local matrix
    vector<vector<double>> mag = sobel(local_mat);

    // Add the top and the bottom mag to the total mag
    if (p==0) {
        mag[H-1] = mag_down; // The top processor/stripe do not have a mag_top
    } else if (p==P-1) {
        mag[0] = mag_top; // The top processor/stripe do not have a mag_down
    } else {
        mag[0] = mag_top;
        mag[H-1] = mag_down;        
    }
    
    if (p==0) {
        for (int i = 0; i < mag.size(); i++) {
            cout << "p" << p << ", i: " << i << ", W: " << mag[i].size() << endl; 
        }        
    }
    // for (int i = 0; i < mag.size(); i++) {
    //     cout << "p" << p << ", i: " << i << ", W: " << mag[i].size() << endl;
    // }


    /* output for graphical representation */
    /* Instead of using gather (which may lead to excessive memory requirements
    on the master process) each process will write its own data portion. This
    introduces a sequentialization: the hard disk can only write (efficiently)
    sequentially. Therefore, we use the following strategy:
    1. The master process writes its portion. (file creation)
    2. The master sends a signal to process 1 to start writing.
    3. Process p waites for the signal from process p-1 to arrive.
    4. Process p writes its portion to disk. (append to file)
    5. process p sends the signal to process p+1 (if it exists).
    */
    FILE *fp;
    if (p==0){ // Master process
        fp = fopen("test.csv", "w");
		for (int i = 0; i < mag.size(); i++) {
            for (int j = 0; j < mag[i].size()-1; j++) {
		        fprintf(fp, "%f,", mag[i][j]);
            }
            fprintf(fp, "%f\n", mag[i][j+1]);
        }
        fclose(fp);		
        MPI_Send("hi", 2, MPI_CHAR, 1, tag, MPI_COMM_WORLD);
    } else {
        
        char message[2]; // Nonesense message
        MPI_Recv(message, 2, MPI_CHAR, p-1, tag, MPI_COMM_WORLD, &status);
        fp = fopen("test.csv", "a");
		for (int i = 0; i < mag.size(); i++) {
            for (int j = 0; j < mag[i].size()-1; j++) {
		        fprintf(fp, "%f,", mag[i][j]);
            }
            fprintf(fp, "%f\n", mag[i][j+1]);
            
        }
        // fprintf(fp, "\n");
        if (p!=P-1) { // If it's not the last processor, keep sending the signals forward
            MPI_Send(message, 2, MPI_CHAR, p+1, tag, MPI_COMM_WORLD);
        }
        fclose(fp);
    }
    

    if (p==0) {
        clock_t secElapsed = (startTime - clock()) / 1000000;
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
    // printf("Process %d finished\n", p);
    // /* That's it */
    MPI_Finalize();
    exit(0);
}