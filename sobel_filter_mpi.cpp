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
#include <filesystem>
#include <libgen.h>
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

/**
 * @brief Print the elements of a vector of doubles.
 *
 * This function iterates through a given vector of doubles and prints each element followed by a space.
 * After printing all elements, it outputs a newline character.
 *
 * @param vec The vector of doubles to be printed.
 * @return 1 upon successful execution.
 */
int printVector(vector<double> vec) {
    for (int i = 0; i < vec.size(); i++) {
        // Output the element at index i
        cout << vec[i] << " ";
    }

    // Output a newline character after printing all elements
    cout << endl;

    return 1;
}

/**
 * Calculates the Sobel operator on a 3x3 matrix.
 * 
 * @param matrix The 3x3 input matrix.
 * @return The Sobel operator value.
 */
double sobel3x3(vector<vector<double>> matrix) {
    if (matrix.size() != 3 || matrix[0].size() != 3) {
        // Print an error message and exit the program if the matrix is not 3x3.
        printf("Error matrix has to be 3X3.");
        MPI_Finalize(); exit(0); return 0;
    }
    double sx, sy, s;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            // Compute the Sobel operator values for each element in the matrix.
            sx += Gx[i][j] * matrix[i][j];
            sy += Gy[i][j] * matrix[i][j];
        }
    }
    // Compute the magnitude of the Sobel operator.
    s = sqrt(pow(sx, 2) + pow(sy, 2));
    return s;
}

vector<double> flattenMat(vector<vector<double>> mat) {
    int H = mat.size(); int W = mat[0].size();
    vector<double> vec;
    for (const auto &row : mat) {
        for (const auto &elem : row) {
            vec.push_back(elem);
        }
    }
    return vec;
}

/**
 * Helper funtion for sobelStar() calculates the magnitude of the gradient of an image at a given pixel using the vertical Sobel operator.
 * 
 * @param image The 1D array of doubles representing the image.
 * @param x The x coordinate of the pixel.
 * @param y The y coordinate of the pixel.
 * @param width The width of the image.
 * @param Gy The 3x3 kernel to use for vertical gradient approximation.
 * @return The magnitude of the gradient of the image at the given pixel.
 */
double sobelStarHelperX(vector<double> image, int width, int x, int y) {

    double mag = 0.0; // this your magnitude

    int xn, yn, index;

    xn = x - 1;
    yn = y - 1;
    index = xn + yn * width;
    mag += image[index] * Gx[0][0];

    xn = x;
    index = xn + yn * width;
    mag += image[index] * Gx[1][0];

    xn = x + 1;
    index = xn + yn * width;
    mag += image[index] * Gx[2][0];

    yn = y + 1;
    index = xn + yn * width;
    mag += image[index] * Gx[2][2];

    xn = x;
    index = xn + yn * width;
    mag += image[index] * Gx[1][2];

    xn = x - 1;
    index = xn + yn * width;
    mag += image[index] * Gx[0][2];
    
    return mag;
}

/**
 * Helper funtion for sobelStar() calculates the magnitude of the gradient of an image at a given pixel using the horizontal Sobel operator.
 * 
 * @param image The 1D array of doubles representing the image.
 * @param x The x coordinate of the pixel.
 * @param y The y coordinate of the pixel.
 * @param width The width of the image.
 * @param Gy The 3x3 kernel to use for horizontal gradient approximation.
 * @return The magnitude of the gradient of the image at the given pixel.
 */
double sobelStarHelperY(vector<double> image, int width, int x, int y) {

    double mag = 0.0; // this is your magnitude

    int xn, yn, index;

    xn = x - 1;
    yn = y - 1;
    index = xn + yn * width;
    mag += image[index] * Gy[0][0];

    yn = y;
    index = xn + yn * width;
    mag += image[index] * Gy[0][1];
    
    yn = y + 1;
    index = xn + yn * width;
    mag += image[index] * Gy[0][2];

    xn = x + 1;
    index = xn + yn * width;
    mag += image[index] * Gy[2][2];

    yn = y;
    index = xn + yn * width;
    mag += image[index] * Gy[2][1];
    
    yn = y - 1;
    index = xn + yn * width;
    mag += image[index] * Gy[2][0];
    
    return mag;
}

/**
 * Applies the best sequential Sobel edge detection algorithm to the given 2D matrix of doubles.
 * 
 * @param mat The 2D matrix of doubles to apply Sobel edge detection on.
 * @return A new 2D matrix of doubles containing the magnitude of the edges detected.
 */
vector<vector<double>> sobelStar(vector<vector<double>> mat) {
    vector<double> vec = flattenMat(mat); // Flattens the matrix
    int H = mat.size();
    int W = mat[0].size();
    vector<vector<double>> mag(H, vector<double> (W - 2, 0));

    double s1, s2;
    for (int h = 1; h < mat.size() - 1; h++) {
        for (int w = 1; w < W - 1; w++) {
            s1 = sobelStarHelperX(vec, W, w, h);
            s2 = sobelStarHelperY(vec, W, w, h);
            mag[h][w] = sqrt(pow(s1, 2) + pow(s2, 2));
        }
    }
    return mag;
}

/**
 * Applies the Sobel operator to a 3x3 neighborhood of the input image matrix.
 * 
 * @param filter: A 2D vector of doubles representing the Sobel operator filter.
 * @param matrix: A 2D vector of doubles representing the input image.
 * @param start_row: The row index of the top-left corner of the 3x3 neighborhood.
 * @param start_col: The column index of the top-left corner of the 3x3 neighborhood.
 * @return The dot product of the Sobel operator filter and the input image neighborhood.
 */
double sobelHelper(vector<vector<double>> filter, vector<vector<double>> matrix, int start_row, int start_col) {
    // Initialize the sum to 0.
    double sum = 0;
    // Compute the dot product between the Sobel operator filter and the input image neighborhood.
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            sum += filter[i][j] * matrix[i + start_row][j + start_col];
        }
    }
    // Return the dot product.
    return sum;   
}

/**
 * Applies the Sobel operator to the input image matrix.
 * 
 * @param mat: A 2D vector of doubles representing the input image.
 * @return A 2D vector of doubles representing the magnitude of the gradient at each pixel.
 */
vector<vector<double>> sobel(vector<vector<double>> mat) {
    // Get the number of rows and columns in the input image.
    int rows = mat.size(); int cols = mat[0].size(); 

    // Create a new matrix to store the magnitude of the gradient at each pixel.
    vector<vector<double>> mag(rows, vector<double> (cols - 2, 0));

    // For each pixel in the input image, apply the Sobel operator to compute the gradient magnitude.
    for (int i = 0; i < rows - 2; i++) {
        for (int j = 0; j < cols - 2; j++) {
            // Compute the gradient in the x and y directions using the Sobel operator.
            double s1 = sobelHelper(Gx, mat, i, j);
            double s2 = sobelHelper(Gy, mat, i, j);
            // Compute the magnitude of the gradient at the current pixel.
            mag[i + 1][j] = sqrt(pow(s1, 2) + pow(s2, 2));
        }
    }
    // Return the matrix of gradient magnitudes.
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

int main(int argc, char* argv[]) {
    /* local variable */
    int p, P, tag, n;
    double start_tot_time, end_tot_time, elapsed_tot_time;
    double start_comm_time, end_comm_time, elapsed_comm_time;
    double start_startup_time, end_startup_time, elapsed_startup_time;

    tag = 42;
    MPI_Status status;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
        
    int matrix_size[2];
    vector<double> flatten_list;
    string image_path = argv[1];
    
    
    // MPI_Barrier(MPI_COMM_WORLD);
    start_tot_time = MPI_Wtime();
    start_comm_time = 0.0;
    
    /* In the case of only one processor */
    if (P==1) {
        char filtered_path[1024] = "images/mat_filtered/filtered_1p_"; 
        strcat(filtered_path, basename(argv[1]));
        vector<vector<double>> mat = getMatrix(image_path);
        vector<vector<double>> mag = sobelStar(mat);
        // printf("hej \n");
        FILE *fp;
        fp = fopen(filtered_path, "w");
        int cnt = 0; int total_runs = mag.size() * mag[0].size(); 
		for (int i = 0; i < mag.size(); i++) {
            // cout << "Run " << cnt << "/" << total_runs << "" << endl;
            // printf("Run %d/%d, %f percent finished...\n", cnt, total_runs, (double) cnt/total_runs * 100);
            fflush(stdout);
            for (int j = 0; j < mag[i].size()-1; j++) {
		        fprintf(fp, "%f,", mag[i][j]);
                cnt++;
            }
            fprintf(fp, "%f\n", mag[i][mag[i].size()-1]);
        }
        fclose(fp);
        end_tot_time = MPI_Wtime();
        elapsed_tot_time = end_tot_time  - start_tot_time;
        printf("Process %d finished. Took %f seconds.\n", p, elapsed_tot_time);
        MPI_Finalize();
        exit(0);
    }
    
    // Master processor reads image
    if (p==0) {
        // cout << image_path << endl;
        vector<vector<double>> mat = getMatrix(image_path);

        int N = (int) mat.size(); 
        int M = (int) mat[0].size();
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

    // Initialize vectors to hold the send counts and displacements for each process
    vector<int> sendcnts(P);
    vector<int> displs(P);

    // Set the displacement of the root process to zero
    displs[0] = 0; 

    // Calculate the width and height of the image
    int W = M;
    int H = (int) N/P;
    // Calculate the remainder of N/P
    int res = N - P * H;

    start_startup_time = MPI_Wtime();
    // Calculate the send counts and displacements for each process except the root process
    for (int i = 0; i < P - 1; i++) {
        sendcnts[i] = W * (i < res ? H + 1 : H);
        displs[i+1] = displs[i] + sendcnts[i];
    }

    // Calculate the send counts and displacements for each process except the root process
    sendcnts[P - 1] = W * H;

    // Initialize a new vector to hold the received elements
    vector<double> local_flattened_list(sendcnts[p]);

    // Scatter the flattened list using MPI_Scatterv
    MPI_Scatterv(flatten_list.data(), sendcnts.data(), displs.data(), MPI_DOUBLE, 
                local_flattened_list.data(), sendcnts[p], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    elapsed_startup_time = MPI_Wtime() - start_startup_time;

    flatten_list.clear();

    H = H + (p < res ? 1 : 0); // Handles the residuals
    vector<vector<double>> local_mat = reshape(local_flattened_list, H, W);

    int W_mag = W - 2;
    int H_mag = H + (p==0 || p==P-1 ? 1 : 0); // Handles the residuals
    vector<double> mag_top(W_mag); 
    vector<double> mag_down(W_mag); 
    vector<vector<double>> mat_temp(3 , vector<double> (3, 0)); // Create a 3x3 matrix with zeros  
    int j, k;

    // Sending downwards and receiving from upwards
    for (int i = 0; i < W-2; i++) {         
        // First iteration, send an receive 3 elements                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
        start_comm_time = MPI_Wtime();

        // initialize the original temp_matrix
        if (p!=0) {
            for (j = 1; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    mat_temp[j][k] = local_mat[j-1][k+i];
                }
            }
        }

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
        elapsed_comm_time += MPI_Wtime() - start_comm_time;
        mag_top[i] = sobel3x3(mat_temp);
    }

    // Sending upwards and receiving from downwards
    for (int i = 0; i < W-2; i++) {         
        // First iteration, send an receive 3 elements                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
        start_comm_time = MPI_Wtime();

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
        elapsed_comm_time += MPI_Wtime() - start_comm_time;
        mag_down[i] = sobel3x3(mat_temp);
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

    /* output for graphical representation */
    /* Instead of using gather (which may lead to excessive memory requirements
    on the master process) each process will write its own data portion. This
    introduces a sequentialization: the hard disk can only write (efficiently)
    sequentially. Therefore, we use the following strategy:
    1. The master process writes its portion. (file creation)
    2. The master sends a signal to process 1 to start writing.
    3. Process p waites for the signal from process p-1 to arrive.
    4. Process p writes its portion to disk. (append to file)
    5. process p sends the signal to process p+1 (if it exists)
    */
    FILE *fp;
    char filtered_path[1024] = "images/mat_filtered/filtered_"; 
    strcat(filtered_path, basename(argv[1]));
    
    if (p==0){ // Master process
        fp = fopen(filtered_path, "w");
		for (int i = 0; i < mag.size(); i++) {
            for (int j = 0; j < mag[i].size()-1; j++) {
		        fprintf(fp, "%f,", mag[i][j]);
            }
            fprintf(fp, "%f\n", mag[i][j+1]);
        }
        fclose(fp);
        start_comm_time = MPI_Wtime();
        MPI_Send("hi", 2, MPI_CHAR, 1, tag, MPI_COMM_WORLD);
        elapsed_comm_time += MPI_Wtime() - start_comm_time;

    } else {
        char message[2]; // Nonesense message
        MPI_Recv(message, 2, MPI_CHAR, p-1, tag, MPI_COMM_WORLD, &status);
        fp = fopen(filtered_path, "a");
		for (int i = 0; i < mag.size(); i++) {
            for (int j = 0; j < mag[i].size()-1; j++) {
		        fprintf(fp, "%f,", mag[i][j]);
            }
            fprintf(fp, "%f\n", mag[i][j+1]);
        }
        if (p!=P-1) { // If it's not the last processor, keep sending the signals forward
            start_comm_time = MPI_Wtime();
            MPI_Send(message, 2, MPI_CHAR, p+1, tag, MPI_COMM_WORLD);
            elapsed_comm_time += MPI_Wtime() - start_comm_time;
        }
        fclose(fp);
    }
    
    // MPI_Barrier(MPI_COMM_WORLD);
    elapsed_tot_time =  MPI_Wtime() - start_tot_time;
    elapsed_comm_time += elapsed_startup_time;
    double elapsed_comp_time = elapsed_tot_time - elapsed_comm_time;

    double max_time, min_time, tot_time;
    double max_comm_time, min_comm_time, tot_comm_time;
    double max_comp_time, min_comp_time, tot_comp_time;
    double max_startup_time, min_startup_time, tot_startup_time;

    /* Inspired by https://stackoverflow.com/questions/5298739/mpi-global-execution-time */
    MPI_Reduce(&elapsed_tot_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    MPI_Reduce(&elapsed_comm_time, &max_comm_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    MPI_Reduce(&elapsed_comp_time, &max_comp_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    MPI_Reduce(&elapsed_startup_time, &max_startup_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    MPI_Finalize();

    if (p==0) {
        // printf("Process %d finished. Min time: %lf Max time: %lf Avg time: %lf, ", p, min_time, max_time, tot_time / P);
        // printf("Min comm time: %lf Max comm time: %lf Avg comm time: %lf, ", max_comm_time, min_comm_time, tot_comm_time / P);
        // printf("Min comp time: %lf Max comp time: %lf Avg comp time: %lf.\n", max_comp_time, min_comp_time, tot_comp_time / P);
        printf("Process %d finished. Total time: %lf, communication time: %lf, computation time: %lf, starup time: %lf \n", p, max_time, max_comm_time, max_comp_time, max_startup_time);
    } else {
        printf("Process %d finished.\n", p);
    }
    // /* That's it */
    exit(0);
}