#include <mpi.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "../helpers.c"
#include <string.h>

int main(int argc, char* argv[]){
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    // Find out rank, size
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); 
    
    //data structures for timing
    double start_time, end_time, total_time;
    struct timeval tz;
    struct timezone tx;

    double *value=NULL;
    int *colind=NULL;
    int *rbegin=NULL;
    double *answer=NULL;
    int n;
    int totalEle;
    int numElements[world_size];
    int rowPerProc;// = n / world_size;
    int * rbeginDiff = NULL;
    int leftover = 0;
    int leftoverEles = 0;
    double * leftoverVal = NULL;
    int * leftoverColind = NULL;
    int * leftoverRbegin = NULL;

    if (world_rank == 0) {
        printf("Setting up...\n");
        //proc 0 handles splitting up all the arrays to individual procs
        
        //number of non zeroes
        int NZ = readfile(argv[1],&n,&value,&colind,&rbegin);

        if (!NZ || n < world_size) {
            printf("No non-zeroes found or number of rows < number of processors!\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
/*        
        if (n <= 6) {
            //DEBUG
            printf("=========== ORIGINAL MATRIX =============\n");
            printf("vals\n");
            printVecDoubles(NZ, value);
            printf("colind\n");
            printVec(NZ, colind);
            printf("rbegin\n");
            printVec(n + 1, rbegin);
            printf("=========== ORIGINAL MATRIX =============\n");
        }
*/
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

        rbeginDiff = (int *) malloc(n * sizeof(int));
        double * valuePTR = value;
        int * colindPTR = colind;

        rowPerProc = n / world_size;
        leftover = n % world_size;

        //create a diff array
        int i;
        for (i = n; i > 0; i--) {
            rbeginDiff[i - 1] = rbegin[i] - rbegin[i - 1];
        }
        
        //printf("rbeginDiff:\n");
        //printVec(n, rbeginDiff);
        //count up element for each proc
        
        int index = 0;
        for (i = 0; i < world_size; i++) {
            int j;
            int numOfEle = 0;
            for (j = 0; j < rowPerProc; j++) {
                //leftover calc relies on value of index
                numOfEle += rbeginDiff[index++];
            }
            numElements[i] = numOfEle;
        }
        
        //printf("num of elements\n");
        //printVec(world_size, numElements);

        //tell each process how much space to make for the elements
        MPI_Scatter(&numElements, 1, MPI_INT, &totalEle, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //send rbegin
        MPI_Scatter(rbeginDiff, rowPerProc, MPI_INT, &(rbegin[1]), rowPerProc, MPI_INT, 0, MPI_COMM_WORLD);

        //recompute rolling sum
        for (i = 1; i < rowPerProc + 1; i++) {
            rbegin[i] += rbegin[i - 1];
        }

        valuePTR += totalEle;
        colindPTR += totalEle;

        //TODO: send value and colind
        for (i = 1; i < world_size; i++) {
            //send value array
            MPI_Send(valuePTR, numElements[i], MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
            //send colind array
            MPI_Send(colindPTR, numElements[i], MPI_INT, i, 2, MPI_COMM_WORLD);
            valuePTR += numElements[i];
            colindPTR += numElements[i];
        }

        if (leftover > 0) {
            int end = 0;
            for (i = 0; i < world_size; i++) {
                end += numElements[i];
            }
            //let p0 handle it
            //pick up the leftover values
            //printf("there's %d extra rows\n", n % world_size);
            int numOfEle = 0;
            
            for (i = 0; i < leftover; i++) {
                numOfEle += rbeginDiff[index++];
            }
            //leftoverVal = (double *) malloc(sizeof(double));
            //save the number of elements
            leftoverEles = numOfEle;
            
            double * leftoverValPTR = value + end;
            leftoverVal = (double *) malloc(numOfEle * sizeof(double));
            memcpy(leftoverVal, leftoverValPTR, numOfEle * sizeof(double));

            int * leftoverColindPTR = colind + end;
            leftoverColind = (int *) malloc(numOfEle * sizeof(int));
            memcpy(leftoverColind, leftoverColindPTR, numOfEle * sizeof(int));
            
            int * leftoverRbeginPTR = rbeginDiff + world_size * rowPerProc;
            leftoverRbegin = (int *) malloc((leftover + 1) * sizeof(int));
            leftoverRbegin[0] = 0;
            memcpy(leftoverRbegin + 1, leftoverRbeginPTR, leftover * sizeof(int));

            //recompute rolling sum
            for (i = 1; i < leftover + 1; i++) {
                leftoverRbegin[i] += leftoverRbegin[i - 1];
            }
/*
            printf("==================== LEFTOVER SANITY CHECK ==================\n");
            printf("%d: vals\n", world_rank);
            printVecDoubles(numOfEle, leftoverVal);
            printf("%d: colind\n", world_rank);
            printVec(numOfEle, leftoverColind);
            printf("%d: rbegin\n", world_rank);
            printVec(leftover + 1, leftoverRbegin);
            printf("==================== LEFTOVER SANITY CHECK ==================\n");
*/
        }
    }
    else {
        //non p0, wait for all the stuff
        //get the dimensions of the matrix
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        rowPerProc = n / world_size;
        leftover = n % world_size;

        //printf("%d: rows: %d\n", world_rank, rowPerProc);

        MPI_Scatter(&numElements, 1, MPI_INT, &totalEle, 1, MPI_INT, 0, MPI_COMM_WORLD);

        //create buffers
        value = (double *) malloc(totalEle * sizeof(double));
        colind = (int *) malloc(totalEle * sizeof(int));
        rbegin = (int *) malloc((rowPerProc + 1) * sizeof(int));
        rbegin[0] = 0;

        //get rbegin
        MPI_Scatter(rbeginDiff, rowPerProc, MPI_INT, &(rbegin[1]), rowPerProc, MPI_INT, 0, MPI_COMM_WORLD);

        //printf("%d: got new rbegin\n", world_rank);

        //recompute rolling sum
        int i;
        for (i = 1; i < rowPerProc + 1; i++) {
            rbegin[i] += rbegin[i - 1];
        }

        MPI_Status stat;
        //get value and colind arrays
        MPI_Recv(value, totalEle, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &stat);
        MPI_Recv(colind, totalEle, MPI_INT, 0, 2, MPI_COMM_WORLD, &stat);        
    }

    //printf("%d: elements: %d n: %d rows: %d\n", world_rank, totalEle, n, rowPerProc);
    //sanity check
    /*
       printf("==================== SANITY CHECK ==================\n");
       printf("%d: vals\n", world_rank);
       printVecDoubles(totalEle, value);
       printf("%d: colind\n", world_rank);
       printVec(totalEle, colind);
       printf("%d: rbegin\n", world_rank);
       printVec(rowPerProc + 1, rbegin);
       printf("==================== SANITY CHECK ==================\n");
       */
    MPI_Barrier(MPI_COMM_WORLD);
    double * result = malloc(rowPerProc * sizeof(double));
    double * leftoverResult = malloc(leftover * sizeof(double));
    double * vector = malloc(n * sizeof(double));
    double * vecAndAlpha = malloc(n * sizeof(double));
    double * temp;
    int i;
    for (i = 0; i < n; i++) {
        //init for k_0
        vector[i] = 1.0 / n;
    }

    //compute PageRank
    int iterations = 0;
    double alpha = 0.85;
    double epsiMax = 0; 

    if (world_rank == 0) {
        printf("Done! Starting Page Rank Computation\n");
        if (leftover > 0) 
            printf("There's leftovers\n");
    }
    //start timer
    if (world_rank == 0) {
        gettimeofday(&tz, &tx);
        start_time = (double)tz.tv_sec + (double) tz.tv_usec / 1000000.0;
    }
    
    while(1) {
        //result <- alpha * matrix * k
        //csrMatVec(rowPerProc, alpha, value, colind, rbegin, vecAndAlpha, result); 
        csrMatVec(rowPerProc, alpha, value, colind, rbegin, vector, result); 
        if (leftover > 0 && world_rank == 0) {
            //matvec for leftover array at p0
            csrMatVec(leftover, alpha, leftoverVal, leftoverColind, leftoverRbegin, vector, leftoverResult); 
        }
        
        iterations++;

        //result = partial y
        //sum y locally
        double sum = sumVec(rowPerProc, result);
        double totalSum;
        if (leftover > 0 && world_rank == 0) {
            sum += sumVec(leftover, leftoverResult);
        }

        //reduce and broadcast somehow is faster
        //MPI_Allreduce(&sum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Reduce(&sum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Bcast(&totalSum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        double gamma = 1.0 - totalSum;
        double frac = gamma / n;
        
        //add frac to local partial result
        addToVec(rowPerProc, result, frac);
        if (leftover > 0 && world_rank == 0) {
            addToVec(leftover, leftoverResult, frac);
        }
       
        double epsi;
        //compute local convergence
        epsi = computeChange(rowPerProc, result, &vector[rowPerProc * world_rank]);
        if (leftover > 0 && world_rank == 0) {
            double leftoverEpsi = computeChange(leftover, leftoverResult, &vector[rowPerProc * world_size]);
            if (epsi < leftoverEpsi)
                epsi = leftoverEpsi;
        }
        
        MPI_Allreduce(&epsi, &epsiMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        //all to all bcast the partial vectors and store in vecAndAlpha
        MPI_Allgather(result, rowPerProc, MPI_DOUBLE, vecAndAlpha, rowPerProc, MPI_DOUBLE, MPI_COMM_WORLD);
  
        if (leftover > 0) {
            //broadcast the extra bits of the vector
            MPI_Bcast(leftoverResult, leftover, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            //copy it over to result array
            memcpy(vecAndAlpha + world_size * rowPerProc, leftoverResult, leftover * sizeof(double));
        }
         
        if (epsiMax < 0.00001) break;
        if (iterations > 100) { printf("iter over 100\n"); break; }

        if (world_rank == 0 && n == 6) {
            printf("Result of iter %d: ", iterations);
            printVecDoubles(n, vecAndAlpha);
        }

        //swap the arrays
        temp = vector;
        vector = vecAndAlpha;
        vecAndAlpha = temp;
    }

    if (world_rank == 0) {
        //stop timer
        gettimeofday(&tz, &tx);
        end_time = (double)tz.tv_sec + (double) tz.tv_usec / 1000000.0;
        total_time = end_time - start_time;

        printf("%d power iter conv at %lf\n", iterations, epsiMax);
        if (n == 6) {
            printf("Final pagerank vector\n");
            printVecDoubles(n, vecAndAlpha);
        }

        //print the output to file 
        printf("Writing to file, please be patient... ");
        output(n, vecAndAlpha);
        printf("Done!\n");
    }

    if (value!=NULL) {free(value);value=NULL;}
    if (colind!=NULL) {free(colind);colind=NULL;}
    if (rbegin!=NULL) {free(rbegin);rbegin=NULL;}
    if (answer!=NULL) {free(answer);answer=NULL;}
    if (result!=NULL) {free(result);result=NULL;}
    if (leftoverResult!=NULL) {free(leftoverResult);leftoverResult=NULL;}
    if (vector!=NULL) {free(vector);vector=NULL;}
    if (vecAndAlpha!=NULL) {free(vecAndAlpha);vecAndAlpha=NULL;}
    if (leftoverVal!=NULL) {free(leftoverVal);leftoverVal=NULL;}
    if (leftoverColind!=NULL) {free(leftoverColind);leftoverRbegin=NULL;}
    if (leftoverRbegin!=NULL) {free(leftoverRbegin);leftoverRbegin=NULL;}

    if (world_rank == 0) {
        printf("Elapsed time = %lf seconds\n", total_time);
    }

    MPI_Finalize();

    return 0;
}
