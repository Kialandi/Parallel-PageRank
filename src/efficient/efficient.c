#include <mpi.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "../helpers.c"
#include <string.h>

#define VALUETAG    500
#define COLINDTAG   501

typedef struct {
    int processor;
    int numOfElements;
    //array of the rows i have to send
    int * rowNumbers;
    int index;
    double * sendBuf;
} sendTo;

typedef struct {
    int processor;
    int numOfElements;
    //the rows i'm getting
    int * rowNumbers;
    double * recvBuf;
} recvFrom;

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
    int ** destinationsOf = NULL;
    int rowPerProc;// = n / world_size;
    int * rbeginDiff = NULL;
    int leftover = 0;
    int leftoverEles = 0;
    double * leftoverVal = NULL;
    int * leftoverColind = NULL;
    int * leftoverRbegin = NULL;
    int end = 0;

    //arrays for each processor
    sendTo ** sendToArr = (sendTo **) malloc(world_size * sizeof(sendTo *));
    recvFrom ** recvFromArr = (recvFrom **) malloc(world_size * sizeof(recvFrom *));


    if (world_rank == 0) {
        printf("Setting up...\n");
        //proc 0 handles splitting up all the arrays to individual procs

        //number of non zeroes
        int NZ = readfile(argv[1],&n,&value,&colind,&rbegin);


        if (!NZ || n < world_size) {
            printf("No non-zeroes found or number of rows < number of processors!\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        int i;
        rowPerProc = n / world_size;
        leftover = n % world_size;

        //need an array of same dimension as number of rows
        destinationsOf = (int **) malloc(n * sizeof(int *));

        for (i = 0; i < n; i++) {
            //total number of destinations for each is just the number of procs
            destinationsOf[i] = (int *) malloc(world_size * sizeof(int));
            int j;
            for (j = 0; j < world_size; j++) {
                //destination of row is j
                destinationsOf[i][j] = -1;
            }
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

              //create a diff array

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
        
        

        int colArrIndex = 0;

        for (i = 0; i < world_size; i++) {
            //go through colind array and decide things
            //each colind we're sending out has to be processed, processor at a
            //time
            int j;
            int myLabel = i;
            for (j = 0; j < numElements[i]; j++) {
                //vector index that I need something from
                //leftover depends on colArrIndex
                int vectorIndex = colind[colArrIndex++];
                addMeTo(world_size, destinationsOf[vectorIndex], myLabel);
                //printf("Adding %d to %d's list of destinations\n", myLabel, vectorIndex);
            }
        }

        int numOfEle = 0;
        if (leftover > 0) {
            for (i = 0; i < leftover; i++) {
                numOfEle += rbeginDiff[index++];
            }
            //save the number of elements
            leftoverEles = numOfEle;


            for (i = 0; i < leftoverEles; i++) {
                int myLabel = 0;

                int vectorIndex = colind[colArrIndex++];
                addMeTo(world_size, destinationsOf[vectorIndex], myLabel);
                //printf("leftover:Adding %d to %d's list of destinations\n", myLabel, vectorIndex);
            }
        }
        /*
           for (i = 0; i < n; i++) {
           printf("Index %d destinations: ", i);
           int j;
           for (j = 0; j < world_size; j++) {
           printf("%d, ", destinationsOf[i][j]);
           }
           printf("\n");
           }
           */
        //tell each process how much space to make for the elements
        MPI_Scatter(&numElements, 1, MPI_INT, &totalEle, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //send rbegin
        MPI_Scatter(rbeginDiff, rowPerProc, MPI_INT, &(rbegin[1]), rowPerProc, MPI_INT, 0, MPI_COMM_WORLD);

        //send destinations array
        for (i = 1; i < world_size; i++) {
            int j;
            for (j = 0; j < rowPerProc; j++) {
                /*printf("%d sending %d to proc %d at tag %d\n", world_rank, j, i, i + j);
                  int z;
                  printf("This arr: ");
                  for (z = 0; z < world_size; z++) {
                  printf("%d, ", destinationsOf[(i * rowPerProc) + j][z]);
                  }
                  printf("\n");*/
                MPI_Send(destinationsOf[(i * rowPerProc) + j], world_size, MPI_INT, i, i + j, MPI_COMM_WORLD);
                //printf("send success\n");
            }
        }

        //recompute rolling sum
        for (i = 1; i < rowPerProc + 1; i++) {
            rbegin[i] += rbegin[i - 1];
        }

        valuePTR += totalEle;
        colindPTR += totalEle;



        //send value and colind
        for (i = 1; i < world_size; i++) {
            //send value array
            MPI_Send(valuePTR, numElements[i], MPI_DOUBLE, i, VALUETAG, MPI_COMM_WORLD);
            //send colind array
            MPI_Send(colindPTR, numElements[i], MPI_INT, i, COLINDTAG, MPI_COMM_WORLD);
            valuePTR += numElements[i];
            colindPTR += numElements[i];
        }

        if (leftover > 0) {
            end = 0;
            for (i = 0; i < world_size; i++) {
                end += numElements[i];
            }
            //let p0 handle it
            //pick up the leftover values

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
        int i;
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

        //set up for destinations
        destinationsOf = (int **) malloc(rowPerProc * sizeof(int *));
        for (i = 0; i < rowPerProc; i++) {
            //total number of destinations for each is just the number of procs
            destinationsOf[i] = (int *) malloc(world_size * sizeof(int));
        }

        //get rbegin
        MPI_Scatter(rbeginDiff, rowPerProc, MPI_INT, &(rbegin[1]), rowPerProc, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Status stat;
        for (i = 0; i < rowPerProc; i++) {
            //printf("%d waiting for row %d at tag %d\n", world_rank, i, world_rank + i);
            MPI_Recv(destinationsOf[i], world_size, MPI_INT, 0, world_rank + i, MPI_COMM_WORLD, &stat);
            //printf("%d received row %d's destinations\n", world_rank, i);
            //printVec(world_size, destinationsOf[i]);     
        }
        //after this loop we have an array of personal rows with a list of their
        //destinations

        //recompute rolling sum
        for (i = 1; i < rowPerProc + 1; i++) {
            rbegin[i] += rbegin[i - 1];
        }

        //get value and colind arrays
        MPI_Recv(value, totalEle, MPI_DOUBLE, 0, VALUETAG, MPI_COMM_WORLD, &stat);
        MPI_Recv(colind, totalEle, MPI_INT, 0, COLINDTAG, MPI_COMM_WORLD, &stat);        
    }
    int i;

    //printf("%d: elements: %d n: %d rows: %d\n", world_rank, totalEle, n, rowPerProc);



    MPI_Barrier(MPI_COMM_WORLD);
    //return 1;



    double * result = malloc(n * sizeof(double));
    double * leftoverResult = malloc(leftover * sizeof(double));
    double * vector = malloc(n * sizeof(double));
    double * vecAndAlpha = malloc(n * sizeof(double));
    double * temp;

    //create an expecting array, worst case receive all rows

    int * senders = (int *) malloc(world_size * sizeof(int));
    int * expectingFrom = (int *) calloc(n, sizeof(int));
    int * recvRow = (int *) calloc(n, sizeof(int));
    int * howmany = (int *) calloc(world_size, sizeof(int));

    int * rowsSendingToProc = (int *) calloc(world_size, sizeof(int));

    
    for (i = 0; i < rowPerProc; i++) {
        //for each row figure out which procs need it
        int j;
        for (j = 0; j < world_size; j++) {
            if (destinationsOf[i][j] == -1) {
                //found the end of the list
                break;
            }
            int recipient = destinationsOf[i][j];

            rowsSendingToProc[recipient]++;
        }

    }
    int leftoverFirstRow = world_size * rowPerProc;
    
    if (leftover > 0 && world_rank == 0) {
        for (i = leftoverFirstRow; i < leftoverFirstRow + leftover; i++) {
            int j;
            for (j = 0; j < world_size; j++) {
                if (destinationsOf[i][j] == -1) {
                    //found the end of the list
                    break;
                }
                int recipient = destinationsOf[i][j];

                rowsSendingToProc[recipient]++;
            }
        }
    }



    for (i = 0; i < world_size; i++) {        
        if (world_rank == i || rowsSendingToProc[i] == 0) {
            sendToArr[i] = NULL;
            //printf("%d doesn't need to send %d anything\n", world_rank, i);
            continue;
        }
        //printf("%d needs to send %d %d rows\n", world_rank, i, rowsSendingToProc[i]);
        sendToArr[i] = (sendTo *) malloc(sizeof(sendTo));

        sendToArr[i]->processor = i;
        sendToArr[i]->numOfElements = rowsSendingToProc[i];
        //index used to figure out if I filled out the row numbers
        sendToArr[i]->index = 0;
        sendToArr[i]->rowNumbers = (int *) malloc(rowsSendingToProc[i] * sizeof(int));
        sendToArr[i]->sendBuf = (double *) malloc(rowsSendingToProc[i] * sizeof(double));
    }

    //all values in the array that are not -1 are destination processors
    for (i = 0; i < rowPerProc; i++) {
        int j;

        //printf("%d has to send row %d to: ",world_rank,  i + world_rank * rowPerProc);

        for (j = 0; j < world_size; j++) {
            //skip yourself
            if (destinationsOf[i][j] == world_rank) continue;
            if (destinationsOf[i][j] == -1) break;
            int recipient = destinationsOf[i][j];
            //printf("%d, ", recipient);
            sendToArr[recipient]->rowNumbers[sendToArr[recipient]->index++] = i + world_rank * rowPerProc;
        }
        //printf("\n");
    }

    if (leftover > 0 && world_rank == 0) {
        //p0 handle the extra rows
        for (i = leftoverFirstRow; i < leftoverFirstRow + leftover; i++) {
            int j;

            //printf("leftovers: %d has to send row %d to: \n",world_rank,  i);
            //    printf("dest: %d\n", destinationsOf[i][0]);

            for (j = 0; j < world_size; j++) {
                //skip yourself
                if (destinationsOf[i][j] == world_rank) continue;
                if (destinationsOf[i][j] == -1) break;

                int recipient = destinationsOf[i][j];
                //printf("%d, ", recipient);

                sendToArr[recipient]->rowNumbers[sendToArr[recipient]->index++] = i;
            }
            //printf("\n");
        }
    }

    //each element we are in charge of
    for (i = 0; i < totalEle; i++) {
        //colind[i] means the ith element's column index
        //expecting a certain row, 1 means yes 0 means no
        expectingFrom[colind[i]] = 1;
    }
    
    if (leftover > 0 && world_rank == 0) {
        for (i = 0; i < leftoverEles; i++) {
            expectingFrom[colind[end + i]] = 1;
        }
    }

    for (i = 0; i < n; i++) {
        //set everything to -1
        recvRow[i] = -1;
    }

    for (i = 0; i < world_size; i++) {
        senders[i] = -1;
    }

    int tail = 0;
    //int lowerBound = world_rank * rowPerProc;
    //int upperBound = world_size * rowPerProc;
    int tally = 0;
    int first = 0;
    for (i = 0; i < n; i++) {
        int owner = i / rowPerProc;
        //skip over yourself since you own it

        if (owner == world_rank) continue;
        if (expectingFrom[i] == 1) {
            //need row i
            if (owner + 1 > world_size) {
                owner = 0;
                if (owner == world_rank) continue;

                if (first == 0) first = tail;

                tally++;
            }
            if (owner == world_rank) continue;

            recvRow[tail++] = i;
            //if (world_rank == 0)
            //printf("%d expect row %d from %d tally: %d\n", world_rank, i, owner, tally);
        }
    }



    for (i = 0; i < tail; i++) {
        //count up how many rows for each proc
        int owner = recvRow[i] / rowPerProc;
        if (owner + 1 > world_size) {
            owner = 0;
            if (0 == world_rank) continue;
        }
        howmany[owner]++;

        //if (world_rank == 1)
        //printf("%d wants row %d from %d\n", world_rank, recvRow[i], owner);

    }

    //sort recvRow by proc

    if (leftover > 0 && tally > 0) {
        int * recvRowTemp = (int *) malloc(n * sizeof(int));
        int recvRowTempIndex = 0;
        //this is the last row not counting leftover p0 is in charge of
        int lastRow = rowPerProc - 1;
        int upperBound = world_size * rowPerProc;
        int flag = 0;
        for (i = 0; i < tail; i++) {

            if (recvRow[i] > lastRow && flag == 0) {
                flag = 1;
                int z;
                for (z = 0; z < tally; z++) {
                    recvRowTemp[recvRowTempIndex++] = recvRow[first + z];
                }
            }
            if (recvRow[i] >= upperBound) break;
            recvRowTemp[recvRowTempIndex++] = recvRow[i];

        }

        //swap
        int * temp = recvRowTemp;
        recvRowTemp = recvRow;
        recvRow = temp;
        free(recvRowTemp);
    }


    //all this is to figure out what rows you want and from what procs
    int recvRowIndex = 0;
    for (i = 0; i < world_size; i++) {
        if (i == world_rank || howmany[i] == 0) { 
            recvFromArr[i] = NULL;
            continue;
        }
        //however many rows you expect from each proc
        //allocate how much space you need for things you expect from
        //this is where you store the incoming arrays

        //make a struct for each processor you expect rows from
        recvFromArr[i] = (recvFrom *) malloc(sizeof(recvFrom));

        recvFromArr[i]->processor = i;
        recvFromArr[i]->numOfElements = howmany[i];
        //i need this many rows, adding labels for them
        recvFromArr[i]->rowNumbers = (int *) malloc(howmany[i] * sizeof(int));
        int j;
        for (j = 0; j < howmany[i]; j++) {
            if (leftover > 0 && i == 0) { 

            }
            //add in the labels
            recvFromArr[i]->rowNumbers[j] = recvRow[recvRowIndex++];
        }

        //store stuff here
        recvFromArr[i]->recvBuf = (double *) malloc(howmany[i] * sizeof(double));
        /*
           if (world_rank == 3) {
           printf("%d expects %d rows from %d\n", world_rank, howmany[i], i);
           printf("they are rows: ");
           for (j = 0; j < howmany[i]; j++) {
           int rowNum = recvFromArr[i]->rowNumbers[j];
           printf("%d, ", rowNum);
           }
           printf("\n");
           }
           */
    }
    /* 
       if(world_rank == 0) {
       for (i = 0; i < world_size; i++) {
       if (sendToArr[i] == NULL) {
       printf("%d sends nothing to %d\n", world_rank ,i);
       continue;
       }
       int numOfRowsSent = sendToArr[i]->numOfElements;
       printf("%d needs to send %d %d rows\n", world_rank, i, numOfRowsSent);
       printf("they are rows: ");
       int j;
       for (j = 0; j < numOfRowsSent; j++) {
       int rowNum = sendToArr[i]->rowNumbers[j];
       printf("%d, ", rowNum);
       }
       printf("\n");
       }
       }*/

    //printf("%d expects %d total rows\n", world_rank ,tail);

    for (i = 0; i < n; i++) {
        //init for k_0
        vector[i] = 1.0 / n;
    }

    //compute PageRank
    int iterations = 0;
    double alpha = 0.85;
    double epsiMax = 0; 
    //location in result vector to write into
    int location = world_rank * rowPerProc;
    double * locationPTR = &(result[location]);
    if (world_rank == 0) {
        printf("Done! Starting Page Rank Computation\n");
        if (leftover > 0) 
            printf("There's leftovers\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //start timer
    if (world_rank == 0) {
        gettimeofday(&tz, &tx);
        start_time = (double)tz.tv_sec + (double) tz.tv_usec / 1000000.0;
    }
    int leftoverLoc = world_size * rowPerProc; 
    
    
    while(1) {
        //result <- alpha * matrix * k
        //csrMatVec(rowPerProc, alpha, value, colind, rbegin, vecAndAlpha, result); 
        if (world_rank == 1) {
            //printf("Iter %d Vector:\n", iterations);
            //printVecDoubles(n, vector);
        }
        //        csrMatVec(rowPerProc, alpha, value, colind, rbegin, vector, result); 
        csrMatVec(rowPerProc, alpha, value, colind, rbegin, vector, &(result[location])); 
        if (leftover > 0 && world_rank == 0) {
            //matvec for leftover array at p0
            csrMatVec(leftover, alpha, leftoverVal, leftoverColind, leftoverRbegin, vector, &result[leftoverLoc]); 
        }
        /*
           if (world_rank == 1) {
           printf("result Iter %d Vector:\n", iterations);
           printVecDoubles(n, result);
        //printVecDoubles(n, vector);
        }*/
        iterations++;

        //result = partial y
        //sum y locally
        double sum = sumVec(rowPerProc, &(result[location]));
        double totalSum;
        if (leftover > 0 && world_rank == 0) {
            sum += sumVec(leftover, &result[leftoverLoc]);
        }


        //reduce and broadcast somehow is faster
        //MPI_Allreduce(&sum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Reduce(&sum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Bcast(&totalSum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        double gamma = 1.0 - totalSum;
        double frac = gamma / n;

        //add frac to local partial result
        addToVec(rowPerProc, &result[location], frac);
        if (leftover > 0 && world_rank == 0) {
            addToVec(leftover, &result[leftoverLoc], frac);
        }

        double epsi;
        //compute local convergence
        epsi = computeChange(rowPerProc, &result[location], &vector[location]);
        if (leftover > 0 && world_rank == 0) {
            double leftoverEpsi = computeChange(leftover, &result[leftoverLoc], &vector[leftoverLoc]);
            if (epsi < leftoverEpsi)
                epsi = leftoverEpsi;
        }


        MPI_Allreduce(&epsi, &epsiMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
             /*   if (world_rank == 0) {
                  printf("epsimax:%lf\n", epsiMax);
                  }*/
                  
        MPI_Request req;
        /*******************************************************************************************
         *
         *                              Sending loop
         *
         * *******************************************************************************************/

        //a bunch of non blocking sends to all the destinations
        for (i = 0; i < world_size; i++) {
            if (sendToArr[i] == NULL) continue;

            int j;
            for (j = 0; j < sendToArr[i]->numOfElements; j++) {
                //sending numOfElements rows to proc i, so copy things over into
                //the buffer
                int realRow = sendToArr[i]->rowNumbers[j];
                //int localRow = realRow - (world_rank * rowPerProc);
                memcpy(&sendToArr[i]->sendBuf[j], &result[realRow] , sizeof(double));
            }
            //non blocking send
            MPI_Isend(sendToArr[i]->sendBuf, sendToArr[i]->numOfElements, MPI_DOUBLE, i, world_rank, MPI_COMM_WORLD, &req);
            //if (world_rank == 3)
            //printf("%d sending %d rows to %d\n", world_rank, sendToArr[i]->numOfElements, i);
        }
    
        /*******************************************************************************************
         *
         *                              Receiving loop
         *
         * *******************************************************************************************/
        //blocking receives for all incoming rows
        MPI_Status stat;
        for (i = 0; i < world_size; i++) {
            if (recvFromArr[i] == NULL) continue;

            //if (world_rank == 2)
                //printf("%d receiving %d rows from %d\n", world_rank, recvFromArr[i]->numOfElements, i);
            MPI_Recv(recvFromArr[i]->recvBuf, recvFromArr[i]->numOfElements, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &stat);

            int j;
            for (j = 0; j < recvFromArr[i]->numOfElements; j++) {
                int spot = recvFromArr[i]->rowNumbers[j];
                memcpy(&result[spot], &recvFromArr[i]->recvBuf[j], sizeof(double));
            }

        }
   
        MPI_Barrier(MPI_COMM_WORLD);
  
        if (epsiMax < 0.00001) break;
        if (iterations > 100) { break; }

        //swap the arrays
        temp = vector;
        vector = result;
        result = temp;
    }

    if (world_rank == 0) {
        //stop timer
        gettimeofday(&tz, &tx);
        end_time = (double)tz.tv_sec + (double) tz.tv_usec / 1000000.0;
        total_time = end_time - start_time;
    }

    MPI_Allgather(&result[location], rowPerProc, MPI_DOUBLE, vecAndAlpha, rowPerProc, MPI_DOUBLE, MPI_COMM_WORLD);
    
    if (leftover > 0) {
        memcpy(&vecAndAlpha[leftoverLoc], &result[leftoverLoc], leftover * sizeof(double));
    }
    if (world_rank == 0) {
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
