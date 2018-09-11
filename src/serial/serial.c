#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

double myabsol(double x) {
    return x < 0 ? x * - 1 : x;
}

void multVec(int num, double * vector, double val, double * answer) {
    int i;
    for (i = 0; i < num; i++) {
        answer[i] = vector[i] * val;
    }
}

void addVecs(int num, double * vector1, double * vector2) {
    int i;
    for (i = 0; i < num; i++) {
        vector1[i] += vector2[i];
    }
}

double computeChange(int num, double * vector1, double * vector2) {
    int i;
    double max = myabsol(vector1[0] - vector2[0]);

    for (i = 1; i < num; i++) {
        double diff = myabsol(vector1[i] - vector2[i]);
        if (diff > max) 
            max = diff;
    }

    //printf("computing change, first max: %lf\n", max);
    return max;
}

double sumVec(int num, double * vector) {
    int i;
    double res = 0;
    for (i = 0; i < num; i++) {
        res += vector[i];
    }

    return res;
}

void printVec(int num, int * vector) {
    int i;
    for (i = 0; i < num; i++) {
        printf("%d, ", vector[i]);
    }
    printf("\n");
}

void printVecDoubles(int num, double * vector) {
    int i;
    for (i = 0; i < num; i++) {
        if (i == num - 1)
            printf("%lf ", vector[i]);
        else
            printf("%lf, ", vector[i]);
    }
    printf("\n");
}

void csrMatVec(int n, double alpha, double * values, int * colindex, int * rbegin, double * vector, double * result) {
    int i;
    int k1, k2, k;
    //printf("num: %d\n", num);
    //printVec(num, rbegin);
    for (i = 0; i < n; i++) {
        //zero out result vector every iteration
        result[i] = 0;
        k1 = rbegin[i];
        k2 = rbegin[i + 1] - 1;
        if (k2 < k1) {
            continue;
        }
        int j;
        for (k = k1; k <= k2; k++) {
            //printf("k: %d, k2 %d", k, k2);
            //printf("%d: inner %d: res: %lf\n", i, j, result[i]);
            //printf("k: %d, k2 %d", k, k2);
            j = colindex[k];
            result[i] += values[k] * vector[j] * alpha;
        }
        //printf("%d: result: %lf\n", i, result[i]);
    }
    /*
       printf("vector: ");
       printVecDoubles(n, vector);
       printf("result: ");
       printVecDoubles(n, result);
       printf("\n");
       */
}

int readfile(char *fileName, int *n, double **value, int **colind, int **rbegin){
    FILE *fp=fopen(fileName,"r");
    char buf[200];
    int i,j,k,l;
    int p,q;
    double w;
    int m;
    if ((fp=fopen(fileName,"r"))==NULL){
        fprintf(stderr,"Open file errer, check filename please.\n");
        return 0;
    }
    fgets(buf,200,fp);
    *n=atoi(buf);
    l=0;
    while(buf[l++]!=' ');
    m=atoi(&buf[l]);
    printf("Matrix size: %d, #Non-zeros: %d\n",*n,m);
    (*value)=(double*)malloc(sizeof(double)*m);
    (*colind)=(int*)malloc(sizeof(int)*m);
    //(*rbegin)=(int*)malloc(sizeof(int)*((*n)+1));
    (*rbegin)=(int*)calloc((*n) + 1, sizeof(int));

    k=-1;
    for(i=0;i<m;i++){
        fgets(buf,200,fp);
        l=0;p=atoi(&buf[l]);
        while(buf[l++]!=' '); q=atoi(&buf[l]);
        while(buf[l++]!=' '); w=atof(&buf[l]);
        (*value)[i]=w;
        (*colind)[i]=q;
        if(p!=k){
            k=p;
            (*rbegin)[p]=i;
        }
    }
    (*rbegin)[*n]=m;
    fclose(fp);

    return m;
}


void output(int n, double *answer){
    FILE *fp=fopen("output.dat","w");
    int i;
    for(i=0;i<n;i++){
        fprintf(fp,"%.16f\n",answer[i]);
    }
    fclose(fp);
}

double max(int num, double * arr) {
    int i = 0;
    double max = arr[0];

    for (i = 1; i < num; i++) {
        if (arr[i] > max) max = arr[i];
    }
    return max;
}

int main(int argc, char* argv[]){
    double *value=NULL;
    int *colind=NULL;
    int *rbegin=NULL;
    double *answer=NULL;
    int n;

    //number of non zeroes
    int NZ = readfile(argv[1],&n,&value,&colind,&rbegin);

    if (!NZ) {
        printf("No non-zeroes found!\n");
        return -1;
    }

    // data structures for timing
    double start_time, end_time, total_time;
    struct timeval tz;
    struct timezone tx;

    //create a vector initialized to 1/n
    double * result = malloc(n * sizeof(double));
    double * vector = malloc(n * sizeof(double));
    double * vecAndAlpha = malloc(n * sizeof(double));
    double * e = malloc(n * sizeof(double));
    double * temp;
    int i;

    for (i = 0; i < n; i++) {
        //1/n in result vector since they're swapped initially on the loop
        result[i] = 1.0 / n;
    }

    //compute PageRank
    int iterations = 0; //= atoi(argv[2]);
    double alpha = 0.85;
    double epsi = 0; 
    int z = 0;

    //start timer
    gettimeofday(&tz, &tx);
    start_time = (double)tz.tv_sec + (double) tz.tv_usec / 1000000.0;

    while(1) {

        temp = result;
        result = vector;
        vector = temp;

        //multVec(n, vector, alpha, vecAndAlpha);

        //printf("b4 csrmatvec\n"); 
        csrMatVec(n, alpha, value, colind, rbegin, vector, result);
        iterations++;
        //at this point, result = y in pseudocode
        //vector is still = x_k

        //printf("after csrmatvec\n"); 
        double gamma = 1.0 - sumVec(n, result);
        //printf("totalsum%d: %lf\n", iterations, 1.0 - gamma);
        

        double frac = gamma / n;
        for (i = 0; i < n; i++) {
            e[i] = frac;
        }

        //add the fraction to results, store in result
        addVecs(n, result, e);

        //at this point, result = x_k+1, vector is x_k
        epsi = computeChange(n, result, vector);
        //printf("epsimax%d: %lf\n", iterations, epsi);

        if (n <= 6) {
            //printf("Result of iter %d: ", iterations);
            //printVecDoubles(n, result);
        }

        if (epsi < 0.00001) break;

        if (iterations > 100) { printf("iter over 100\n"); break; }
    }

    //stop timer
    gettimeofday(&tz, &tx);
    end_time = (double)tz.tv_sec + (double) tz.tv_usec / 1000000.0;
    total_time = end_time - start_time;

    printf("%d power iter conv at %lf\n", iterations, epsi);
    if (n <= 6) {
        printf("Final pagerank vector\n");
        printVecDoubles(n, result);
    }
    //print the output to file
    output(n, result);

    if (value!=NULL) {free(value);value=NULL;}
    if (colind!=NULL) {free(colind);colind=NULL;}
    if (rbegin!=NULL) {free(rbegin);rbegin=NULL;}
    if (result!=NULL) {free(result);result=NULL;}
    if (vector!=NULL) {free(vector);vector=NULL;}
    if (e!=NULL) {free(e);e=NULL;}
    if (vecAndAlpha!=NULL) {free(vecAndAlpha);vecAndAlpha=NULL;}

    printf("Elapsed time = %lf seconds\n", total_time);

    return 0;
}
