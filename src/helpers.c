
void addMeTo(int size, int * destinations, int label) {
    int i;
    for (i = 0; i < size; i++) {
        if (destinations[i] == label) break;
        if (destinations[i] != -1) 
            continue;
        destinations[i] = label;
        break;
    }
}


double myabsol(double x) {
    return x < 0 ? x * - 1 : x;
}

void multVec(int num, double * vector, double val, double * answer) {
    int i;
    for (i = 0; i < num; i++) {
        answer[i] = vector[i] * val;
    }
}

void addToVec(int num, double * vector, double add) {
    int i;
    for (i = 0; i < num; i++) {
        vector[i] += add;
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
        if (i == num - 1) 
            printf("%d", vector[i]);
        else
            printf("%d, ", vector[i]);
    }
    printf("\n");
}

void printVecDoubles(int num, double * vector) {
    int i;
    for (i = 0; i < num; i++) {
        if (i == num - 1) 
            printf("%lf", vector[i]);
        else
            printf("%lf, ", vector[i]);
    }
    printf("\n");
}

void csrMatVec(int n, double alpha, double * values, int * colindex, int * rbegin, double * vector, double * result) {
    int i;
    int k1, k2, k;
    
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
            j = colindex[k];
            result[i] += values[k] * vector[j] * alpha;
        }
    }
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

