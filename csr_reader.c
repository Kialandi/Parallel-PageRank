#include <stdio.h>
#include <stdlib.h>

void readfile(char *fileName, int *n, double **value, int **colind, int **rbegin){
	FILE *fp=fopen(fileName,"r");
	char buf[200];
	int i,j,k,l;
	int p,q;
	double w;
	int m;
	if ((fp=fopen(fileName,"r"))==NULL){
	  fprintf(stderr,"Open file errer, check filename please.\n");
	}
	fgets(buf,200,fp);
	*n=atoi(buf);
	l=0;
	while(buf[l++]!=' ');
	m=atoi(&buf[l]);
	printf("Matrix size: %d, #Non-zeros: %d\n",*n,m);
	(*value)=(double*)malloc(sizeof(double)*m);
	(*colind)=(int*)malloc(sizeof(int)*m);
	(*rbegin)=(int*)malloc(sizeof(int)*((*n)+1));
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
	close(fp);
}


void output(int n, double *answer){
	FILE *fp=fopen("output.dat","w");
	int i;
	for(i=0;i<n;i++){
	  fprintf(fp,"%.16f\n",answer[i]);
	}
	close(fp);
}


int main(int argc, char* argv[]){
	double *value=NULL;
	int *colind=NULL;
	int *rbegin=NULL;
	double *answer=NULL;
	int n;
	readfile(argv[1],&n,&value,&colind,&rbegin);
	// compute PageRank
	output(n,answer);
	if (value!=NULL) {free(value);value=NULL;}
	if (colind!=NULL) {free(colind);colind=NULL;}
	if (rbegin!=NULL) {free(rbegin);rbegin=NULL;}
	if (answer!=NULL) {free(answer);answer=NULL;}
}