#define MAX 60  //nb_channel
#define KMAX 73 //pslength
#include<time.h>
#include<stdlib.h>
#include<stdio.h>

typedef struct _COMPLEX_FLOAT_
{
	float real;
	float imag;
} COMPLEX_FLOAT;

void AddNormalizeCorrelation(
		//int pslength, int nb_channel, int sum_count,
		float X_real[MAX*MAX],float X_imag[MAX*MAX],int UsesFreq[KMAX],
		float R_real[KMAX*MAX*MAX],float R_imag[KMAX*MAX*MAX],char startt[1],char stopt[1])
{

#pragma HLS INTERFACE s_axilite port=X_real
#pragma HLS INTERFACE s_axilite port=X_imag
#pragma HLS INTERFACE s_axilite port=UsesFreq
#pragma HLS INTERFACE s_axilite port=R_real
#pragma HLS INTERFACE s_axilite port=R_imag
#pragma HLS INTERFACE s_axilite port=return
#pragma HLS INTERFACE axis port=startt
#pragma HLS INTERFACE axis port=stopt

	int sum_count = 1;
    int i,j,k;
    float tmp1[MAX],tmp2[MAX];
    COMPLEX_FLOAT X[MAX][MAX];
    COMPLEX_FLOAT Rxx[KMAX][MAX][MAX];
    COMPLEX_FLOAT x_i,x_j;
//#pragma HLS ARRAY_PARTITION variable=tmp1 complete dim=0
//#pragma HLS ARRAY_PARTITION variable=tmp2 complete dim=0
//#pragma HLS ARRAY_PARTITION variable=Rxx complete dim=0

    roop1:
    for(k=0; k<KMAX; k++){
    	for(i=0; i<MAX; i++){
    		for(j=0; j<MAX; j++){
    			Rxx[k][i][j].real = 0.0f;
    			Rxx[k][i][j].imag = 0.0f;
    		}
    		Rxx[k][i][i].real = 1.0f;
    	}
    }

    //Input
    roop2:
    for (i=0;i<MAX;i++){
    	for(j=0;j<MAX;j++){
    		X[i][j].real = X_real[i*MAX+j];
    		X[i][j].imag = X_imag[i*MAX+j];
    	}
    }

    startt[0] = (char)X[MAX-1][MAX-1].imag;

    const float sum_count_1 = 1.0f / (float)sum_count;
    roop3:
    for(k=0; k<KMAX; k++){
        //if(UsesFreq[k]){
            for(i=0; i<MAX; i++){
            	x_i = X[i][k];
                for(j=0; j<MAX; j++){
//#pragma HLS UNROLL
                	x_j = X[j][k];
                	tmp1[j]  =  x_i.real * x_j.real;
                	tmp1[j] +=  x_i.imag * x_j.imag;
                	tmp1[j] +=  Rxx[k][i][j].real;
                	tmp1[j] *=  sum_count_1;
                    Rxx[k][i][j].real = tmp1[j];
                    tmp2[j]  = -x_i.real * x_j.imag;
                    tmp2[j] +=  x_i.imag * x_j.real;
                    tmp2[j] +=  Rxx[k][i][j].imag;
                    tmp2[j] *=  sum_count_1;
                    Rxx[k][i][j].imag = tmp2[j];
                }
            }
        //}
    }

    stopt[0] = (char)Rxx[KMAX-1][MAX-1][MAX-1].imag;

    //Output
    roop4:
    for(k=0;k<KMAX;k++){
    	for(i=0;i<MAX;i++){
        	for(j=0;j<MAX;j++){
        		R_real[k*MAX*MAX+i*MAX+j] = Rxx[k][i][j].real;
        		R_imag[k*MAX*MAX+i*MAX+j] = Rxx[k][i][j].imag;
        	}
        }
    }
}

int main(){
	clock_t start, end;
	int i;
	int UsesFreq[KMAX];
	float X_real[MAX*MAX],X_imag[MAX*MAX],R_real[KMAX*MAX*MAX],R_imag[KMAX*MAX*MAX];
	char startt[1], stopt[1];
	srand((unsigned)time(NULL));
	//start = clock();
	for(i = 0;i<KMAX;i++){
		//UsesFreq[i] = rand()%100 + 1;
		UsesFreq[i] = 1;
	}
	for(i = 0;i<MAX*MAX;i++){
		//X_real[i] = rand()%100 + 1;
		//X_imag[i] = rand()%100 + 1;
		X_real[i] = 1;
		X_imag[i] = 2;
	}
	start = clock();
	AddNormalizeCorrelation(X_real,X_imag,UsesFreq,R_real,R_imag,startt,stopt);
	end = clock();
	printf("time=%f[s]\n", (double)(end - start) / CLOCKS_PER_SEC);
	return 0;
}
