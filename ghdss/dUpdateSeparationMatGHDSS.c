// #define THRESHOLD (1000)
#define ha_iUpdate  (0)     // SS_METHOD
#define ha_fWmyu    (0)     //
#include <math.h>
#include "dup.h"
typedef struct _HA_Complex_
{
    float re;
    float im;
} HA_Complex;

void CalcLCStepSize_GHDSS();
void CalcSSStepSize_GHDSS();

HA_Complex  J_SS[ARY_NUM_SRC][ARY_NUM_MIC];     // μSS
HA_Complex  J_LC[ARY_NUM_SRC][ARY_NUM_MIC];     // μLC
HA_Complex  A[ARY_NUM_MIC][ARY_NUM_SRC];        // 演算用配列
HA_Complex  W[ARY_NUM_SRC][ARY_NUM_MIC];        // 演算用配列（分離行列）
HA_Complex  X[ARY_NUM_SRC];
HA_Complex  Y[ARY_NUM_SRC];

void	dUpdateSeparationMatGHDSS(
//	HA_Complex inputSpec[ARY_NUM_MIC][ARY_NUM_FFT],
	float inputSpec_re[ARY_NUM_MIC*ARY_NUM_FFT],
	float inputSpec_im[ARY_NUM_MIC*ARY_NUM_FFT],
//	HA_Complex      OutTF[ARY_NUM_SRC][ARY_NUM_FFT][ARY_NUM_MIC]
	float OutTF_re[ARY_NUM_SRC*ARY_NUM_FFT*ARY_NUM_MIC],
	float OutTF_im[ARY_NUM_SRC*ARY_NUM_FFT*ARY_NUM_MIC],
	float outsepMat_re[ARY_NUM_SRC*ARY_NUM_FFT*ARY_NUM_MIC],
    float outsepMat_im[ARY_NUM_SRC*ARY_NUM_FFT*ARY_NUM_MIC],
	float THRESHOLD[1], int ifreqstart[1], int ifreqnum[1],
	char startt[1], char stopt[1]
	)
{
#pragma HLS INTERFACE s_axilite port=inputSpec_re
#pragma HLS INTERFACE s_axilite port=inputSpec_im
#pragma HLS INTERFACE s_axilite port=OutTF_re
#pragma HLS INTERFACE s_axilite port=OutTF_im
#pragma HLS INTERFACE s_axilite port=outsepMat_re
#pragma HLS INTERFACE s_axilite port=outsepMat_im
#pragma HLS INTERFACE s_axilite port=THRESHOLD
#pragma HLS INTERFACE s_axilite port=ifreqstart
#pragma HLS INTERFACE s_axilite port=ifreqnum
#pragma HLS INTERFACE axis port=startt
#pragma HLS INTERFACE axis port=stopt
#pragma HLS INTERFACE s_axilite port=return
#pragma HLS array_partition variable=W dim=2 type=cyclic factor=4
#pragma HLS array_partition variable=X type=cyclic factor=4

int num_src=ARY_NUM_SRC;
int num_mic=ARY_NUM_MIC;
//    int 		Nsrc, Nfreq, Nmic, ifreq, ind;
    int 		ifreq, ind;
	float 	Threshold_val, Max_val, fTmp;
	HA_Complex tmp, tmp2, tmp3;

	HA_Complex  sepMat[ARY_NUM_SRC][ARY_NUM_FFT][ARY_NUM_MIC];                 // 演算用配列
	int		m,f,s;	
	int 	iUp = 0;

	int             mat_id[ARY_NUM_SRC];
	HA_Complex      csaOutputSpec[ARY_NUM_SRC][ARY_NUM_FFT];


	// Dummy: set sepMat
	for(s=0;s<ARY_NUM_SRC;s++) 
		for(f=0;f<ARY_NUM_FFT;f++)
			for(m=0;m<ARY_NUM_MIC;m++) {
				sepMat[s][f][m].re = 1.0;
				sepMat[s][f][m].im = 2.0;
			}
			startt[0]=(char)sepMat[9][ARY_NUM_FFT-1][ARY_NUM_MIC-1].im;
	// Dummy: set matid
	for(s=0;s<ARY_NUM_SRC;s++) 
		mat_id[s] = s;
/*
	for(s=0;s<10;s++) 
			for(m=0;m<ARY_NUM_MIC;m++) {
				J_SS[s][m].re=0.001;
				J_SS[s][m].im=0.002;
				J_LC[s][m].re=0.001;
				J_LC[s][m].im=0.002;
			}
*/

	//Treshold_val

	Max_val = 32768.0f;
	Max_val = 10 * log10(Max_val);
	fTmp 	= (Max_val - THRESHOLD[0]) / 10.0f;
	//TAKIGAHIRA
	//Threshold_val = pow((float)10.0 , fTmp) ;
	{
		union {
			double d_val;
			int i_val[2];
		} u;
		u.i_val[1] = (int)(fTmp*3468673+1072632447);
		u.i_val[0] = 0;
		Threshold_val = (float) u.d_val;
	}

	//■　初期化＆チェック
	//Nmic 	= ARY_NUM_MIC;
	//Nfreq 	= ARY_NUM_FFT;
	//Nsrc 	= ARY_NUM_SRC;

	/////////////////////////////
	// メインループ
mloop:	for (ifreq=ifreqstart[0]; ifreq<ifreqstart[0]+ifreqnum[0];ifreq++) {
//	#pragma HLS unroll facter=4
		//printf("--------- ifreq: %d --------\n", ifreq);
		// 範囲チェック。とりあえずパス。

		// ■ カレントの周波数帯のデータのみを切り出して，Matrix型に変換
		//cdMat A (Nmic, Nsrc);	// 8 x 3
		//cdMat X (Nmic, 1);	// 8 x 1
		//cdMat W (Nsrc, Nmic);	// 3 x 8

		// psSavedSeparation->micary_tdTFDB を変換
		// _Convert_ComplexSpec() 相当
		for(m=0; m<num_mic; m++){
			for(s=0; s<num_src; s++){
				A[m][s].re = OutTF_re[m*num_mic*num_src+num_mic*ifreq+s];
				A[m][s].im = OutTF_im[m*num_mic*num_src+num_mic*ifreq+s];
				//LOGI("m=%d, s=%d : A[m][s](%15.8f, %15.8f)\n", A[m][s].re, A[m][s].im);
				//LOGI("m=%d, s=%d : A[m][s](%f, %f)\n", A[m][s].re, A[m][s].im);
			}
		}
		//LOGI("A[0][0](%f, %f)\n", A[0][0].re, A[0][0].im);

		// csInputSpecを変換
		// _Convert_Complex_InputSpec()相当
		for(m=0; m<num_mic; m++){
			X[m].re = inputSpec_re[m+ifreq*ARY_NUM_MIC];
			X[m].im = inputSpec_im[m+ifreq*ARY_NUM_MIC];
		}
		//LOGI("X[0](%f, %f)\n", X[0].re, X[0].im);

		// psSavedSeparation->micary_smSeparation を変換
		// _Convert_ComplexSpec()相当
		for(s=0; s<num_src; s++){
			for(m=0; m<num_mic; m++){
				// KAMA
				//W[s][m].re = sepMat[s][ifreq][m].re;
				//W[s][m].im = sepMat[s][ifreq][m].im;
				W[s][m].re = sepMat[mat_id[s]][ifreq][m].re;
				W[s][m].im = sepMat[mat_id[s]][ifreq][m].im;
				//LOGI("s=%d, m=%d : W[s][m](%f, %f)\n", s, m, W[s][m].re, W[s][m].im);
			}
		}

		// ---------- Execute GHDSS --------------- //
		//	Y = prod(W, X); // B = prod(A) は、配列の次元ごとに積を返します
		//
		//	(A+jB)(C+jD)
		// = (AC-BD)+j(AD+BC)
loop10:		for(s=0; s<num_src; s++){
			tmp.re = 0.0f;
			tmp.im = 0.0f;
	loop11:		for(m=0; m<num_mic; m++){
	//#pragma HLS unroll factor=4
#pragma HLS pipeline
				//float	aa = W[s][m].re;
				//float	bb = W[s][m].im;
				//float	cc = X[m].re;
				//float	dd = X[m].im;

				//Y[s].re += (aa*cc - bb*dd);
				//Y[s].im += (aa*dd + bb*cc);
				tmp2 = W[s][m];
				tmp3 = X[m];
				tmp.re += tmp2.re * tmp3.re;
				tmp.re -= tmp2.im * tmp3.im;
				tmp.im += tmp2.re * tmp3.im;
				tmp.im += tmp2.im * tmp3.re;
				Y[s] = tmp;
				//if(ifreq == 256){ LOGI("%d : %f, %f\n", ifreq, X[m].re, X[m].im); }
				//if(ifreq == 256){ LOGI("%d : %f, %f\n", ifreq, W[s][m].re, W[s][m].im); }
			}
		}

	    // --- Copy : Y -> Separated spectrum(src->freqの順) --- //
		//csSeparationSpecをYで更新（出力用）

		for(s=0; s<num_src; s++){
			csaOutputSpec[s][ifreq].re = Y[s].re;
			csaOutputSpec[s][ifreq].im = Y[s].im;
		}

		// --------------------------------------- //
		// ------ Update separated matrix --------- //
		for(m=0; m<num_mic; m++){
			//fTmp = pow(X[m].re , (double)2.0) + pow(X[m].im , (double)2.0) ;
			fTmp  = (X[m].re * X[m].re); // expand pow function
			fTmp += (X[m].im * X[m].im); // expand pow function
			fTmp = sqrtf(fTmp); // change function from sqrt() to sqrtf()

			if(fTmp < Threshold_val){
				iUp = 1 ;
				break ;
			}
		}
		//	printf("iUp = %d\n", iUp);
		if(iUp == 0){
			/////////////////////
			CalcSSStepSize_GHDSS();
			if(ha_iUpdate == 0){	//段階的
				// Update W
//				W = W - J_SS;
				for(s=0; s<num_src; s++){
					for(m=0; m<num_mic; m++){
						W[s][m].re	-= J_SS[s][m].re;
						W[s][m].im	-= J_SS[s][m].im;
					}
				}
			}

			/////////////////////
			// ■　LC
			CalcLCStepSize_GHDSS();
		    if(ha_iUpdate == 0){	// 段階的
				// Update W
//				W = W - J_LC;
				for(s=0; s<num_src; s++){
					for(m=0; m<num_mic; m++){
						W[s][m].re	-= J_LC[s][m].re;
						W[s][m].im	-= J_LC[s][m].im;
						W[s][m].re	*= (1.0 - ha_fWmyu);
						W[s][m].im	*= (1.0 - ha_fWmyu);
					}
				}
		    }
		    else{
				//W = W - J_SS;
				//W = W - J_LC;
				for(s=0; s<num_src; s++){
					for(m=0; m<num_mic; m++){
						W[s][m].re	*= (1.0 - ha_fWmyu);
						W[s][m].im	*= (1.0 - ha_fWmyu);
						W[s][m].re	-= J_LC[s][m].re;
						W[s][m].im	-= J_LC[s][m].im;
					}
				}
			}
		}
		stopt[0]=(char) W[s-1][m-1].im;
		for(s=0; s<num_src; s++){
			for(m=0; m<num_mic; m++){
				// KAMA
				//sepMat[s][ifreq][m].re = W[s][m].re;
				//sepMat[s][ifreq][m].im = W[s][m].im;
				sepMat[mat_id[s]][ifreq][m].re = W[s][m].re;
				sepMat[mat_id[s]][ifreq][m].im = W[s][m].im;
			}
		}
	}	//end while
	for(s=0; s<ARY_NUM_SRC; s++)
        for(f=0; f<ARY_NUM_FFT; f++)
            for(m=0; m<ARY_NUM_MIC; m++){
                outsepMat_re[s*ARY_NUM_FFT*ARY_NUM_MIC+f*ARY_NUM_MIC+m]=sepMat[s][f][m].re;
                outsepMat_im[s*ARY_NUM_FFT*ARY_NUM_MIC+f*ARY_NUM_MIC+m]=sepMat[s][f][m].im;
            }
//	return;
}
