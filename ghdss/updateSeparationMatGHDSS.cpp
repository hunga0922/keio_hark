#include "GHDSS.hpp"
#include "CalcLCStepSize_GHDSS.h"
#include "CalcSSStepSize_GHDSS.h"

int num_src=ARY_NUM_SRC;
int num_mic=ARY_NUM_MIC;
void dUpdateSeparationMatGHDSS(
//	HA_Complex inputSpec[ARY_NUM_MIC][ARY_NUM_FFT],
	float inputSpec_re[ARY_NUM_MIC*ARY_NUM_FFT],
	float inputSpec_im[ARY_NUM_MIC*ARY_NUM_FFT],
//	HA_Complex      OutTF[ARY_NUM_SRC][ARY_NUM_FFT][ARY_NUM_MIC]
	float OutTF_re[ARY_NUM_SRC*ARY_NUM_FFT*ARY_NUM_MIC],
	float OutTF_im[ARY_NUM_SRC*ARY_NUM_FFT*ARY_NUM_MIC],
	float outsepMat_re[ARY_NUM_SRC*ARY_NUM_FFT*ARY_NUM_MIC],
	float outsepMat_im[ARY_NUM_SRC*ARY_NUM_FFT*ARY_NUM_MIC],
	float THRESHOLD[1], char startt[1], char stopt[1]
	) {
#pragma HLS INTERFACE s_axilite port=inputSpec_re
#pragma HLS INTERFACE s_axilite port=inputSpec_im
#pragma HLS INTERFACE s_axilite port=OutTF_re
#pragma HLS INTERFACE s_axilite port=OutTF_im
#pragma HLS INTERFACE s_axilite port=outsepMat_re
#pragma HLS INTERFACE s_axilite port=outsepMat_im
#pragma HLS INTERFACE s_axilite port=THRESHOLD
#pragma HLS INTERFACE axis port=startt
#pragma HLS INTERFACE axis port=stopt
#pragma HLS INTERFACE s_axilite port=return

//    int 		Nsrc, Nfreq, Nmic, ifreq, ind;
	int ifreq, ind;
	float Threshold_val, Max_val, fTmp;
	HA_Complex tmp, tmp2, tmp3;
	HA_Complex csaOutputSpec[ARY_NUM_SRC][ARY_NUM_FFT];
	HA_Complex sepMat[ARY_NUM_SRC][ARY_NUM_FFT][ARY_NUM_MIC];                 // 演算用配列
	HA_Complex  Y[ARY_NUM_SRC];
	HA_Complex  X[ARY_NUM_SRC];
	HA_Complex  W[ARY_NUM_SRC][ ARY_NUM_MIC];		// 演算用配列（分離行列）
	HA_Complex  J_SS[ARY_NUM_SRC][ARY_NUM_MIC];     // μSS
	HA_Complex  J_LC[ARY_NUM_SRC][ARY_NUM_MIC];     // μLC
	HA_Complex  A[ARY_NUM_MIC][ARY_NUM_SRC];        // 演算用配列

	int m,f,s;
	int iUp = 0;
	int mat_id[ARY_NUM_SRC];

	// Dummy: set sepMat
	for(s=0;s<ARY_NUM_SRC;s++) {
		for(f=0;f<ARY_NUM_FFT;f++) {
			for(m=0;m<ARY_NUM_MIC;m++) {
				sepMat[s][f][m].re = 1.0;
				sepMat[s][f][m].im = 2.0;
			}
		}
	}
	startt[0]=(char)sepMat[9][ARY_NUM_FFT-1][ARY_NUM_MIC-1].im;

	// Dummy: set matid
	for(s=0;s<ARY_NUM_SRC;s++) {
		mat_id[s] = s;
	}

	//Treshold_val

	Max_val = 32768.0f;
	Max_val = 10 * log10(Max_val);
	fTmp 	= (Max_val - THRESHOLD[0]) / 10.0f;

	{
		union {
			double d_val;
			int i_val[2];
		} u;
		u.i_val[1] = (int)(fTmp*3468673+1072632447);
		u.i_val[0] = 0;
		Threshold_val = (float) u.d_val;
	}


	/////////////////////////////
	// メインループ
	for (ifreq=0; ifreq<ARY_NUM_FFT;ifreq++) {
		#pragma HLS PIPELINE
		for(m=0; m<num_mic; m++) {
// Clang-format off
			#pragma HLS loop_tripcount min=8 max=8
			//Clang-format on
			for(s=0; s<num_src; s++) {
// Clang-format off
				#pragma HLS loop_tripcount min=8 max=8
				//Clang-format on
				A[m][s].re = OutTF_re[m*num_mic*num_src+num_mic*ifreq+s];
				A[m][s].im = OutTF_im[m*num_mic*num_src+num_mic*ifreq+s];

			}
		}

		// csInputSpecを変換
		// _Convert_Complex_InputSpec()相当
		for(m=0; m<num_mic; m++){
// Clang-format off
			#pragma HLS loop_tripcount min=8 max=8
			//Clang-format on
			X[m].re = inputSpec_re[m+ifreq*ARY_NUM_MIC];
			X[m].im = inputSpec_im[m+ifreq*ARY_NUM_MIC];
		}
		for(s=0; s<num_src; s++){
// Clang-format off
			#pragma HLS loop_tripcount min=8 max=8
			//Clang-format on
			for(m=0; m<num_mic; m++){
// Clang-format off
				#pragma HLS loop_tripcount min=8 max=8
				//Clang-format on
				W[s][m].re = sepMat[mat_id[s]][ifreq][m].re;
				W[s][m].im = sepMat[mat_id[s]][ifreq][m].im;
			}
		}

		for(s=0; s<ARY_NUM_SRC; s++){
// Clang-format off
			#pragma HLS loop_tripcount min=8 max=8
			//Clang-format on
			tmp.re = 0.0f;
			tmp.im = 0.0f;
			float tmp_re, tmp_im;
	loop11:		for(m=0; m<ARY_NUM_MIC; m++){
// Clang-format off
				#pragma HLS loop_tripcount min=8 max=8
				//Clang-format on
				tmp2 = W[s][m];
				tmp3 = X[m];
				tmp_re += tmp2.re * tmp3.re;
				tmp_re -= tmp2.im * tmp3.im;
				tmp_im += tmp2.re * tmp3.im;
				tmp_im += tmp2.im * tmp3.re;
			}
			tmp.re = tmp_re;
			tmp.im = tmp_im;
			Y[s] = tmp;
		}

	    // --- Copy : Y -> Separated spectrum(src->freqの順) --- //
		//csSeparationSpecをYで更新（出力用）

		for(s=0; s<num_src; s++){
// Clang-format off
			#pragma HLS loop_tripcount min=8 max=8
			//Clang-format on
			csaOutputSpec[s][ifreq].re = Y[s].re;
			csaOutputSpec[s][ifreq].im = Y[s].im;
		}

		// ------ Update separated matrix --------- //
		for(m=0; m<num_mic; m++){
// Clang-format off
			#pragma HLS loop_tripcount min=8 max=8
			#pragma HLS BIND_OP variable=fTmp op=mul impl=dsp
			//Clang-format on
			fTmp  = (X[m].re * X[m].re); // expand pow function
			fTmp += (X[m].im * X[m].im); // expand pow function
			fTmp = sqrtf(fTmp); // change function from sqrt() to sqrtf()

			if(fTmp < Threshold_val){
				iUp = 1 ;
				break ;
			}
		}

		if(iUp == 0){
			/////////////////////
			CalcSSStepSize_GHDSS(X, Y, J_SS);
			if(ha_iUpdate == 0){	//段階的
				// Update W
				for(s=0; s<num_src; s++){
// Clang-format off
					#pragma HLS loop_tripcount min=8 max=8
					//Clang-format on
					for(m=0; m<num_mic; m++){
// Clang-format off
						#pragma HLS loop_tripcount min=8 max=8
						//Clang-format on
						W[s][m].re	-= J_SS[s][m].re;
						W[s][m].im	-= J_SS[s][m].im;
					}
				}
			}

			/////////////////////
			// ■　LC
			CalcLCStepSize_GHDSS(A, W, J_LC);
		    if(ha_iUpdate == 0){	// 段階的
				for(s=0; s<num_src; s++){
// Clang-format off
					#pragma HLS loop_tripcount min=8 max=8
					//Clang-format on
					for(m=0; m<num_mic; m++){
// Clang-format off
						#pragma HLS loop_tripcount min=8 max=8
						//Clang-format on
						W[s][m].re	-= J_LC[s][m].re;
						W[s][m].im	-= J_LC[s][m].im;
						W[s][m].re	*= (1.0 - ha_fWmyu);
						W[s][m].im	*= (1.0 - ha_fWmyu);
					}
				}
		    }
		    else{
				for(s=0; s<num_src; s++){
// Clang-format off
					#pragma HLS loop_tripcount min=8 max=8
					//Clang-format on
					for(m=0; m<num_mic; m++){
// Clang-format off
						#pragma HLS loop_tripcount min=8 max=8
						// Clang-format on
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
// Clang-format off
			#pragma HLS loop_tripcount min=8 max=8
			//Clang-format on
			for(m=0; m<num_mic; m++){
// Clang-format off
				#pragma HLS loop_tripcount min=8 max=8
				//Clang-format on
				// KAMA
				sepMat[mat_id[s]][ifreq][m].re = W[s][m].re;
				sepMat[mat_id[s]][ifreq][m].im = W[s][m].im;
			}
		}
	}	//end while
	for(s=0; s<ARY_NUM_SRC; s++) {
        for(f=0; f<ARY_NUM_FFT; f++) {
            for(m=0; m<ARY_NUM_MIC; m++) {
                outsepMat_re[s*ARY_NUM_FFT*ARY_NUM_MIC+f*ARY_NUM_MIC+m]=sepMat[s][f][m].re;
                outsepMat_im[s*ARY_NUM_FFT*ARY_NUM_MIC+f*ARY_NUM_MIC+m]=sepMat[s][f][m].im;
            }
        }
	}
}
