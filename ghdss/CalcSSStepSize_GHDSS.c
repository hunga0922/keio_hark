#include <math.h>
#include "dup.h"

typedef struct _HA_Complex_
{
    float re;
    float im;
} HA_Complex;
HA_Complex  X[ARY_NUM_SRC]; 
HA_Complex  Y[ARY_NUM_SRC]; 
HA_Complex  J_SS[ARY_NUM_SRC][ARY_NUM_MIC]; 

void	CalcSSStepSize_GHDSS()
{
	int s, m;

	HA_Complex	E[ARY_NUM_SRC][ARY_NUM_MIC];

	HA_Complex	Jss[ARY_NUM_SRC][ARY_NUM_MIC];
	HA_Complex	Jpri[ARY_NUM_SRC][ARY_NUM_MIC];

    double fTmp,fYrms;
    double fTmp2 ;
	int num_mic=ARY_NUM_MIC;
	int num_src=ARY_NUM_SRC;

	HA_Complex	tmp, tmp2, tmp3;
	HA_Complex  PY[ARY_NUM_SRC]; 
	HA_Complex  QY[ARY_NUM_SRC]; 

	// GetPQY	

	float 	fabsY, fabsYs, ftmp ;
	HA_Complex 	cmpExpY ;	// <- これはFloatで演算している

	for(s=0; s<num_src; s++){
	#pragma HLS pipeline II=1
//		fabsY 	= sqrtf((Y[s].re * Y[s].re) + (Y[s].im * Y[s].im)) + _MINIMUM_ZERO_VALUE; // expand pow function, change function from sqrt() to sqrtf()
		fabsY 	= sqrtf((Y[s].re * Y[s].re) + (Y[s].im * Y[s].im)) ; // expand pow function, change function from sqrt() to sqrtf()
		//fabsYs 	= fabsY * 1.0;			// ha_fSSScal = 1.0
		fabsYs 	= fabsY;			// ha_fSSScal = 1.0

		cmpExpY.re = Y[s].re / fabsY ;
		cmpExpY.im = Y[s].im / fabsY ;
	    
		PY[s].re = tanhf(fabsYs) * cmpExpY.re; // change function from tanh() to tanhf()
		PY[s].im = tanhf(fabsYs) * cmpExpY.im; // change function from tanh() to tanhf()

		ftmp = coshf(fabsYs); // change function from cosh() to coshf()
		ftmp = ftmp * ftmp ;
		ftmp = fabsYs / ftmp;		// ha_fSSScal = 1.0

		QY[s].re = (tanhf(fabsYs) + ftmp) * cmpExpY.re; // change function from tanh() to tanhf()
		QY[s].im = (tanhf(fabsYs) + ftmp) * cmpExpY.im; // change function from tanh() to tanhf()
	}
	////-- merged formula --//// E = (y*y') - diag(y*y')
	for(s=0; s<num_src; s++){
		tmp2 = PY[s];
		for(m=0; m<num_src; m++){ // 対角を含めRyyつまり(y*y')を代入
			tmp3 = Y[m];
			tmp.re  = tmp2.re * tmp3.re;
			tmp.re += tmp2.im * tmp3.im;
			E[s][m].re = tmp.re;
			tmp.im  = tmp2.im * tmp3.re;
			tmp.im -= tmp2.re * tmp3.im;
			E[s][m].im = tmp.im;
		}
		E[s][s].re = 0.0f;
		E[s][s].im = 0.0f; // 対角だけゼロで上書き
	}

	////-- merged formula --//// Jss = (E * QY) * X'
	for(m=0; m<num_mic; m++){
		//X[m].re =  X[m].re;
		X[m].im = -X[m].im;	// 複素共役
	}
	for(s=0; s<num_src; s++){ // tmp は Zns1[s] に相当
		tmp.re = 0.0f;
		tmp.im = 0.0f;
		for(m=0; m<num_src; m++){
			tmp2 = E[s][m];
			tmp3 = QY[m];
			tmp.re += tmp2.re * tmp3.re;
			tmp.re -= tmp2.im * tmp3.im;
			tmp.im += tmp2.re * tmp3.im;
			tmp.im += tmp2.im * tmp3.re;
		}
		for(m=0; m<num_mic; m++){
			tmp2 = X[m];
			Jss[s][m].re = (tmp.re * tmp2.re) - (tmp.im * tmp2.im);
			Jss[s][m].im = (tmp.re * tmp2.im) + (tmp.im * tmp2.re);
		}
	}

	for(s=0; s<num_src; s++){
	#pragma HLS pipeline II=1

		for(m=0; m<num_mic; m++){
			Jpri[s][m].re = Jss[s][m].re;
			Jpri[s][m].im = -Jss[s][m].im;
		}
	}

	{ // <- とりあえず2なのでこれを先に移植
		fTmp = 0.0f;
		for(s=0; s<num_src; s++){
		#pragma HLS pipeline II=1

			for(m=0; m<num_src; m++){
				tmp = E[s][m];
				fTmp += tmp.re * tmp.re;
				fTmp += tmp.im * tmp.im;
			}
		}

        /* 2014.12.12 ktakagi@sif
        **
        ** 音源IDが異なり＆音源方向が同じ＆発話区間が重複した場合に
        ** 計算のオーバーフローが生じ、分離音が破綻した
        ** 計算をdoubleで実行するように修正
        ** 
        */
        fTmp2 = 0.0;
        {
            double tmp_re;
            double tmp_im;
	    	for(s=0; s<num_src; s++){
		#pragma HLS pipeline II=1

	    		for(m=0; m<num_mic; m++){
	    			//Znsnm[s][m].re = sqrtf((Jss[s][m].re * Jss[s][m].re + Jss[s][m].im * Jss[s][m].im) * (Jpri[s][m].re * Jpri[s][m].re + Jpri[s][m].im * Jpri[s][m].im));
	    			//Znsnm[s][m].im = 0.0f;
	    			//LOGI("s=%d, m=%d : Znsnm(%f, %f)\n", s, m, Znsnm[s][m].re, Znsnm[s][m].im);
	    			tmp2 = Jss[s][m];
	    			tmp_re  = tmp2.re * tmp2.re;
	    			tmp_re += tmp2.im * tmp2.im;
	    			tmp2 = Jpri[s][m];
	    			tmp_im  = tmp2.re * tmp2.re;
	    			tmp_im += tmp2.im * tmp2.im;
	    			fTmp2 += sqrt(tmp_re * tmp_im);
	    		}
	    	}
		}
		/*
		if(fabsf(fTmp2) < 0){ // change function from fabs() to fabsf()
			//fTmp2 = fTmp2 + _MINIMUM_ZERO_VALUE;
			fTmp2 += _MINIMUM_ZERO_VALUE;
		} */

		fTmp2 *= 2.0f;
		fTmp /= fTmp2;			// <- ゼロ割注意！！
		for(s=0; s<num_src; s++){
		#pragma HLS pipeline II=1
			for(m=0; m<num_mic; m++){
				tmp.re = Jss[s][m].re * fTmp;
				tmp.im = Jss[s][m].im * fTmp;
				J_SS[s][m].re = tmp.re;
				J_SS[s][m].im = tmp.im;
				//LOGI("s=%d, m=%d : J_SS[s][m](%f, %f)\n", s, m, J_SS[s][m].re, J_SS[s][m].im);
			}
		}
	} // end iSSMethod = (2 == SS_METHOD_ADAPTIVE)
	fYrms = 0.0f;
	for(s=0; s<num_src; s++){
	#pragma HLS pipeline II=1
		tmp = Y[s];
		fYrms += tmp.re * tmp.re;
		fYrms += tmp.im * tmp.im;
	}

	fYrms /= num_src;
	fYrms = sqrtf(fYrms); // change function from sqrt() to sqrtf()

	if(fYrms < 0.0f){		// <- 後でパラメータ化すること
		////_InitCompMat(pcamyu_J_SS);		//NoiseFloorより小さかったら，更新値をなかったことにする？ ゼロ埋め
		for(s=0; s<num_src; s++){
		#pragma HLS pipeline II=1
			for(m=0; m<num_mic; m++){
				J_SS[s][m].re = 0.0f;
				J_SS[s][m].im = 0.0f;
			}
		}
	}
}

