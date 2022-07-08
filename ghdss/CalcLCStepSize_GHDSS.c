#include <math.h>
#include "dup.h"

typedef struct _HA_Complex_
{
    float re;
    float im;
} HA_Complex;
HA_Complex  A[ARY_NUM_MIC][ARY_NUM_SRC];        // 演算用配列
HA_Complex  W[ARY_NUM_SRC][ARY_NUM_MIC];        // 演算用配列（分離行列）

HA_Complex  J_LC[ARY_NUM_SRC][ARY_NUM_MIC];     // μLC

void	CalcLCStepSize_GHDSS()
{
	//double fJlcp;
	int num_src=ARY_NUM_SRC;
	int num_mic=ARY_NUM_MIC;
	float fJlcp;
	//double foptLCMyu;
	float foptLCMyu;
	int row, col;
	int	s,m,k;
	//double dTemp ;
	float dTemp ;

	HA_Complex	tmp, tmp2, tmp3;

	////HA_Complex	Atmp[ARY_NUM_MIC][ARY_NUM_SRC]; // TAKIGAHIRA
	////HA_Complex	Zsns[ARY_NUM_SRC][ARY_NUM_SRC]; // TAKIGAHIRA
	HA_Complex	E[ARY_NUM_SRC][ARY_NUM_SRC];
	////HA_Complex	Znsnm[ARY_NUM_SRC][ARY_NUM_MIC]; // TAKIGAHIRA
	////HA_Complex	Znsns[ARY_NUM_SRC][ARY_NUM_SRC]; // TAKIGAHIRA
	HA_Complex	Jlc[ARY_NUM_SRC][ARY_NUM_MIC];

	//////// TAKIGAHIRA : merge formula
	//// 複素共役の値は符号だけの違いであるため、Atmpを使用せず代入による式変形が可能
	//// この関数内でAとしてAtmpの値を使いまわす場所全てで単純置き換えが可能なのでAtmp配列を削減するために統合する
	//// Zsns配列は単位行列の保存のためだけに存在し、単位行列の減算のみで使用している
	//// そのため、 E = E - Znsns （プログラム内ではZnsns配列ではなくZsns配列）は単純化が可能なので統合する
	////
	////for(col=0; col<caA.size2(); col++){
	////	for(row=0; row<caA.size1(); row++){
	////		Atmp(row, col) = complex<double> (caA(row, col).real(), caA(row, col).imag() * (-1.0)) ;
	////	}
	////}
	//OLD-1//for(m=0; m<num_mic; m++){
	//OLD-1//	for(s=0; s<num_src; s++){
	//OLD-1//		Atmp[m][s].re = A[m][s].re;
	//OLD-1//		//Atmp[m][s].im = A[m][s].im * (-1.0) ;	// 複素共役
	//OLD-1//		Atmp[m][s].im = -A[m][s].im ;	// 複素共役
	//OLD-1//		//LOGI("m=%d, s=%d : Atmp[m][s](%f, %f)\n", Atmp[m][s].re, Atmp[m][s].im);
	//OLD-1//	}
	//OLD-1//}
	//// ■ Znsns = eye(Nsrc) ;	単位行列の作成
	////if(_EyeCompMat(&Znsns) == false){
	////	return 0;
	////}
	//OLD-1//for(s=0; s<num_src; s++){
	//OLD-1//	for(m=0; m<num_src; m++){
	//OLD-1//		if(s==m){ Zsns[s][m].re = 1.0; Zsns[s][m].im = 0.0; }
	//OLD-1//		else{     Zsns[s][m].re = 0.0; Zsns[s][m].im = 0.0; }
	//OLD-1//	}
	//OLD-1//}
	//OLD-2//for(s=0; s<num_src; s++){
	//OLD-2//	for(m=0; m<num_src; m++){
	//OLD-2//		Zsns[s][m].re = 0.0;
	//OLD-2//		Zsns[s][m].im = 0.0;
	//OLD-2//	}
	//OLD-2//	Zsns[s][s].re = 1.0;
	//OLD-2//}
	//// ■ E = W*A ;
	////E = prod(caW, Atmp);
	//OLD-1//for(s=0; s<num_src; s++){
	//OLD-1//	for(k=0; k<num_src; k++){
	//OLD-1//		HA_Complex	tmp;
	//OLD-1//		tmp.re	= 0.0;
	//OLD-1//		tmp.im	= 0.0;
	//OLD-1//		for(m=0; m<num_mic; m++){
	//OLD-1//			tmp.re	+= (W[s][m].re * Atmp[m][k].re) - (W[s][m].im * Atmp[m][k].im);
	//OLD-1//			tmp.im 	+= (W[s][m].re * Atmp[m][k].im) + (W[s][m].im * Atmp[m][k].re);
	//OLD-1//		}
	//OLD-1//		E[s][k].re = tmp.re;
	//OLD-1//		E[s][k].im = tmp.im;
	//OLD-1//	}
	//OLD-1//}
	//OLD-2//for(s=0; s<num_src; s++){
	//OLD-2//	for(k=0; k<num_src; k++){
	//OLD-2//		HA_Complex	tmp;
	//OLD-2//		tmp.re	= 0.0;
	//OLD-2//		tmp.im	= 0.0;
	//OLD-2//		for(m=0; m<num_mic; m++){
	//OLD-2//			tmp.re	+= (W[s][m].re * A[m][k].re) + (W[s][m].im * A[m][k].im);
	//OLD-2//			tmp.im 	+= (W[s][m].im * A[m][k].re) - (W[s][m].re * A[m][k].im);
	//OLD-2//		}
	//OLD-2//		E[s][k].re = tmp.re;
	//OLD-2//		E[s][k].im = tmp.im;
	//OLD-2//	}
	//OLD-2//}
	//for(s=0; s<num_src; s++){
	//	for(m=0; m<num_mic; m++){
	//		LOGI("s=%d, m=%d : W[s][m](%f, %f)\n", s, m, W[s][m].re, W[s][m].im);
	//	}
	//}
	//// ■ E = E -Znsns ;
	////E = E - Znsns;
	//OLD-1//for(s=0; s<num_src; s++){
	//OLD-1//	for(m=0; m<num_src; m++){
	//OLD-1//		E[s][m].re -= Zsns[s][m].re;
	//OLD-1//		E[s][m].im -= Zsns[s][m].im;
	//OLD-1//		LOGI("s=%d, m=%d : E[s][m](%f, %f)\n", s, m, E[s][m].re, E[s][m].im);
	//OLD-1//	}
	//OLD-1//}
	//OLD-2//for(s=0; s<num_src; s++){
	//OLD-2//	E[s][s].re -= Zsns[s][s].re;
	//OLD-2//	E[s][s].im -= Zsns[s][s].im;
	//OLD-2//}
	//OLD-3//for(s=0; s<num_src; s++){
	//OLD-3//	E[s][s].re -= 1.0f;
	//OLD-3//	E[s][s].im -= 0.0f;
	//OLD-3//}
	////-- merged formula --//// ■ E = W*A - eye(Nsrc);
	for(s=0; s<num_src; s++){
		for(k=0; k<num_src; k++){
			tmp.re = 0.0f;
			tmp.im = 0.0f;
			for(m=0; m<num_mic; m++){
				tmp2 = W[s][m];
				tmp3 = A[m][k];
				tmp.re	+= tmp2.re * tmp3.re;
				tmp.re	+= tmp2.im * tmp3.im;
				tmp.im 	+= tmp2.im * tmp3.re;
				tmp.im 	-= tmp2.re * tmp3.im;
			}
			E[s][k].re = tmp.re; 
			E[s][k].im = tmp.im; 
		}
		E[s][s].re -= 1.0f;
		//E[s][s].im -= 0.0f;
	}
	//for(s=0; s<num_src; s++){
	//	for(m=0; m<num_src; m++){
	//		//LOGI("s=%d, m=%d : E[s][m](%f, %f)\n", s, m, E[s][m].re, E[s][m].im);
	//	}
	//}

	// LC_Const = FULL は以上まで
	// LC_Const = DIAG
#ifdef ENABLE_SWITCHING_GHDSS_MODE
	//以下、FULLなのでとりあえずスキップ
	switch(){
	case LC_CONST_DIAG:
	{
		////if(iLCConst == _LC_CONST_DIAG){
		////	if(_DiagCompMat(&E, E) == false){	// Eを対角行列に。
		////		return 0;
		////	}
		////} //end LC_CONST = DIAG
		for(s=0; s<num_src; s++){ // diag(E)
			tmp.re = E[s][s].re; // 対角成分E[s][s]のバックアップ
			tmp.im = E[s][s].im;
			for(m=0; m<num_src; m++){
				E[s][m].re = 0.0f; // Eのクリア
				E[s][m].im = 0.0f;
			}
			E[s][s].re = tmp.re; // 対角成分E[s][s]のリストア
			E[s][s].im = tmp.im;
		}
	}
		break;
	case LC_CONST_FULL: // default is FULL
	default:
		break;
	}
#endif

	//////// TAKIGAHIRA : merge formula
	//// 複素共役転置の値は符号と配列のxy交換だけの違いであるため、Znsnmを使用せず代入による式変形が可能
	//// 既にAtmpの時点で複素共役をとっているので複素共役の複素共役となり行列に対する転置のみ行う
	//// この関数内でA'としてZnsnmの値を使いまわす事がなかったのでZnsnm配列を削減するためにJlcの算出ループと統合する
	////
	// ■ Znsnm = A' ;  共役転置複素数
	//cdMat Znsnm (iNsrc, iNmic);
	//Znsnm = trans(conj(Atmp));
	//OLD-1//for(s=0; s<num_src; s++){
	//OLD-1//	for(m=0; m<num_mic; m++){
	//OLD-1//		Znsnm[s][m].re =  Atmp[m][s].re;
	//OLD-1//		Znsnm[s][m].im = -Atmp[m][s].im;	// 複素共役
	//OLD-1//	}
	//OLD-1//}
	//OLD-2//for(s=0; s<num_src; s++){
	//OLD-2//	for(m=0; m<num_mic; m++){ // TAKIGAHIRA : [注]転置行列のため[m][s]を[s][m]へ代入
	//OLD-2//		Znsnm[s][m].re = A[m][s].re;
	//OLD-2//		Znsnm[s][m].im = A[m][s].im;
	//OLD-2//		//LOGI("s=%d, m=%d : Znsnm[s][m](%f, %f)\n", s, m, Znsnm[s][m].re, Znsnm[s][m].im);
	//OLD-2//	}
	//OLD-2//}
	//// ■ Jlc = E * Znsnm ;
	////cdMat Jlc (iNsrc, iNmic);
	////Jlc = prod(E, Znsnm);
	//OLD-1//for(s=0; s<num_src; s++){
	//OLD-1//	for(m=0; m<num_mic; m++){
	//OLD-1//		HA_Complex	tmp;
	//OLD-1//		tmp.re	= 0.0;
	//OLD-1//		tmp.im	= 0.0;
	//OLD-1//		for(k=0; k<num_src; k++){
	//OLD-1//			tmp.re += (E[s][k].re * Znsnm[k][m].re) - (E[s][k].im * Znsnm[k][m].im);
	//OLD-1//			tmp.im += (E[s][k].re * Znsnm[k][m].im) + (E[s][k].im * Znsnm[k][m].re);
	//OLD-1//		}
	//OLD-1//		Jlc[s][m].re = tmp.re;
	//OLD-1//		Jlc[s][m].im = tmp.im;
	//OLD-1//	}
	//OLD-1//}
	////-- merged formula --//// Jlc = E * A';
	for(s=0; s<num_src; s++){
		for(m=0; m<num_mic; m++){
			tmp.re = 0.0f;
			tmp.im = 0.0f;
			for(k=0; k<num_src; k++){
				tmp2 = E[s][k];
				tmp3 = A[m][k];
				tmp.re += tmp2.re * tmp3.re;
				tmp.re -= tmp2.im * tmp3.im;
				tmp.im += tmp2.re * tmp3.im;
				tmp.im += tmp2.im * tmp3.re;
			}
			Jlc[s][m].re = tmp.re;
			Jlc[s][m].im = tmp.im;
		}
	}
	//for(s=0; s<num_src; s++){
	//	for(m=0; m<num_mic; m++){
	//		LOGI("s=%d, m=%d : Jlc[s][m](%f, %f)\n", s, m, Jlc[s][m].re, Jlc[s][m].im);
	//	}
	//}

	// --- ステップサイズと更新量の決定 ---
#ifdef ENABLE_SWITCHING_GHDSS_MODE
	switch(iLCMethod){
	case LC_METHOD_FIX:
	{
		////for(row = 0 ; row < Jlc.size1() ; row++){		//rowでループ
		////	for(col = 0 ; col < Jlc.size2() ;col++){
		////		if(Jlc(row, col).real() * (double)fLCMyu > __MAX_FLOAT_VALUE
		////					&&  Jlc(row, col).imag() * (double)fLCMyu  > __MAX_FLOAT_VALUE){
		////			printf("Data is Inf.\n") ;
		////			return 0 ;
		////		}
		////		complex<double> caTemp((HA_Float)(Jlc(row, col).real() * (double)fLCMyu )
		////							, (HA_Float) (Jlc(row, col).imag() * (double)fLCMyu) );
		////		pcamyu_J_LC->insert_element(row, col, caTemp);
		////	}
		////}
		// need input fLCMyu. ex) fLCMyu = 0.001f;
		for(s=0; s<num_src; s++){
			for(m=0; m<num_mic; m++){
				tmp.re = Jlc[s][m].re * fLCMyu;
				tmp.im = Jlc[s][m].im * fLCMyu;
				/*
				if(tmp.re > _MAX_FLOAT_VALUE && tmp.im > _MAX_FLOAT_VALUE){
					LOGI("Data is Inf.\n") ;
					return 0;
				} */

				J_LC[s][m].re = tmp.re;
				J_LC[s][m].im = tmp.im;
			}
		}
	}
		break;
	case LC_METHOD_ADAPTIVE: // default is ADAPTIVE
	default:
#endif
	{ // <- とりあえずこっち
		//////// TAKIGAHIRA : merge formula
		//// ABS計算で虚数部が0.0fとなり、またABS計算のルート後にfJlcpの計算で実数部を二乗する事が分かっているので無駄な計算を省く
		//// 更に、fJlcpのSum計算以外でZnsnm配列自体が不要となる事からメモリ削減のためSumのループと統合する、虚数部は計算しない
		////
		////if(_AbsCompMat(&Znsnm, Jlc) == false){
		////	return 0;
		////}
		//OLD-2//for(s=0; s<num_src; s++){
		//OLD-2//	for(m=0; m<num_mic; m++){
		//OLD-2//		//dTemp = pow((double)Jlc[s][m].re, 2.0) + pow((double)Jlc[s][m].im, 2.0);
		//OLD-2//		dTemp = (Jlc[s][m].re * Jlc[s][m].re) + (Jlc[s][m].im * Jlc[s][m].im);
		//OLD-2//		Znsnm[s][m].re = sqrtf(dTemp); // change function from sqrt() to sqrtf()
		//OLD-2//		Znsnm[s][m].im = 0.0;
		//OLD-2//	}
		//OLD-2//}
		////fJlcp = 0.0;
		////for(row=0; row<Znsnm.size1(); row++){
		////	for(col=0; col<Znsnm.size2(); col++){
		////		fJlcp = (double)(pow(Znsnm(row, col).real(), (double)2.0) + fJlcp);
		////	}
		////}
		//OLD-1//fJlcp = 0.0f;
		//OLD-1//for(s=0; s<num_src; s++){
		//OLD-1//	for(m=0; m<num_mic; m++){
		//OLD-1//		//fJlcp = (double)(pow(Znsnm[s][m].re, (double)2.0) + fJlcp);
		//OLD-1//		//fJlcp = ((Znsnm[s][m].re * Znsnm[s][m].re) + fJlcp);
		//OLD-1//		fJlcp += (Znsnm[s][m].re * Znsnm[s][m].re);
		//OLD-1//	}
		//OLD-1//}
		//OLD-2//fJlcp = 0.0f;
		//OLD-2//for(s=0; s<num_src; s++){
		//OLD-2//	for(m=0; m<num_mic; m++){
		//OLD-2//		//fJlcp = (double)(pow(Znsnm[s][m].re, (double)2.0) + fJlcp);
		//OLD-2//		//fJlcp = ((Znsnm[s][m].re * Znsnm[s][m].re) + fJlcp);
		//OLD-2//		fJlcp += (Znsnm[s][m].re * Znsnm[s][m].re);
		//OLD-2//	}
		//OLD-2//}
		////fJlcp = (double)(fJlcp + 0.0000000001);
		//OLD-1//fJlcp = (fJlcp + 0.0000000001);
		////-- merged formula --//// fJlcp = sum(abs(Jlc) .* abs(Jlc))  ==> fJlcp = ((fJlcp + avoid 0div value) * 2)
		fJlcp = 0.0f;
		for(s=0; s<num_src; s++){
		#pragma HLS pipeline ii=1

			for(m=0; m<num_mic; m++){
				tmp2 = Jlc[s][m];
				tmp.re = tmp2.re * tmp2.re;
				tmp.im = tmp2.im * tmp2.im;
				fJlcp += (tmp.re + tmp.im);
			}
		}
		fJlcp += 0.0000000001f;
		/*
		if(fabsf(fJlcp) < _MINIMUM_ZERO_VALUE){
			fJlcp += _MINIMUM_ZERO_VALUE;
		} */
		fJlcp *= 2;

		//////// TAKIGAHIRA : merge formula
		//// ABS計算で虚数部が0.0fとなり、またABS計算のルート後にfJlcpの計算で実数部を二乗する事が分かっているので無駄な計算を省く
		//// 更に、foptLCMyuのSum計算以外でZnsns配列自体が不要となる事からメモリ削減のためSumのループと統合する、虚数部は計算しない
		////
		////if(_AbsCompMat(&Znsns, E) == false){
		////	return 0;
		////}
		//OLD-1//for(s=0; s<num_src; s++){
		//OLD-1//	for(m=0; m<num_src; m++){
		//OLD-1//		//dTemp = pow((double)E[s][m].re, 2.0) + pow((double)E[s][m].im, 2.0);
		//OLD-1//		dTemp = (E[s][m].re * E[s][m].re) + (E[s][m].im * E[s][m].im);
		//OLD-1//		Znsns[s][m].re = sqrtf(dTemp); // change function from sqrt() to sqrtf()
		//OLD-1//		Znsns[s][m].im = 0.0;
		//OLD-1//		//LOGI("s=%d, m=%d : Znsns[s][m](%f, %f), E[s][m](%f, %f)\n", Znsns[s][m].re, Znsns[s][m].im, E[s][m].re, E[s][m].im);
		//OLD-1//	}
		//OLD-1//}
		////foptLCMyu = 0.0;
		////for(row=0; row<Znsns.size1(); row++){
		////	for(col=0; col<Znsns.size2(); col++){
		////		foptLCMyu = (double) (pow(Znsns(row, col).real(), (double)2.0) + foptLCMyu );
		////	}
		/////}
		//OLD-1//foptLCMyu = 0.0f;
		//OLD-1//for(s=0; s<num_src; s++){
		//OLD-1//	for(m=0; m<num_src; m++){
		//OLD-1//		//foptLCMyu = (double) (pow(Znsns[s][m].re, (double)2.0) + foptLCMyu );
		//OLD-1//		//foptLCMyu = ((Znsns[s][m].re * Znsns[s][m].re) + foptLCMyu );
		//OLD-1//		foptLCMyu += (Znsns[s][m].re * Znsns[s][m].re);
		//OLD-1//	}
		//OLD-1//}
		//// fJlcpに対する0除算防止値の加算と予め2倍にしておくことでfoptLCMyu計算の除算を減らせることから計算位置を移動
		//OLD-1//if(fabsf(fJlcp) < _MINIMUM_ZERO_VALUE){ // change function from fabs() to fabsf()
		//OLD-1//	foptLCMyu = foptLCMyu / (fJlcp+_MINIMUM_ZERO_VALUE) / 2 ;
		//OLD-1//}
		//OLD-1//else{
		//OLD-1//	foptLCMyu = foptLCMyu / fJlcp / 2 ;
		//OLD-1//}
		//OLD-2//if(fabsf(fJlcp) < _MINIMUM_ZERO_VALUE){
		//OLD-2//	fJlcp += _MINIMUM_ZERO_VALUE;
		//OLD-2//}
		//OLD-2//fJlcp *= 2;
		//OLD-2//foptLCMyu = foptLCMyu / fJlcp;
		////-- merged formula --//// foptLCMyu = sum(abs(E) .* abs(E)) / fJlcp  [ただし、fJlcp=((fJlcp+avoid_0div_value)*2)済み]
		foptLCMyu = 0.0f;
	#pragma HLS pipeline ii=1

		for(s=0; s<num_src; s++){
			for(m=0; m<num_src; m++){
				tmp2 = E[s][m];
				tmp.re = tmp2.re * tmp2.re;
				tmp.im = tmp2.im * tmp2.im;
				foptLCMyu += (tmp.re + tmp.im);
			}
		}
		foptLCMyu = foptLCMyu / fJlcp;

		////for(row=0; row<Jlc.size1(); row++){
		////	for(col = 0 ; col <  Jlc.size2() ; col++ ){
		////		if(Jlc(row, col).real() * foptLCMyu > __MAX_FLOAT_VALUE && Jlc(row, col).imag() * foptLCMyu > __MAX_FLOAT_VALUE){
		////			printf("Data is Inf.\n") ;
		////			return 0;
		////		}
		////		complex<double> caTemp((HA_Float)( Jlc(row, col).real() * foptLCMyu ),(HA_Float) (Jlc(row, col).imag() * foptLCMyu) );
		////		pcamyu_J_LC->insert_element(row, col, caTemp);
		////	}
		////}
		for(s=0; s<num_src; s++){
		#pragma HLS pipeline ii=1

			for(m=0; m<num_mic; m++){
				tmp.re = Jlc[s][m].re * foptLCMyu;
				tmp.im = Jlc[s][m].im * foptLCMyu;
//				if(tmp.re > _MAX_FLOAT_VALUE && tmp.im > _MAX_FLOAT_VALUE){
//					LOGI("Data is Inf.\n") ;
//					return ;
	//			}
				J_LC[s][m].re = tmp.re;
				J_LC[s][m].im = tmp.im;
			}
		}
	} // end LC_Method = Adaptive
#ifdef ENABLE_SWITCHING_GHDSS_MODE
		break;
	}
#endif
	//for(s=0; s<num_src; s++){
	//	for(m=0; m<num_mic; m++){
	//		LOGI("s=%d, m=%d : J_LC[s][m](%f, %f)\n", s, m, J_LC[s][m].re, J_LC[s][m].im);
	//	}
	//}
	return;
}
