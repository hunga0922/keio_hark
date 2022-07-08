#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define MMAX 16
#define NMAX 16
typedef struct _COMPLEX_FLOAT_
{
    float real;
    float imag;
} COMPLEX_FLOAT;


void Evd(float a_real[MMAX*NMAX], float a_imag[MMAX*NMAX], float d[NMAX],
	float v_real[MMAX*NMAX], float v_imag[MMAX*NMAX],
	char startt[1], char stopt[1])
{
#pragma HLS INTERFACE s_axilite port=a_real
#pragma HLS INTERFACE s_axilite port=a_imag
#pragma HLS INTERFACE s_axilite port=d
#pragma HLS INTERFACE s_axilite port=v_real
#pragma HLS INTERFACE s_axilite port=v_imag
#pragma HLS INTERFACE axis port=startt
#pragma HLS INTERFACE axis port=stopt
#pragma HLS INTERFACE s_axilite port=return


    static float b[16];
    static float z[16];
    int n;
    float sm;
    float smp[16];
    //double sm;
//#pragma HLS PIPELINE
    float thresh;
    float g;
    //double g;
    COMPLEX_FLOAT h;
    COMPLEX_FLOAT t;
    float c;
    //double c;
    COMPLEX_FLOAT v[MMAX][NMAX];
    COMPLEX_FLOAT s;
    COMPLEX_FLOAT tau;
    //COMPLEX_FLOAT theta;
    
    float norm_a;
    float t_tmp;
    float norm_t;
    //double norm_a;
    //double t_tmp;
    //double norm_t;
    
    int i,j;
    int ip, iq;
    
    float d_tmp;
    float v_tmp;
    float g_real, g_imag;
    float h_real, h_imag;
    float tmp_real, tmp_imag;
    float tmp;
	COMPLEX_FLOAT a[MMAX][NMAX];

	n=16;
    float tmp1[n],tmp2[n];
//#pragma HLS ARRAY_PARTITION variable=tmp1 complete dim=0
//#pragma HLS ARRAY_PARTITION variable=tmp2 complete dim=0
#pragma HLS ARRAY_PARTITION variable=a complete dim=2
//#pragma HLS ARRAY_PARTITION variable=smp complete dim=1
//#pragma HLS ARRAY_PARTITION variable=v


	roop0:
	for(i=0;i<MMAX;i++) {
        for(j=0;j<NMAX;j++) {
            a[i][j].real = a_real[i*NMAX+j];
            a[i][j].imag = a_imag[i*NMAX+j]; } }
	startt[0] = (char)a[i-1][j-1].imag;

    
    //double d_tmp;
    //double v_tmp;
    //double g_real, g_imag;
    //double h_real, h_imag;
    //double tmp_real, tmp_imag;
    
    //_FPU_SETCW(wctype);
    
    //n = amatrix_info.mic_num;
   // n = localize_music_info.tf_channel_num;
   //n=16;
    
    /* 単位行列の作成 */
    roop1:
    for(j=0; j<n; j++){
        for(i=0; i<n; i++){
            v[j][i].real = 0.0f;
            v[j][i].imag = 0.0f;
        }
        v[j][j].real = 1.0f;
    }

    roop2:
    for(ip=0; ip<n; ip++){
        //b[ip] = d[ip] = a(ip, ip).real();
        //z[ip] = 0;
        b[ip] = a[ip][ip].real;
        d[ip] = a[ip][ip].real;
        z[ip] = 0;
    }

    roop3:
    for(i=0; i<50; i++){
        sm = 0.0f;
        roopA:
        for(ip=0; ip<n-1; ip++){
//#pragma HLS UNROLL
            for(iq=ip+1; iq<n; iq++){
#pragma HLS PIPELINE it=1
                sm += sqrtf(a[ip][iq].real * a[ip][iq].real + a[ip][iq].imag * a[ip][iq].imag);
            }
        }
        // デバック
        //LOGI("i = %d : %g\n", i, sm);
        
        //if (sm == 0.0) 
        //if (sm <= __FLT_EPSILON__) 
        if(sm <= 1.0e-7){
            //LOGI("sm == 0.0\n");
            return;
        }
        
        //if(i<4){ thresh = 0.2f * sm / (float)(n * n); }
        //else{    thresh = 0.0f; }
        thresh = ((i>=4) ? 0.0f : (0.2f * sm / (float)(n * n)));
        //LOGI("thresh = %g\n", thresh);
        
        roopB:
        for(ip=0; ip<n-1; ip++){
            for(iq=ip+1; iq<n; iq++){
                float abs_a;
#pragma HLS PIPELINE it=1
                abs_a = sqrtf(a[ip][iq].real * a[ip][iq].real + a[ip][iq].imag * a[ip][iq].imag);
                g = 100.0f * abs_a;
                if ((i > 4) && (g <= 5.0e-9))
                {
                    a[ip][iq].real = 0.0f;
                    a[ip][iq].imag = 0.0f;
                    
                }
                else if (abs_a > thresh){
                    tmp = d[iq] - d[ip]; // hの実部(tmp)しか使っていない

                    if(g <= 1.0e-8){ // この時点では h は実部しか持たない
                        t.real =  a[ip][iq].real / tmp; // hの実部(tmp)しか使っていない
                        t.imag = -a[ip][iq].imag / tmp; // hの実部(tmp)しか使っていない

                    }
                    else{
                        norm_a  = a[ip][iq].real * a[ip][iq].real;
                        norm_a += a[ip][iq].imag * a[ip][iq].imag;
                        t_tmp  = -tmp + sqrtf(tmp * tmp + 4.0f * norm_a); // h.real=>(d[iq] - d[ip])=>tmp
                        t.real =  t_tmp * a[ip][iq].real / (2.0f * norm_a);
                        t.imag = -t_tmp * a[ip][iq].imag / (2.0f * norm_a);
                        
                    }
                    
                    norm_t = t.real * t.real + t.imag * t.imag;
                    c = 1.0 / sqrtf(1.0 + norm_t);
                    s.real = t.real * c;
                    s.imag = t.imag * c;
                    tau.real = s.real / (float)(1.0f + c);
                    tau.imag = s.imag / (float)(1.0f + c);
                    tmp  = t.real * a[ip][iq].real; // hの実部(tmp)しか使っていない
                    tmp -= t.imag * a[ip][iq].imag; // hの実部(tmp)しか使っていない
                    z[ip] -= tmp; // hの実部(tmp)しか使っていない
                    z[iq] += tmp; // hの実部(tmp)しか使っていない
                    d[ip] -= tmp; // hの実部(tmp)しか使っていない
                    d[iq] += tmp; // hの実部(tmp)しか使っていない
                    a[ip][iq].real = 0.0f;
                    a[ip][iq].imag = 0.0f;
                    
                    //デバック
                    //LOGI("h = %g , %g\n", h.real, h.imag);
                    //LOGI("%g , %g\n", s.real, s.imag);
                    //LOGI("%g , %g\n", tau.real, tau.imag);
                    
                    roopa:
                    for(j=0; j<=ip-1; j++){
//#pragma HLS UNROLL
                        ////rotate(a, j, ip, j, iq, c, s, tau);
                        //std::complex<T> g = a(j, ip);
                        //std::complex<T> h = a(j, iq);
                        g_real = a[j][ip].real;
                        g_imag = a[j][ip].imag;
                        h_real = a[j][iq].real;
                        h_imag = a[j][iq].imag;
                        //a(j, ip) = g - s * (h + g * conj(tau));
                        tmp_real  = h_real;
                        tmp_real += g_real * tau.real;
                        tmp_real += g_imag * tau.imag;
                        tmp_imag  = h_imag;
                        tmp_imag -= g_real * tau.imag;
                        tmp_imag += g_imag * tau.real;
                        tmp1[j]  = g_real;
                        tmp1[j] -= s.real * tmp_real;
                        tmp1[j] += s.imag * tmp_imag;
                        a[j][ip].real = tmp1[j];
                        tmp2[j]  = g_imag;
                        tmp2[j] -= s.real * tmp_imag;
                        tmp2[j] -= s.imag * tmp_real;
                        a[j][ip].imag = tmp2[j];
                        
                        //a(j, iq) = h + conj(s) * (g - h * tau);
                        tmp_real  = g_real;
                        tmp_real -= h_real * tau.real;
                        tmp_real += h_imag * tau.imag;
                        tmp_imag  = g_imag;
                        tmp_imag -= h_real * tau.imag;
                        tmp_imag -= h_imag * tau.real;
                        tmp1[j]  = h_real;
                        tmp1[j] += s.real * tmp_real;
                        tmp1[j] += s.imag * tmp_imag;
                        a[j][iq].real = tmp1[j];
                        tmp2[j]  = h_imag;
                        tmp2[j] += s.real * tmp_imag;
                        tmp2[j] -= s.imag * tmp_real;
                        a[j][iq].imag = tmp2[j];
                    }
                    roopb:
                    for(j=ip+1; j<=iq-1; j++){
#pragma HLS LOOP_TRIPCOUNT min=14 max=14
#pragma HLS UNROLL
                        ////rotate(a, ip, j, j, iq, c, s, tau);
                        //std::complex<T> g = a(ip, j);
                        //std::complex<T> h = a(j, iq);
                        g_real = a[ip][j].real;
                        g_imag = a[ip][j].imag;
                        h_real = a[j][iq].real;
                        h_imag = a[j][iq].imag;
                        //a(ip, j) = g - conj(s) * (conj(h) + g * tau);
                        tmp_real  =  h_real;
                        tmp_real +=  g_real * tau.real;
                        tmp_real -=  g_imag * tau.imag;
                        tmp_imag  = -h_imag;
                        tmp_imag +=  g_real * tau.imag;
                        tmp_imag +=  g_imag * tau.real;
                        tmp1[j]  = g_real;
                        tmp1[j] -= s.real * tmp_real;
                        tmp1[j] -= s.imag * tmp_imag;
                        a[ip][j].real = tmp1[j];
                        tmp2[j]  = g_imag;
                        tmp2[j] -= s.real * tmp_imag;
                        tmp2[j] += s.imag * tmp_real;
                        a[ip][j].imag = tmp2[j];
                        //a(j, iq) = h + conj(s) * (conj(g) - h * tau);
                        tmp_real  =  g_real;
                        tmp_real -=  h_real * tau.real;
                        tmp_real +=  h_imag * tau.imag;
                        tmp_imag  = -g_imag;
                        tmp_imag -= h_real * tau.imag;
                        tmp_imag -= h_imag * tau.real;
                        tmp1[j]  = h_real;
                        tmp1[j] += s.real * tmp_real;
                        tmp1[j] += s.imag * tmp_imag;
                        a[j][iq].real = tmp1[j];
                        tmp2[j]  = h_imag;
                        tmp2[j] += s.real * tmp_imag;
                        tmp2[j] -= s.imag * tmp_real;
                        a[j][iq].imag = tmp2[j];
                    }
                    roopc:
                    for(j=iq+1; j<n; j++){
#pragma HLS LOOP_TRIPCOUNT min=14 max=14
//#pragma HLS UNROLL
                        ////rotate(a, ip, j, iq, j, c, s, tau);
                        //std::complex<T> g = a(ip, j);
                        //std::complex<T> h = a(iq, j);
                        g_real = a[ip][j].real;
                        g_imag = a[ip][j].imag;
                        h_real = a[iq][j].real;
                        h_imag = a[iq][j].imag;
                        //a(ip, j) = g - conj(s) * (h + g * tau);
                        tmp_real  = h_real;
                        tmp_real += g_real * tau.real;
                        tmp_real -= g_imag * tau.imag;
                        tmp_imag  = h_imag;
                        tmp_imag += g_real * tau.imag;
                        tmp_imag += g_imag * tau.real;
                        tmp1[j]  = g_real;
                        tmp1[j] -= s.real * tmp_real;
                        tmp1[j] -= s.imag * tmp_imag;
                        a[ip][j].real = tmp1[j];
                        tmp2[j]  = g_imag;
                        tmp2[j] -= s.real * tmp_imag;
                        tmp2[j] += s.imag * tmp_real;
                        a[ip][j].imag = tmp2[j];
                        //a(iq, j) = h + s * (g - h * conj(tau));
                        tmp_real  = g_real;
                        tmp_real -= h_real * tau.real;
                        tmp_real -= h_imag * tau.imag;
                        tmp_imag  = g_imag;
                        tmp_imag += h_real * tau.imag;
                        tmp_imag -= h_imag * tau.real;
                        tmp1[j]  = h_real;
                        tmp1[j] += s.real * tmp_real;
                        tmp1[j] -= s.imag * tmp_imag;
                        a[iq][j].real = tmp1[j];
                        tmp2[j]  = h_imag;
                        tmp2[j] += s.real * tmp_imag;
                        tmp2[j] += s.imag * tmp_real;
                        a[iq][j].imag = tmp2[j];
                    }
                    roopd:
                    for(j=0; j<n; j++){
                        g_real = v[j][ip].real;
                        g_imag = v[j][ip].imag;
                        h_real = v[j][iq].real;
                        h_imag = v[j][iq].imag;
                        tmp_real  = h_real;
                        tmp_real += g_real * tau.real;
                        tmp_real += g_imag * tau.imag;
                        tmp_imag  = h_imag;
                        tmp_imag -= g_real * tau.imag;
                        tmp_imag += g_imag * tau.real;
                        tmp1[j]  = g_real;
                        tmp1[j] -= s.real * tmp_real;
                        tmp1[j] += s.imag * tmp_imag;
                        v[j][ip].real = tmp1[j];
                        tmp2[j]  = g_imag;
                        tmp2[j] -= s.real * tmp_imag;
                        tmp2[j] -= s.imag * tmp_real;
                        v[j][ip].imag = tmp2[j];
                        tmp_real  = g_real;
                        tmp_real -= h_real * tau.real;
                        tmp_real += h_imag * tau.imag;
                        tmp_imag  = g_imag;
                        tmp_imag -= h_real * tau.imag;
                        tmp_imag -= h_imag * tau.real;
                        tmp1[j]  = h_real;
                        tmp1[j] += s.real * tmp_real;
                        tmp1[j] += s.imag * tmp_imag;
                        v[j][iq].real = tmp1[j];
                        tmp2[j]  = h_imag;
                        tmp2[j] += s.real * tmp_imag;
                        tmp2[j] -= s.imag * tmp_real;
                        v[j][iq].imag = tmp2[j];
                    }
                }
            }
        }
        roopC:
        for(ip=0; ip<n; ip++){
            //LOGI("z[%d] = %g\n", ip, z[ip]);
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    
    // sorting eigenvalues in ascending order
    roop4:
    for(i=0; i<n-1; i++){
        float min = d[i];
        int imin = i;
        for(j=i+1; j<n; j++){
            if(min > d[j]){
                min = d[j];
                imin = j;
            }
        }
        if(imin != i){
            //std::swap(d[i], d[imin]);
            d_tmp   = d[i];
            d[i]    = d[imin];
            d[imin] = d_tmp;
            //ublas::column(v, i).swap(ublas::column(v, imin));
            for(j=0; j<n; j++){
                v_tmp           = v[j][i].real;
                v[j][i].real    = v[j][imin].real;
                v[j][imin].real = t_tmp;
                
                v_tmp           = v[j][i].imag;
                v[j][i].imag    = v[j][imin].imag;
                v[j][imin].imag = t_tmp;
            }
        }
    }
	stopt[0] = t_tmp;

	roop5:
	for(i=0;i<MMAX;i++)
        for(j=0;j<NMAX;j++) {
            v_real[i*NMAX+j]=v[i][j].real;
            v_imag[i*NMAX+j]=v[i][j].imag; }

}
