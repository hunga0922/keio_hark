
///////////////////////////////////////////////////////////////
//  Include Files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MMAX 16
#define NMAX 16

typedef struct
{   
    //// cancel //// int device_type;
    int  length;                     //
    int  advance;                    //  _ADVANCE_VAL_
    int  channel_count;              //  _DEV_CH_VAL_
    int  logical_channel_count;      //  _CH_VAL_
    int  channel_offset;             //  _CH_OFFSET_VAL_     for Easy implementation of ChannelSelector like function.
    int  sampling_rate;              //
    int  nbits;                      //  _BIT_VAL_           Number of sample bits for audio device.
} Param_AudioStreamFromMic_t;
typedef struct _COMPLEX_FLOAT_
{
    float real;
    float imag;
} COMPLEX_FLOAT;

float ABS(float a)
{
    return (a < 0.0f ? -a : a > 0.0f ? a : 0.0f);
}


void Svd(float a_real[MMAX*NMAX], float a_imag[MMAX*NMAX],
	float s[NMAX], float u_real[MMAX*NMAX], float u_imag[MMAX*NMAX],
	float v_real[MMAX*NMAX], float v_imag[MMAX*NMAX])
{
	COMPLEX_FLOAT a[MMAX][NMAX];
	COMPLEX_FLOAT u[MMAX][NMAX];
	COMPLEX_FLOAT v[MMAX][NMAX];
    int m;
    int n;
    int p;  // transformation Conj(Tran(u)) is applied to the p vectors given in columns n, n+1, ..., n+p-1 of matrix a
    int nu; // considered number of cols of u
    int nv; // considered number of cols of v
    
    float b[16];
    float c[16];
    float t[16];
    
    COMPLEX_FLOAT q;
    COMPLEX_FLOAT r;
    int i, j, k, k1, L, L1, nM1, np;
    
    float cs, eps, eta, f, g, h;
    float sn;
    float tol, w, x, y, z;
    float tmp, tmp2, tmp_real, tmp_imag;
//  Input
	for(i=0;i<MMAX;i++) {
		for(j=0;j<NMAX;j++) {
			a[i][j].real = a_real[i*NMAX+j];
			a[i][j].imag = a_imag[i*NMAX+j]; } }
//    m = localize_music_info.tf_channel_num;
//   n = localize_music_info.tf_channel_num;
    m = 16;
    n = 16;
    p = 0;
    nu = m;
    nv = n;
    
    eta = 1.5E-07F;
    tol = 1.5E-31F;
    np  = n + p;
    nM1 = n - 1;
    L   = 0;
    // HOUSEHOLDER REDUCTION
    c[0] = 0.0F;
    k = 0;
    
    while(1){
        k1 = k + 1;
        // ELIMINATION OF a(i,k), i = k, ..., m-1
        z = 0.0F;
        
        for(i = k; i < m; i++){
    	    //z += norm(a(i,k));
    	    z += a[i][k].real * a[i][k].real;
    	    z += a[i][k].imag * a[i][k].imag;
    	}
    	
        b[k] = 0.0F;
        if(z > tol){
            //z = (float)sqrt(z);
            z = sqrtf(z);
            b[k] = z;
            
            //w = abs(a(k,k));
            //w = sqrt(a[k][k].real * a[k][k].real + a[k][k].imag * a[k][k].imag);
            w = sqrtf(a[k][k].real * a[k][k].real + a[k][k].imag * a[k][k].imag);
            
            if(w != 0.0F){
                //q = a(k,k) / w;
                //q.real = a[k][k].real / w;
                //q.imag = a[k][k].imag / w;
                tmp = 1.0f / w;
                q.real = a[k][k].real * tmp;
                q.imag = a[k][k].imag * tmp;
            }
            else{
                q.real = 1.0f;
                q.imag = 0.0f;
            }
            //a(k,k) = q * (z + w);
            //a[k][k].real = q.real * (z + w);
            //a[k][k].imag = q.imag * (z + w);
            tmp = (z + w);
            a[k][k].real = q.real * tmp;
            a[k][k].imag = q.imag * tmp;
            tmp = 1.0f / (tmp * z);
            
            if(k != np - 1){
                for(j=k1; j<np; j++){
                    //q = std::complex<T>(0.0F, 0.0F);
                    q.real = 0.0f;
                    q.imag = 0.0f;
                    
                    for (i=k; i<m; i++){
                        //q += conj(a(i,k)) * a(i,j);
                        q.real += a[i][k].real * a[i][j].real;
                        q.real += a[i][k].imag * a[i][j].imag;
                        q.imag += a[i][k].real * a[i][j].imag;
                        q.imag -= a[i][k].imag * a[i][j].real;
                    }
                    
                    //q /= z * (z + w);
                    ////q.real = q.real / (z * (z + w));
                    ////q.imag = q.imag / (z * (z + w));
                    q.real *= tmp;
                    q.imag *= tmp;
                    for(i=k; i<m; i++){
                        //a(i,j) -= q * a(i,k);
                        tmp_real  = q.real * a[i][k].real;
                        tmp_real -= q.imag * a[i][k].imag;
                        tmp_imag  = q.real * a[i][k].imag;
                        tmp_imag += q.imag * a[i][k].real;
                        a[i][j].real -= tmp_real;
                        a[i][j].imag -= tmp_imag;
                    }
                }
            }
            // PHASE TRANSFORMATION
            //q = -conj(a(k,k)) / abs(a(k,k));
            //w = sqrt(a[k][k].real * a[k][k].real + a[k][k].imag * a[k][k].imag);
            w = sqrtf(a[k][k].real * a[k][k].real + a[k][k].imag * a[k][k].imag);
            //q.real = -a[k][k].real / w;
            //q.imag =  a[k][k].imag / w;
            tmp = 1.0f / w;
            q.real = -a[k][k].real * tmp;
            q.imag =  a[k][k].imag * tmp;
            
            for(j=k1; j<np; j++){
                //a(k,j) *= q;
                tmp_real  = a[k][j].real * q.real;
                tmp_real -= a[k][j].imag * q.imag;
                tmp_imag  = a[k][j].real * q.imag;
                tmp_imag += a[k][j].imag * q.real;
                a[k][j].real = tmp_real;
                a[k][j].imag = tmp_imag;
            }
        }
        // ELIMINATION OF a(k][j], j = k+2, ..., n-1
        if (k == nM1) break;
        z = 0.0F;
        for(j=k1; j<n; j++){
            //z += norm(a(k,j));
            z += a[k][j].real * a[k][j].real;
            z += a[k][j].imag * a[k][j].imag;
        }
        
        c[k1] = 0.0F;
        if(z > tol){
            //z = (float)sqrt(z);
            z = sqrtf(z);
            c[k1] = z;
            
            //w = abs(a(k,k1));
            //w = sqrt(a[k][k1].real * a[k][k1].real + a[k][k1].imag * a[k][k1].imag);
            w = sqrtf(a[k][k1].real * a[k][k1].real + a[k][k1].imag * a[k][k1].imag);
            
            if(w != 0.0F){
                //q = a(k,k1) / w;
                //q.real = a[k][k1].real / w;
                //q.imag = a[k][k1].imag / w;
                tmp = 1.0f / w;
                q.real = a[k][k1].real * tmp;
                q.imag = a[k][k1].imag * tmp;
            }
            else{
                q.real = 1.0f;
                q.imag = 0.0f;
            }
            //a(k,k1) = q * (z + w);
            tmp = (z + w);
            ////a[k][k1].real = q.real * (z + w);
            ////a[k][k1].imag = q.imag * (z + w);
            a[k][k1].real = q.real * tmp;
            a[k][k1].imag = q.imag * tmp;
            tmp = 1.0f / (tmp * z);
            
            for(i=k1; i<m; i++){
                //q = std::complex<T>(0.0F, 0.0F);
                q.real = 0.0f;
                q.imag = 0.0f;
                
                for(j=k1; j<n; j++){
	                //q += conj(a(k,j)) * a(i,j);
	                q.real += a[k][j].real * a[i][j].real;
	                q.real += a[k][j].imag * a[i][j].imag;
	                q.imag += a[k][j].real * a[i][j].imag;
	                q.imag -= a[k][j].imag * a[i][j].real;
                }
                
                //q /= z * (z + w);
                ////q.real = q.real / (z * (z + w));
                ////q.imag = q.imag / (z * (z + w));
                q.real *= tmp;
                q.imag *= tmp;
                for(j=k1; j<n; j++){
	                //a(i,j) -= q * a(k,j);
	                tmp_real  = q.real * a[k][j].real;
	                tmp_real -= q.imag * a[k][j].imag;
	                tmp_imag  = q.real * a[k][j].imag;
	                tmp_imag += q.imag * a[k][j].real;
	                a[i][j].real -= tmp_real;
	                a[i][j].imag -= tmp_imag;
	            }
            }
            // PHASE TRANSFORMATION
            //q = -conj(a(k,k1)) / abs(a(k,k1));
            //w = sqrt(a[k][k1].real * a[k][k1].real + a[k][k1].imag * a[k][k1].imag);
            w = sqrtf(a[k][k1].real * a[k][k1].real + a[k][k1].imag * a[k][k1].imag);
            q.real = -a[k][k1].real / w;
            q.imag =  a[k][k1].imag / w;
            for(i=k1; i<m; i++){
                //a(i,k1) *= q;
                tmp_real  = a[i][k1].real * q.real;
                tmp_real -= a[i][k1].imag * q.imag;
                tmp_imag  = a[i][k1].real * q.imag;
                tmp_imag += a[i][k1].imag * q.real;
                a[i][k1].real = tmp_real;
                a[i][k1].imag = tmp_imag;
            }
        }
        k = k1;
    } /* End of while(1) */
    
    // TOLERANCE FOR NEGLIGIBLE ELEMENTS
    eps = 0.0F;
    for(k=0; k<n; k++){
        s[k] = b[k];
        t[k] = c[k];
        if(s[k] + t[k] > eps){
            eps = s[k] + t[k];
        }
    }
    eps *= eta;
    
    // INITIALIZATION OF u AND v
    if(nu > 0){
        for(i=0; i<m; i++){
            for(j=0; j<nu; j++){
                //u(i,j) = std::complex<T>(0.0F, 0.0F);
                u[i][j].real = 0.0f;
                u[i][j].imag = 0.0f;
            }
            //u(j,j) = std::complex<T>(1.0F, 0.0F);
            //u[j][j].real = 1.0f; // skip because already set! u[j][j].imag = 0.0f;
            u[i][i].real = 1.0f; // swap loop i,j;
        }
    }
    
    if(nv > 0){
        for(i=0; i<n; i++){
            for(j=0; j<nv; j++){
                //v(i,j) = std::complex<T>(0.0F, 0.0F);
                v[i][j].real = 0.0f;
                v[i][j].imag = 0.0f;
            }
            //v(j,j) = std::complex<T>(1.0F, 0.0F);
            //v[j][j].real = 1.0f; // skip because already set! v[j][j].imag = 0.0f;
            v[i][i].real = 1.0f; // swap loop i,j;
        }
    }
    //QR DIAGONALIZATION
    for(k=nM1; k>=0; k--){
        //デバック
        //LOGI("k = %d\n", k);

        //TEST FOR SPLIT
        while(1){
            for(L=k; L>=0; L--){
                //デバック
	            //LOGI("%d : %g, %g, eps = %g\n", L, t[L], s[L-1], eps);
                if (ABS(t[L]) <= eps) goto Test;
                if(L > 0){
                    if (ABS(s[L - 1]) <= eps) break;
                }
            }

            //CANCELLATION OF E(L)
            cs = 0.0F;
            sn = 1.0F;
            L1 = L - 1;
            for(i=L; i<=k; i++){
                f = sn * t[i];
                t[i] *= cs;
                
                if (ABS(f) <= eps) goto Test;
                
                h = s[i];
                //w = (float)sqrt(f * f + h * h);
                w = sqrtf(f * f + h * h);
                s[i] = w;
                //cs = h / w;
                //sn = -f / w;
                tmp = 1.0f / w;
                cs = h * tmp;
                sn = -f * tmp;
                
                if(nu > 0){
                    for(j=0; j<n; j++){
                        //x = real(u(j,L1));
                        //y = real(u(j,i));
                        x = u[j][L1].real;
                        y = u[j][i].real;
                        
                        //u(j,L1) = std::complex<T>(x * cs + y * sn, 0.0F);
                        tmp  = x * cs;
                        tmp += y * sn;
                        u[j][L1].real = tmp;
                        u[j][L1].imag = 0.0f;
                        
                        //u(j,i) = std::complex<T>(y * cs - x * sn, 0.0F);
                        tmp  = y * cs;
                        tmp -= x * sn;
                        u[j][i].real = tmp;
                        u[j][i].imag = 0.0f;
                    }
	            }
	            
                if(np == n) continue;
                
                for(j=n; j<np; j++){
	                //q = a(L1,j);
	                q = a[L1][j];
	                //r = a(i,j);
	                r = a[i][j];
	                
	                //a(L1,j) = q * cs + r * sn;
	                tmp  = q.real * cs;
	                tmp += r.real * sn;
	                a[L1][j].real = tmp;
	                tmp  = q.imag * cs;
	                tmp += r.imag * sn;
	                a[L1][j].imag = tmp;
	                
	                //a(i,j) = r * cs - q * sn;
	                tmp  = r.real * cs;
	                tmp -= q.real * sn;
	                a[i][j].real = tmp;
	                tmp  = r.imag * cs;
	                tmp -= q.imag * sn;
	                a[i][j].imag = tmp;
	            }
            }

            //TEST FOR CONVERGENCE
Test:	    w = s[k];
            if (L == k) break;

            //ORIGIN SHIFT
            x = s[L];
            y = s[k - 1];
            g = t[k - 1];
            h = t[k];
            f = ((y - w) * (y + w) + (g - h) * (g + h)) / (2.0F * h * y);
            //g = (float)sqrt(f * f + 1.0F);
            g = sqrtf(f * f + 1.0F);
            if(f < 0.0F){ g = -g; }
            f = ((x - w) * (x + w) + (y / (f + g) - h) * h) / x;
  
            //QR STEP
            cs = 1.0F;
            sn = 1.0F;
            L1 = L + 1;
            for(i=L1; i<=k; i++){
                g = t[i];
                y = s[i];
                h = sn * g;
                g = cs * g;
                //w = (float)sqrt(h * h + f * f);
                w = sqrtf(h * h + f * f);
                t[i - 1] = w;
                //cs = f / w;
                //sn = h / w;
                tmp = 1.0f / w;
                cs = f * tmp;
                sn = h * tmp;
                f  = x * cs;
                f += g * sn;
                g  = g * cs;
                g -= x * sn;
                h  = y * sn;
                y  = y * cs;
                if(nv > 0){
	                for(j=0; j<n; j++){
                        //x = real(v(j,i - 1));
                        //w = real(v(j,i));
                        x = v[j][i-1].real;
                        w = v[j][i].real;
                        
                        //v(j,i - 1) = std::complex<T>(x * cs + w * sn, 0.0F);
                        tmp  = x * cs;
                        tmp += w * sn;
                        v[j][i-1].real = tmp;
                        v[j][i-1].imag = 0.0f;
                        //v(j,i) = std::complex<T>(w * cs - x * sn, 0.0F);
                        tmp  = w * cs;
                        tmp -= x * sn;
                        v[j][i].real = tmp;
                        v[j][i].imag = 0.0f;
                    }
                }
                //w = (float)sqrt(h * h + f * f);
                w = sqrtf(h * h + f * f);
                s[i - 1] = w;
                //cs = f / w;
                //sn = h / w;
                tmp = 1.0f / w;
                cs = f * tmp;
                sn = h * tmp;
                f = cs * g + sn * y;
                x = cs * y - sn * g;
                if(nu > 0){
                    for(j=0; j<n; j++){
                        //y = real(u(j,i - 1));
                        //w = real(u(j,i));
                        y = u[j][i-1].real;
                        w = u[j][i].real;
                        
                        //u(j,i - 1) = std::complex<T>(y * cs + w * sn, 0.0F);
                        tmp  = y * cs;
                        tmp += w * sn;
                        u[j][i-1].real = tmp;
                        u[j][i-1].imag = 0.0f;
                        //u(j,i) = std::complex<T>(w * cs - y * sn, 0.0F);
                        tmp  = w * cs;
                        tmp -= y * sn;
                        u[j][i].real = tmp;
                        u[j][i].imag = 0.0f;
                    }
                }
	            
	            if (n == np) continue;
	            
	            for(j=n; j<np; j++){
                    //q = a(i - 1,j);
                    //r = a(i,j);
                    q = a[i-1][j];
                    r = a[i][j];
                    
                    //a(i - 1,j) = q * cs + r * sn;
                    tmp  = q.real * cs;
                    tmp += r.real * sn;
                    a[i-1][j].real = tmp;
                    tmp  = q.imag * cs;
                    tmp += r.imag * sn;
                    a[i-1][j].imag = tmp;
                    //a(i,j) = r * cs - q * sn;
                    tmp  = r.real * cs;
                    tmp -= q.real * sn;
                    a[i][j].real = tmp;
                    tmp  = r.imag * cs;
                    tmp -= q.imag * sn;
                    a[i][j].imag = tmp;
                }
            }
	        t[L] = 0.0F;
	        t[k] = f;
	        s[k] = x;
	    } /* End of while(1) */
	    
        //CONVERGENCE
        if (w >= 0.0F) continue;
        
        s[k] = -w;
        
        if (nv == 0) continue;
        
        for(j=0; j<n; j++){
	        //v(j,k) = -v(j,k);
	        v[j][k].real = -v[j][k].real;
	        v[j][k].imag = -v[j][k].imag;
	    }
    }
    
    //SORT SINGULAR VALUES
    for(k=0; k<n; k++){	/* sort descending */
        g = -1.0F;
        j = k;
        for(i=k; i<n; i++){	/* sort descending */
            if (s[i] <= g) continue;
            
            g = s[i];
            j = i;
        }
        
        if (j == k) continue;
        
        s[j] = s[k];
        s[k] = g;
        if(nv > 0){
            for(i=0; i<n; i++){
                //q = v(i,j);
                //v(i,j) = v(i,k);
                //v(i,k) = q;
                q = v[i][j];
                v[i][j] = v[i][k];
                v[i][k] = q;
            }
        }
        
        if(nu > 0){
            for(i=0; i<n; i++){
                //q = u(i,j);
                //u(i,j) = u(i,k);
                //u(i,k) = q;
                q = u[i][j];
                u[i][j] = u[i][k];
                u[i][k] = q;
            }
        }
        
        if (n == np) continue;
        
        for(i=n; i<np; i++){
            //q = a(j,i);
            //a(j,i) = a(k,i);
            //a(k,i) = q;
            q = a[j][i];
            a[j][i] = a[k][i];
            a[k][i] = q;
        }
    }

    //BACK TRANSFORMATION
    if(nu > 0){
        for(k=nM1; k>=0; k--){
            if (b[k] == 0.0F) continue;
            
            //w = sqrt(a[k][k].real * a[k][k].real + a[k][k].imag * a[k][k].imag);
            w = sqrtf(a[k][k].real * a[k][k].real + a[k][k].imag * a[k][k].imag);
            //q = -a(k,k) / abs(a(k,k));
            //q.real = -a[k][k].real / w;
            //q.imag = -a[k][k].imag / w;
            tmp2 = 1.0f / (w * b[k]);
            tmp = 1.0f / w;
            q.real = -a[k][k].real * tmp;
            q.imag = -a[k][k].imag * tmp;
            
            for(j=0; j<nu; j++){
                //u(k,j) *= q;
                tmp_real  = u[k][j].real * q.real;
                tmp_real -= u[k][j].imag * q.imag;
                tmp_imag  = u[k][j].real * q.imag;
                tmp_imag += u[k][j].imag * q.real;
                u[k][j].real = tmp_real;
                u[k][j].imag = tmp_imag;
            }
            
            for(j=0; j<nu; j++){
                //q = std::complex<T>(0.0F, 0.0F);
                q.real = 0.0f;
                q.imag = 0.0f;
                for(i=k; i<m; i++){
                    //q += conj(a(i,k)) * u(i,j);
                    q.real += a[i][k].real * u[i][j].real;
                    q.real += a[i][k].imag * u[i][j].imag;
                    q.imag += a[i][k].real * u[i][j].imag;
                    q.imag -= a[i][k].imag * u[i][j].real;
                }
                
                //w = sqrt(a[k][k].real * a[k][k].real + a[k][k].imag * a[k][k].imag);
                //w = sqrtf(a[k][k].real * a[k][k].real + a[k][k].imag * a[k][k].imag); // このループ内でaは更新しないので不要
                //q /= abs(a(k,k)) * b[k];
                //q.real = q.real / (w * b[k]);
                //q.imag = q.imag / (w * b[k]);
                q.real = q.real * tmp2;
                q.imag = q.imag * tmp2;
                
                for(i=k; i<m; i++){
                    //u(i,j) -= q * a(i,k);
                    tmp  = q.real * a[i][k].real;
                    tmp -= q.imag * a[i][k].imag;
                    u[i][j].real -= tmp;
                    tmp  = q.real * a[i][k].imag;
                    tmp += q.imag * a[i][k].real;
                    u[i][j].imag -= tmp;
                }
            }
        }
    }
    
    if(nv > 0 && n > 1){
        for(k=n-2; k>=0; k--){
            k1 = k + 1;
            if (c[k1] == 0.0F) continue;
            
            //w = sqrt(a[k][k1].real * a[k][k1].real + a[k][k1].imag * a[k][k1].imag);
            w = sqrtf(a[k][k1].real * a[k][k1].real + a[k][k1].imag * a[k][k1].imag);
            //q = -conj(a(k,k1)) / abs(a(k,k1));
            //q.real = -a[k][k1].real / w;
            //q.imag =  a[k][k1].imag / w;
            tmp2 = 1.0f / (w * c[k1]);
            tmp = 1.0f / w;
            q.real = -a[k][k1].real * tmp;
            q.imag =  a[k][k1].imag * tmp;
            
            for(j=0; j<nv; j++){
                //v(k1,j) *= q;
                tmp_real  = v[k1][j].real * q.real;
                tmp_real -= v[k1][j].imag * q.imag;
                tmp_imag  = v[k1][j].real * q.imag;
                tmp_imag += v[k1][j].imag * q.real;
                v[k1][j].real = tmp_real;
                v[k1][j].imag = tmp_imag;
            }
            
            for(j=0; j<nv; j++){
                //q = std::complex<T>(0.0F, 0.0F);
                q.real = 0.0f;
                q.imag = 0.0f;
                
                for(i=k1; i<n; i++){
                    //q += a(k,i) * v(i,j);
                    q.real += a[k][i].real * v[i][j].real;
                    q.real -= a[k][i].imag * v[i][j].imag;
                    q.imag += a[k][i].real * v[i][j].imag;
                    q.imag += a[k][i].imag * v[i][j].real;
                }
                
                //w = sqrt(a[k][k1].real * a[k][k1].real + a[k][k1].imag * a[k][k1].imag);
                //w = sqrtf(a[k][k1].real * a[k][k1].real + a[k][k1].imag * a[k][k1].imag); // このループ内でaは更新しないので不要
                //q = q / (abs(a(k,k1)) * c[k1]);
                //q.real = q.real / (w * c[k1]);
                //q.imag = q.imag / (w * c[k1]);
                q.real = q.real * tmp2;
                q.imag = q.imag * tmp2;
                
                for(i=k1; i<n; i++){
                    //v(i,j) -= q * conj(a(k,i));
                    tmp  =  q.real * a[k][i].real;
                    tmp +=  q.imag * a[k][i].imag;
                    v[i][j].real -= tmp;
                    tmp  = -q.real * a[k][i].imag;
                    tmp +=  q.imag * a[k][i].real;
                    v[i][j].imag -= tmp;
                }
            }
        }
    }
	for(i=0;i<MMAX;i++) 
		for(j=0;j<NMAX;j++) {
			u_real[i*NMAX+j]=u[i][j].real;
			u_imag[i*NMAX+j]=u[i][j].imag;
			v_real[i*NMAX+j]=v[i][j].real;
			v_imag[i*NMAX+j]=v[i][j].imag; }

}
