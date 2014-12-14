/* gauss-jordan elimination in finite field (Zp) 
   (semi-sparse version)
   by Martin Bujnak, mar2008
   last edit by Pavel Trutman, oct 2014
*/

#include "mex.h"
#include "math.h"
#include <string.h>

void gcd(int a, int b, int *g, int *c, int *d)
{
    int q;
    int u[3];
    int v[3];
    int t[3];
    
    u[0] = 1;
    u[1] = 0;
    u[2] = abs(a);

    v[0] = 0;
    v[1] = 1;
    v[2] = b;
    
    while (v[2]) {

       q = u[2]/v[2];
       
       t[0] = u[0] - v[0]*q;
       t[1] = u[1] - v[1]*q;
       t[2] = u[2] - v[2]*q;
       
       u[0] = v[0];
       u[1] = v[1];
       u[2] = v[2];
       
       v[0] = t[0];
       v[1] = t[1];
       v[2] = t[2];
    }

    if (a < 0) u[0] = -u[0];
    
    *c = u[0];
    *d = u[1];
    *g = u[2];
}

int invZp(int x, int prime)
{
    int g, c, d;
    gcd(x, prime, &g, &c, &d);
    return c;
}

int max(int a, int b) {
  if (a >= b) {
    return a;
  }
  else {
    return b;
  }
}

void mexFunction(int nlhs , mxArray *plhs[] , int nrhs , const mxArray *prhs[]){

    /*
       B = gjzp(A, prime);
    */
    
    /* in */
    int r = 0;      /* row */
    int c = 0;      /* col */
    int k;
    int l;
    int b;
    int dstofs;
    int srcofs;    
    int ofs = 0;
    int pofs = 0;
    int pivot_i;
    int endcol = 0;
    int *tmp;
    int *B;
    double *Bres;
    int *r_nonzero;
    int *crow;
    double *A = mxGetPr(prhs[0]);
    int rcnt = mxGetM(prhs[0]);
    int ccnt = mxGetN(prhs[0]);
    int prime  = (int)*mxGetPr(prhs[1]);
    
    /* out */
    plhs[0] = mxCreateDoubleMatrix(rcnt, ccnt , mxREAL);
    Bres = mxGetPr(plhs[0]);
    B = (int*)malloc(rcnt*ccnt*sizeof(int));

    /* copy transposed  */
    pofs = 0;
    for (k=0;k<ccnt;k++) {
        
        ofs = k;
        for (l=0;l<rcnt;l++) {
            
            *(B+ofs) = *(A+pofs);
            pofs++;
            ofs += ccnt;
        }
    }

    tmp = (int*)malloc(ccnt*sizeof(int));

    r_nonzero = (int*)malloc(2*rcnt*sizeof(int));
    crow = r_nonzero;

    /* find first and the last nonzero element for each row */
    for (k=0; k < rcnt; k++) {
	
        l = 0;
        ofs = k*ccnt;
        while (*(B+ofs) == 0) {
            ofs++;
            l++;
        };
        crow[0] = l;

        l = ccnt;
        ofs = (k+1)*ccnt - 1;
        while (l > crow[0] && *(B+ofs) == 0) {
            ofs--;
            l--;
        };
        crow[1] = l;

        crow += 2;
    }

    /* gj */
    ofs = 0;
    pofs = 0;
    while (r < rcnt && c < ccnt) {
        
        /* find pivot (in Zp == find first nonzero element) */
        int pivot = 0;
        int pivot_r = -1;
        int pofs = ofs;
        for (k = r; k < rcnt; k++) {
            
            if (*(B+pofs)) {

                pivot = *(B+pofs);
                pivot_r = k;
                break;
            }
            pofs += ccnt;
        }
        
        if (pivot == 0) {
            
            /* shift to next col */
            c++;
            ofs++;
            
        } else {
            
            /* exchange pivot and selected rows */
			endcol = max(*(r_nonzero+2*pivot_r+1), *(r_nonzero+2*r+1));
			k = (endcol - c)*sizeof(int);
            
            memcpy(tmp, B+ccnt*r+c, k);
            memcpy(B+ccnt*r+c, B+ccnt*pivot_r+c, k);
            memcpy(B+ccnt*pivot_r+c, tmp, k);

			endcol = *(r_nonzero+2*pivot_r+1);
           
            /* swap row limits */
            pofs = *(r_nonzero+2*pivot_r);
            *(r_nonzero+2*pivot_r) = *(r_nonzero+2*r);
            *(r_nonzero+2*r) = pofs;
            
            pofs = *(r_nonzero+2*pivot_r+1);
            *(r_nonzero+2*pivot_r+1) = *(r_nonzero+2*r+1);
            *(r_nonzero+2*r+1) = pofs;
           
            /* process rows */
            
            pivot_i = invZp(pivot, prime);
            
            /* divide pivot row */
            srcofs = ofs;
			*(B+srcofs) = 1;
			srcofs++;
            for (l = c+1; l < endcol; l++) {

                *(B+srcofs) = (*(B+srcofs)*pivot_i) % prime;
                srcofs++;
            }
            
            /* zero bottom */
            pofs = ofs + ccnt;
            for (k = r + 1; k < rcnt; k++) {
                
                if (*(B+pofs) != 0) {
                    
                    /* nonzero row */
                    b = *(B+pofs);

                    dstofs = pofs + 1;
                    srcofs = ofs + 1;

                    for (l = c + 1; l < endcol; l++) {
                        
                        int v = (*(B+dstofs) - *(B+srcofs) * b) % prime;
                        *(B+dstofs) = v;
                        dstofs++;
                        srcofs++;
                    }
                   
                    *(B+pofs) = 0;
                    if (*(r_nonzero+2*k+1) < endcol) *(r_nonzero+2*k+1) = endcol;
                }
                pofs += ccnt;
            }
            
            /* zero top */
            pofs = c;
            for (k = 0; k < r; k++) {
                
                if (*(B+pofs) != 0) {
                    
                    /* nonzero row */
                    b = *(B+pofs);

                    dstofs = pofs + 1;
                    srcofs = ofs + 1;

                    for (l = c + 1; l < endcol; l++) {
                        
                        int v = (*(B+dstofs) - *(B+srcofs) * b) % prime;
                        *(B+dstofs) = v;
                        dstofs++;
                        srcofs++;
                    }
                    
                    *(B+pofs) = 0;
                    if (*(r_nonzero+2*k+1) < endcol) *(r_nonzero+2*k+1) = endcol;
                }
                pofs += ccnt;
            }            

            r++;
            c++;
            ofs += ccnt + 1;
            
        }
    }
    
    /* transpose back to matlab format */
    pofs = 0;
    for (k=0;k<ccnt;k++) {
        
        ofs = k;
        for (l=0;l<rcnt;l++) {
            
            if (*(B+ofs) < 0) *(Bres+pofs) = *(B+ofs) + prime;
            else *(Bres+pofs) = *(B+ofs);
            pofs++;
            ofs += ccnt;
        }
    }
    
    free(tmp);
    free(r_nonzero);
    free(B);
}
