/*
 * =============================================================
 * nlemans_weight_sym.c
 *Input: a
 *X: .................... Input image
 *h: .................... Smooth Parameter
 *nwin: ................. half of patch size
 *bloc .................. half of block size
 *  patch size: [2*nwin+1, 2*nwin+1],
 * %  search window: [2*nbloc+1, 2*nbloc+1]
 *
 * This is a MEX-file for MATLAB.
 * =============================================================
 */

/* Revision: 1.0, change the patch and bloc size to half  */
/*Copyright(C) Xiaoqun zhang 16/04/2007*/

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>



#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
/**************************************************************************/
void quickSort(double *a,int left,int right)
{
    int i,j;
    double temp;
    
    i=left;
    j=right;
    temp=a[left];
    if(left>right)
        return;
    while(i!=j)
    {
        while(a[j]>=temp && j>i)
            j--;
        if(j>i)
            a[i++]=a[j];
        while(a[i]<=temp && j>i)
            i++;
        if(j>i)
            a[j--]=a[i];
        
    }
    a[i]=temp;
    quickSort(a,left,i-1);
    quickSort(a,i+1,right);
}
/**************************************************************************/
double  weightswidth(double *a, double v, int wd)
{
    int i;
    double ar, cs, cs2, mt;
    ar = 1;
    for(cs=0,cs2=v*(v),i=0;i<wd;i++)
    {
        cs+=a[i];
        cs2+=a[i]*(a[i]);
        if(cs>0.1)
        {
            mt=(double)(cs2)/cs;
            if(mt<a[i])
                break;
            else
                ar=mt;
        }
    }
    //ar=MAX(ar,0.1);
    return(ar);
}
/**************************************************************************/
void avsqrt_image(double *noisy, double *av,int nwp,int nx,int ny,int nw,int nwnw,int nxsy)
{
    int    x,y,k,adrp,i,xp,yp;
    double sum;
    for(k=0,x=nwp;x<nx+nwp;x++)
    {
        for(y=nwp;y<ny+nwp;y++,k++)
        {
            for(sum=0,i=0,xp=x-nw;xp<=x+nw;xp++)
            {
                for(yp=y-nw;yp<=y+nw;yp++,i++)
                {
                    adrp  = yp*nxsy+xp;
                    sum  += noisy[adrp];
                    //printf("sum=%f\n",sum);
                    //printf("noisy[adrp]=%f\n",noisy[adrp]);
                }
            }
            sum   = MIN( MAX( (sum/nwnw),0 ),255 );
            av[k] = sqrt(sum);
            //printf("av=%f\n",av[k]);
        }
    }
    
}
/**************************************************************************/
void owf_denoising(double *noisy,
        double *denoisy,
        double v,
        int    nw,
        int    np,
        int    nx,
        int    ny,
        double nlh
        )
{
    double *synoisy,*avimage;
    int    nwp= nw+np;
    int    nxsy, nysy, nxnysy;
    int    x,y,xp,yp,i,j,k;
    double *cw,*w;
    int    *dadr,*dd;
    int    nw2,np2,nwnw,npnp,nxny;
    int    adr,adrp;
    double sum,dist,dist2,ar,wrho;
    double *rho,*srho;
    /***************************************************************/
    nxny    = nx*ny;
    nxsy    = nx+2*nwp;
    nysy    = ny+2*nwp;
    nxnysy  = nysy*nxsy;
    synoisy = (double *)malloc(nxnysy*sizeof(double));
    avimage =  (double *)malloc(nxny*sizeof(double));
    printf("nlh=%f\n",nlh);
    
    
    for (k=0,i=0;i<nxsy;i++)
    {
        if (i<nwp)
            x=nwp-i;
        else if (i>nx+nwp-1)
            x=2*nx+nwp-i-2;
        else x=i-nwp;
        //printf("%d\n",x);
        for (j=0;j<nysy;j++,k++)
        {
            if (j<nwp )
                y=nwp-j;
            else if (j>nx+nwp-1)
                y=2*nx+nwp-j-2;
            else y=j-nwp;
            synoisy[k]=noisy[y*nx+x];
            //printf("   %f",synoisy[k]);
            //printf("k=%d\n",k);
        }
    }
    
    
    nw2   = 2*nw+1;
    nwnw  = nw2*nw2;
    //printf("nwnw=%d\n",nwnw);
    np2   = 2*np+1;
    npnp  = np2*np2;
    //printf("npnp=%d\n",npnp);
    w     = (double*)malloc(npnp*sizeof(double));
    dadr  = (int*)malloc(npnp*sizeof(int));
    cw    = (double*)malloc((np+1)*sizeof(double));
    //printf("npnp=%d\n",npnp);
    
    for(i=1;i<=np;i++)
    {
        for  (cw[i]=0,j=i;j<=np;j++)
        {
            cw[i]+=1./( np*(2*j+1)*(2*j+1) );
            
        }
    }
    cw[0]  = cw[1];
    
    for(i=0, x=-np; x<=np;x++)
        for(y=-np;y<=np;y++,i++)
        {
            dadr[i]=y*nxsy+x;
            j=MAX(abs(x),abs(y));
            w[i]=cw[j];
        }
    free(cw);
    
    //main loop
    rho    = (double*)malloc(nwnw*sizeof(double));
    srho   = (double*)malloc(nwnw*sizeof(double));
    avsqrt_image(synoisy,avimage,nwp, nx, ny, nw, nwnw, nxsy);
    for(k=0,x=nwp;x<nx+nwp;x++)
    {
        for(y=nwp;y<ny+nwp;y++,k++)
        {
            adr = y*nxsy+x;
            v   = avimage[k];
            if (v==0)
                denoisy[k]=0;
            else
            {
                //printf("v=%f\n",v);
                //compute the similar function
                sum=0;
                ar=0;
                for(i=0,xp=x-nw;xp<=x+nw;xp++)
                {
                    for(yp=y-nw;yp<=y+nw;yp++,i++)
                    {
                        adrp=yp*nxsy+xp;
                        for(j=npnp,dist2=0.,dd=dadr;j--;dd++)
                        {
                            dist  =  (double)synoisy[adr+*dd]-synoisy[adrp+*dd];
                            dist2 += (double)w[j]*dist*dist;
                        }
                        //rho[i]  = srho[i]=MAX((sqrt(dist2)-1.414*v),0);
                        dist2 =  MAX((sqrt(dist2)-1.414*v),0);
                        //wrho  =  (double)exp(-dist2/(0.2*v));
                        wrho  =  (double)exp(-dist2/(nlh*v));
                        sum   += (double)wrho;
                        ar  +=  (double)wrho*synoisy[adrp];
                        /*printf("wrho=%f\n",wrho);
                         * printf("synoisy[adrp]=%f\n",synoisy[adrp]);
                         * printf("ar=%f\n",ar);
                         * printf("sum=%f\n",sum);
                         **/
                    }
                }
                //printf("sum=%f\n",sum);
                dist =(double)ar/sum;
                //printf("dist=%f\n",dist);
                /*quickSort(srho,0,nwnw-1);
                 * ar=weightswidth(srho, v, nwnw);
                 *
                 * for(dist=0.,sum=0.,i=0,xp=x-nw;xp<=x+nw;xp++)
                 * {
                 * for(yp=y-nw;yp<=y+nw;yp++,i++)
                 * {
                 * dist2 = (double)(ar-rho[i]);
                 * wrho  =  MAX(dist2 ,0. );
                 * dist  += (double)wrho*synoisy[yp*nxsy+xp];
                 * sum   += (double)wrho;
                 * }
                 * }*/
                
                denoisy[k] = MIN( MAX( (dist),0 ),255 );
            }
        }
    }
    
    free(avimage);
    free(rho);
    free(srho);
    free(synoisy);
    free(w);
    free(dadr);
}


/**************************************************************************/
/*****************************************/
/**************The gateway routine.****  */
/*****************************************/

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /*input for C function */
    double  *in;         //input image
    double  *yest;       //output image
    double  sigma,nlh;       //variance
    int     nx,ny;       //parameters of images
    int     nwin, npat;  //parameters of windows
    if (nrhs < 1 )
        mexErrMsgTxt("At least one inputs required.");
    if (nrhs > 6 )
        mexErrMsgTxt("Too many inputs required.");
    if (nlhs != 1)
        mexErrMsgTxt("One output required.");
    
    /*Input Image */
    nx    = mxGetM(prhs[0]);
    ny    = mxGetN(prhs[0]);
    in    = mxGetPr(prhs[0]);
//sigma  = mxGetScalar(prhs[1]);
    
    if (nrhs<2)
        nwin=3;
    else
        nwin=mxGetScalar(prhs[1]);
    if (nrhs<3)
        npat=6;
    else
        npat=mxGetScalar(prhs[2]);
    if (nrhs<4)
        nlh=0.5;
    else
        nlh=mxGetScalar(prhs[3]);
    /*output image*/
    plhs[0]        = mxCreateDoubleMatrix(nx,ny, mxREAL);
    yest           = mxGetPr(plhs[0]);        /*output data*/
    
    owf_denoising(in, yest,sigma,nwin,npat,nx,ny,nlh);
}




