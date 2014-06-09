// function builds 2 dimensional softsegmentation image (nPixes x 1) interpolating values from a given 3d histogram and an input image
// first input parameter is a 3 channel (double) image expected as 2d array nPixels x 3 with values in range <0,255>
// second input parameter is a quantization used to build histogram - it is required to reach istogram bins in 3D 
// third input parameter is a 3d histogram, send as a row vector dim x 1,  
// product of the quanization (dim) must be equal to number of rows in the histogram array (3rd param)
// more than one histogram may be sent at the time (hist is then dim x classes), if so then the output image will also be multidimensional (nPixels x classes)
// this routine is used to interpolate a posterior probability of pixel form a soft-voted histogram 
#include "mex.h"
#include "matrix.h" // mxArray definition
#include <math.h>  // floor()
#include <string.h> // memset()

#if defined(NAN_EQUALS_ZERO)
#define IsNonZero(d) ((d)!=0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d)!=0.0)
#endif

void mexFunction(
		 int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]
		 )
{
    /* Declare variable */
    const int NCHANNELS = 3;
    mwSize nPixels,nImCh,quan_r,quan_c,classes;
    double *img,*pi,*temp,*binC[NCHANNELS],*out;
    mxArray *output;
    int dim;
    mwIndex quan[NCHANNELS];
    mwIndex lower[NCHANNELS];
    double sizbin[NCHANNELS] = {255, 255, 255};
    mwIndex p,clas;
    
    /* Check for proper number of input and output arguments */    
    if (!(nrhs == 3))
      mexErrMsgTxt("3 input arguments required - image, quantization, posterior");
    if(nlhs > 1)
      mexErrMsgTxt("Too many output arguments.");
    /* Check data type of input argument  */
    if (!(mxIsDouble(prhs[0])))
      mexErrMsgTxt("Image must be of type double.");

    // handle image data
    nPixels  = mxGetM(prhs[0]);
    nImCh  = mxGetN(prhs[0]);
    img = mxGetPr(prhs[0]);
    pi = mxGetPi(prhs[0]);
        
    if (pi != NULL || img == NULL)
      mexErrMsgTxt("Incorrect image no real part or imaginary part present");
    if (nImCh != 3)
      mexErrMsgTxt("Image second dimension != 3 - only 3 channel images supported");
    
    // handle quantization matrix
    quan_r  = mxGetM(prhs[1]);
    quan_c  = mxGetN(prhs[1]);
    temp = mxGetPr(prhs[1]);
    pi = mxGetPi(prhs[1]);
    
    if (!(mxIsDouble(prhs[1])))
        mexErrMsgTxt("quantiz must be of type double.");
    if (pi != NULL || temp == NULL)
      mexErrMsgTxt("Incorrect quantization matrix");
    if (quan_r != 1)
     	mexErrMsgTxt("quantization matrix must be 1x3 matrix or scalar");
    if (quan_c == NCHANNELS) // matrix
    {
      quan[0] = int(temp[0]);
      quan[1] = int(temp[1]);
      quan[2] = int(temp[2]);
    }
    else if (quan_c == 1) // scalar
    {
      quan[0] = int(temp[0]);
      quan[1] = int(temp[0]);
      quan[2] = int(temp[0]);
    }
    
    
    dim = quan[0]*quan[1]*quan[2];
    
    //handle histogram matrix
    double* hist;

    classes = mxGetN(prhs[2]);
    if(!(mxIsDouble(prhs[2])))
        mexErrMsgTxt("posterior must be of type double.");
    if(mxGetM(prhs[2]) != dim) 
      mexErrMsgTxt("Number of rows in posterior must match product of quantization ");
    
    hist = mxGetPr(prhs[2]);
    if(!hist || mxGetPi(prhs[2]))
      mexErrMsgTxt("Incorrect posterior matrix");
    
    // compute bin centers
    for(int i=0; i < NCHANNELS; i++)
    {
      sizbin[i] = sizbin[i]/quan[i];
      binC[i] = new double[quan[i]];
      for(int j=0; j<quan[i]; j++)
        binC[i][j] = sizbin[i]/2 + j*sizbin[i];
    }
    
    // allocate matrix for histogram
    output = mxCreateDoubleMatrix(nPixels,classes,mxREAL);
    if (!output)
      mexErrMsgTxt("Not enough memory to allocate histogram matrix");
    out = mxGetPr(output);
    
    // clear the output image!!!
    memset(out, 0, nPixels*classes*sizeof(double)); //much faster!
    
    // loop through all pixels and calculate output image interpolating values of 8 neighbouring histogram bins 
    for(p=0; p < nPixels; p++)
    {
      for(int c=0; c< NCHANNELS; c++)
        lower[c] = mwIndex(floor((img[p+c*nPixels]-binC[c][0])/sizbin[c]));
      
      for(clas=0; clas < classes; clas++)
      {

        if (0 <= lower[0]) {
          if (0 <= lower[1]) {
            if (0 <= lower[2]) 
              out[p+clas*nPixels] += hist[lower[0] + lower[1]*quan[0] + lower[2]*quan[1]*quan[0] + clas*dim]            *  (1-fabs(img[p]-binC[0][lower[0]])/sizbin[0])*  (1-fabs(img[p+1*nPixels]-binC[1][lower[1]])/sizbin[1])*  (1-fabs(img[p+2*nPixels]-binC[2][lower[2]])/sizbin[2]);
            if (lower[2] < quan[2]-1) 
              out[p+clas*nPixels] += hist[lower[0] + lower[1]*quan[0] + (lower[2]+1)*quan[1]*quan[0] + clas*dim]        *  (1-fabs(img[p]-binC[0][lower[0]])/sizbin[0])*  (1-fabs(img[p+1*nPixels]-binC[1][lower[1]])/sizbin[1])*(1-fabs(img[p+2*nPixels]-binC[2][lower[2]+1])/sizbin[2]);
          }
          if (lower[1] < quan[1]-1) {
            if (0 <= lower[2]) 
              out[p+clas*nPixels] += hist[lower[0] + (lower[1]+1)*quan[0] + lower[2]*quan[1]*quan[0] + clas*dim]        *  (1-fabs(img[p]-binC[0][lower[0]])/sizbin[0])*(1-fabs(img[p+1*nPixels]-binC[1][lower[1]+1])/sizbin[1])*  (1-fabs(img[p+2*nPixels]-binC[2][lower[2]])/sizbin[2]);
            if (lower[2] < quan[2]-1) 
              out[p+clas*nPixels] += hist[lower[0] + (lower[1]+1)*quan[0] + (lower[2]+1)*quan[1]*quan[0] + clas*dim]    *  (1-fabs(img[p]-binC[0][lower[0]])/sizbin[0])*(1-fabs(img[p+1*nPixels]-binC[1][lower[1]+1])/sizbin[1])*(1-fabs(img[p+2*nPixels]-binC[2][lower[2]+1])/sizbin[2]);
          }
        }
        if (lower[0] < quan[0]-1) {
          if (0 <= lower[1]) {
            if (0 <= lower[2])  
              out[p+clas*nPixels] += hist[1+lower[0] + lower[1]*quan[0] + lower[2]*quan[1]*quan[0] + clas*dim]          *(1-fabs(img[p]-binC[0][lower[0]+1])/sizbin[0])*  (1-fabs(img[p+1*nPixels]-binC[1][lower[1]])/sizbin[1])*  (1-fabs(img[p+2*nPixels]-binC[2][lower[2]])/sizbin[2]);
            if (lower[2] < quan[2]-1) 
              out[p+clas*nPixels] += hist[1+lower[0] + lower[1]*quan[0] + (lower[2]+1)*quan[1]*quan[0] + clas*dim]      *(1-fabs(img[p]-binC[0][lower[0]+1])/sizbin[0])*  (1-fabs(img[p+1*nPixels]-binC[1][lower[1]])/sizbin[1])*(1-fabs(img[p+2*nPixels]-binC[2][lower[2]+1])/sizbin[2]);
          }
          if (lower[1] < quan[1]-1) {
            if (0 <= lower[2]) 
              out[p+clas*nPixels] += hist[1+lower[0] + (lower[1]+1)*quan[0] + lower[2]*quan[1]*quan[0] + clas*dim]      *(1-fabs(img[p]-binC[0][lower[0]+1])/sizbin[0])*(1-fabs(img[p+1*nPixels]-binC[1][lower[1]+1])/sizbin[1])*  (1-fabs(img[p+2*nPixels]-binC[2][lower[2]])/sizbin[2]);
            if (lower[2] < quan[2]-1) 
              out[p+clas*nPixels] += hist[1+lower[0] + (lower[1]+1)*quan[0] + (lower[2]+1)*quan[1]*quan[0] + clas*dim]  *(1-fabs(img[p]-binC[0][lower[0]+1])/sizbin[0])*(1-fabs(img[p+1*nPixels]-binC[1][lower[1]+1])/sizbin[1])*(1-fabs(img[p+2*nPixels]-binC[2][lower[2]+1])/sizbin[2]);
          }
        }
      }

    }
    
  // assign created histogram to first output
  plhs[0] = output;   
    
  //clean up
  for(int i=0; i < NCHANNELS; i++)
    delete [] binC[i];
  
}
 

