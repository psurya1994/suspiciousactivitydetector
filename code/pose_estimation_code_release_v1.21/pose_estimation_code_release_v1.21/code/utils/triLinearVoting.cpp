// function builds 3 dimensional histogram using values from Nx3 input array according to a given quantization. 
// It performs a soft-voting into 8 bins neighbouring a pixel in 3d color space according to the manhatan distance from a point to bin centers
// first input is a 3 channel (double) image expected as 2d array nPixels x 3 with values in range <0,255>
// second input is the quantization it can be a scalar (then uniform quantization of 3d color cube is done) or a vector [quan1 quan2 quan3] 
// third input (optional) is a k column vector of masks each nPixels x 1
// if it's not specified then 1's are assumed for every entry
// first output is the k-histogras as array (Dim x k) where Dim is quan(1)*quan(2)*quan(3) and k corresponding to the number of masks given
// if no mask was passed in then 1 histogram is returned
// second output is (Dim x 3) array with bin center associated with histogram
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
    mwSize nPixels,nImCh,quan_r,quan_c;
    double *img,*pi,*temp,*binC[NCHANNELS],*hist;
    mxArray *histogram;
    int dim;
    mwIndex quan[NCHANNELS];
    mwIndex lower[NCHANNELS];
    double sizbin[NCHANNELS] = {255, 255, 255};
    
    /* Check for proper number of input and output arguments */    
    if (!(nrhs == 2 || nrhs == 3))
      mexErrMsgTxt("2[3] input arguments required - image, quantization, [weights]");
    if(nlhs > 2)
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
      mexErrMsgTxt("Incorrect image no real part present");
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
    
    //handle weights matrix
    mwIndex p;
    double* weights;
    mwIndex * pp;
    mwIndex weightsnotgiven = 0;
    mwSize weights_c;
 
    if(nrhs == 3)
    {
      if (!(mxIsDouble(prhs[2])))
        mexErrMsgTxt("weights must be of type double.");
//      if(mxGetN(prhs[2]) != 1)
//        mexErrMsgTxt("Incorrect size of weights matrix - only 1 channel supported");
      if(mxGetM(prhs[2]) != nPixels)
        mexErrMsgTxt("Incorrect size of weights matrix - number of entries in each weight channel must match number of pixels in the image");
      weights = mxGetPr(prhs[2]);
      weights_c = mxGetN(prhs[2]);
      if(!weights || mxGetPi(prhs[2]))
        mexErrMsgTxt("Incorrect weights matrix");
      pp = &p;
    }
     else //weights matrix not given
    {
      weights = new double(1.0);
      pp = &weightsnotgiven;
      weights_c = 1;
    }
 
    
    
    // compute bin centers
    for(int i=0; i < NCHANNELS; i++)
    {
      sizbin[i] = sizbin[i]/quan[i];
      binC[i] = new double[quan[i]];
      for(int j=0; j<quan[i]; j++)
        binC[i][j] = sizbin[i]/2 + j*sizbin[i];
    }
    
    // allocate matrix for histogram
    dim = quan[0]*quan[1]*quan[2];
    histogram = mxCreateDoubleMatrix(dim,weights_c,mxREAL);
    if (!histogram)
      mexErrMsgTxt("Not enough memory to allocate histogram matrix");
    hist = mxGetPr(histogram);
    
    // zero the histogram!!!
    memset(hist, 0, dim*sizeof(double)); //much faster!
    
    // loop through all pixels and built a histogram using soft voting with manhattan distance
    for(p=0; p < nPixels; p++)
    {
      for(mwIndex c=0; c< NCHANNELS; c++)
        lower[c] = mwIndex(floor((img[p+c*nPixels]-binC[c][0])/sizbin[c]));
      
      // calculating seaprate histogram for every column mask provided
      for(mwIndex wc=0; wc < weights_c; wc++)
      {
        
        if (0 <= lower[0]) {
          if (0 <= lower[1]) {
            if (0 <= lower[2]) 
              hist[lower[0] + lower[1]*quan[0] + lower[2]*quan[1]*quan[0] + wc*dim]            +=weights[*pp+wc*nPixels]*  (1-fabs(img[p]-binC[0][lower[0]])/sizbin[0])*  (1-fabs(img[p+1*nPixels]-binC[1][lower[1]])/sizbin[1])*  (1-fabs(img[p+2*nPixels]-binC[2][lower[2]])/sizbin[2]);
            if (lower[2] < quan[2]-1) 
              hist[lower[0] + lower[1]*quan[0] + (lower[2]+1)*quan[1]*quan[0] + wc*dim]        +=weights[*pp+wc*nPixels]*  (1-fabs(img[p]-binC[0][lower[0]])/sizbin[0])*  (1-fabs(img[p+1*nPixels]-binC[1][lower[1]])/sizbin[1])*(1-fabs(img[p+2*nPixels]-binC[2][lower[2]+1])/sizbin[2]);
          }
          if (lower[1] < quan[1]-1) {
            if (0 <= lower[2]) 
              hist[lower[0] + (lower[1]+1)*quan[0] + lower[2]*quan[1]*quan[0] + wc*dim]        +=weights[*pp+wc*nPixels]*  (1-fabs(img[p]-binC[0][lower[0]])/sizbin[0])*(1-fabs(img[p+1*nPixels]-binC[1][lower[1]+1])/sizbin[1])*  (1-fabs(img[p+2*nPixels]-binC[2][lower[2]])/sizbin[2]);
            if (lower[2] < quan[2]-1) 
              hist[lower[0] + (lower[1]+1)*quan[0] + (lower[2]+1)*quan[1]*quan[0] + wc*dim]    +=weights[*pp+wc*nPixels]*  (1-fabs(img[p]-binC[0][lower[0]])/sizbin[0])*(1-fabs(img[p+1*nPixels]-binC[1][lower[1]+1])/sizbin[1])*(1-fabs(img[p+2*nPixels]-binC[2][lower[2]+1])/sizbin[2]);
          }
        }
        if (lower[0] < quan[0]-1) {
          if (0 <= lower[1]) {
            if (0 <= lower[2])  
              hist[1+lower[0] + lower[1]*quan[0] + lower[2]*quan[1]*quan[0] + wc*dim]          +=weights[*pp+wc*nPixels]*(1-fabs(img[p]-binC[0][lower[0]+1])/sizbin[0])*  (1-fabs(img[p+1*nPixels]-binC[1][lower[1]])/sizbin[1])*  (1-fabs(img[p+2*nPixels]-binC[2][lower[2]])/sizbin[2]);
            if (lower[2] < quan[2]-1) 
              hist[1+lower[0] + lower[1]*quan[0] + (lower[2]+1)*quan[1]*quan[0] + wc*dim]      +=weights[*pp+wc*nPixels]*(1-fabs(img[p]-binC[0][lower[0]+1])/sizbin[0])*  (1-fabs(img[p+1*nPixels]-binC[1][lower[1]])/sizbin[1])*(1-fabs(img[p+2*nPixels]-binC[2][lower[2]+1])/sizbin[2]);
          }
          if (lower[1] < quan[1]-1) {
            if (0 <= lower[2]) 
              hist[1+lower[0] + (lower[1]+1)*quan[0] + lower[2]*quan[1]*quan[0] + wc*dim]      +=weights[*pp+wc*nPixels]*(1-fabs(img[p]-binC[0][lower[0]+1])/sizbin[0])*(1-fabs(img[p+1*nPixels]-binC[1][lower[1]+1])/sizbin[1])*  (1-fabs(img[p+2*nPixels]-binC[2][lower[2]])/sizbin[2]);
            if (lower[2] < quan[2]-1) 
              hist[1+lower[0] + (lower[1]+1)*quan[0] + (lower[2]+1)*quan[1]*quan[0] + wc*dim]  +=weights[*pp+wc*nPixels]*(1-fabs(img[p]-binC[0][lower[0]+1])/sizbin[0])*(1-fabs(img[p+1*nPixels]-binC[1][lower[1]+1])/sizbin[1])*(1-fabs(img[p+2*nPixels]-binC[2][lower[2]+1])/sizbin[2]);
          }
        }

      }
    }
    
  // assign created histogram to first output
  plhs[0] = histogram;   
    
  // if 2 output arguments required 
  // calculate list of bincenters in 3 dimensional cube corresponding to histogram bins and return it as second output
  mxArray* binCenters;  
  if(nlhs > 1)  
  {
    binCenters = mxCreateDoubleMatrix(dim,NCHANNELS,mxREAL);
    double* temp = mxGetPr(binCenters);
    
    mwIndex index = 0;
    for(int b = 0; b < quan[2]; b++)
      for(int g = 0; g < quan[1]; g++)
        for(int r = 0; r < quan[0]; r++)
        {
        temp[index] = binC[0][r % quan[0]];
        temp[dim+index] = binC[1][g % quan[1]];
        temp[2*dim+index] = binC[2][b % quan[2]];
        index++;
        }
    plhs[1] = binCenters;
  }

  //clean up
  for(int i=0; i < NCHANNELS; i++)
    delete [] binC[i];
  
   if(nrhs < 3)
     delete weights;
}
 

