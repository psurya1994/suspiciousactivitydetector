#include "mex.h"
#include "cv.h" 
#include "highgui.h"

// compiled with opencv2.0 gives different results then opencv2.1
// e.g. buffy_s5e2/buffy_s5e2_frames/001668.jpg

// matlab interface to OpenCV cvHaarDetectObjects function  
// author: Marcin Eichner


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char *classifiername;
    int classifiername_length,ncols, nrows,c,r,d;
    uchar *image_pointer; // equivalent to uint8
    uchar temp;
    const int flags = 0; // no CV_HAAR_DO_CANNY_PRUNING
    // init default values:
    double scale_factor = 1.1;
    int min_neighbors = 2,min_size_x = 30,min_size_y = 30;
    
    if(nrhs < 2 || nrhs > 7) 
      mexErrMsgTxt("ERROR: 2-7 arguments expected: <string xml model path>, <uint8 grayscale image>,\n"
              "[min_size_x], [min_size_y], [scale_factor], [min_neighbors]\n"
              "for optional arguments 3-7 check the cvHaarDetectObjects opencv function help");
    if (!mxIsChar(prhs[0]))
      mexErrMsgTxt("arg1 must be a string - absolute/relative path to the classifier");
    classifiername_length = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
    classifiername = mxArrayToString(prhs[0]);
    
    // HANDLE INPUT
    
    // read in the classifier
    CvHaarClassifierCascade * cascade = (CvHaarClassifierCascade*) cvLoad( classifiername, 0, 0, 0);
    if(!cascade)
       mexErrMsgTxt("ERROR: Could not load the classifier" );

    /* Check if the prhs image is in double format*/
    if (!(mxIsUint8(prhs[1]))) 
      mexErrMsgTxt("ERROR: arg2 (image) must be of type uint8");

    int image_ndim = mxGetNumberOfDimensions(prhs[1]);
    if (mxGetNumberOfDimensions(prhs[1]) > 2) 
      mexErrMsgTxt("ERROR: arg2 (image) must be a imgscale image");

    image_pointer =  (uchar*)mxGetPr(prhs[1]);
    ncols = mxGetN(prhs[1]); 
    nrows = mxGetM(prhs[1]); 
    
    IplImage* img = cvCreateImage( cvSize(ncols, nrows), IPL_DEPTH_8U, 1 );
    
    // Load the column wise vector into the IplImage
    // IplImage data is read in a rowwise manner
    // Appropriate conversion is carried out here
    for(c=0;c<ncols;c++)
       for(r=0;r<nrows;r++)
        {
            temp = (uchar)image_pointer[r+(nrows*c)];
            ((uchar *)(img->imageData + r*img->widthStep))[c]=temp;
        }
    
    if (nrhs > 2)
      min_size_x = (int)round(mxGetScalar(prhs[2]));
    if (nrhs > 3)
      min_size_y = (int)round(mxGetScalar(prhs[3]));
    if (nrhs > 4)
      scale_factor = mxGetScalar(prhs[4]);
    if (nrhs > 5)
      min_neighbors = (int)round(mxGetScalar(prhs[5]));
   
    CvMemStorage* storage = cvCreateMemStorage(0);
    
    // PROCESS IMAGE
      
    // equalize hist image before running face detection;
    cvEqualizeHist(img, img);

    //mexPrintf("%f %d %d %d %d",scale_factor,min_neighbors,flags,min_size_x,min_size_y);
    // run the classifier
    CvSeq* detections = cvHaarDetectObjects(img, cascade, storage, scale_factor, min_neighbors, flags,
                                            cvSize(min_size_x, min_size_y));

    if (detections->total == 0)
    {
      plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
      cvReleaseImage(&img);
      cvReleaseHaarClassifierCascade(&cascade);
      cvReleaseMemStorage(&storage);
      return;
    }
    else
    {
        plhs[0] = mxCreateDoubleMatrix(detections->total, 4, mxREAL);
        double* output = mxGetPr(plhs[0]);
        CvRect* rect;
        for(d = 0; d < detections->total; d++)
        {
           rect = (CvRect*)cvGetSeqElem(detections, d);
           output[d] = rect->x+1;
           output[d+detections->total*1] = rect->y+1;
           output[d+detections->total*2] = rect->width;
           output[d+detections->total*3] = rect->height;
        }
    }
    cvReleaseImage(&img);
    cvReleaseHaarClassifierCascade(&cascade);
    cvReleaseMemStorage(&storage);
//     cvDestroyAllWindows();
}

