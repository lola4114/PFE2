// I = argmaxmin_mex(X, DIM, max_NOT_MIN)
// if max_NOT_MIN = 0 then argmin, else argmax

#include "mex.h"   

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{       
    double *X, *I; //X is the input vector/matrix and I is the output
    int DIM, max_NOT_MIN, X_rows, X_cols;
  
    mxArray *XData;
    
    //Copy input pointer X
    XData = prhs[0];
    
    //Get matrix X
    X = mxGetPr(XData);
    X_rows = mxGetM(XData);
    X_cols = mxGetN(XData);
    
    DIM = mxGetScalar(prhs[1]);
    max_NOT_MIN = mxGetScalar(prhs[2]);
    
    // mexPrintf("Row: %d, Col: %d, DIM: %d, max_NOT_MIN: %d\n",X_rows, X_cols, DIM, max_NOT_MIN);
    
    if (DIM==1){ 
        plhs[0] = mxCreateDoubleMatrix(1, X_cols, mxREAL); //mxReal is our data-type
    }else{
        if(DIM==2){
            plhs[0] = mxCreateDoubleMatrix(X_rows, 1, mxREAL); //mxReal is our data-type
        }else{
            mexPrintf("\nError: argmax for 3D-array is still unsupported.\n");
            return;
        }
    }
    I = mxGetPr(plhs[0]);

    // In order to maximize the speed, two specular functions are used: arg_min and arg_max
    if (max_NOT_MIN == 0)   
        arg_min(X,X_rows,X_cols,I,DIM);
    else
        arg_max(X,X_rows,X_cols,I,DIM);

 }
 
 void arg_max(double *X, int X_rows, int X_cols, double *I, int DIM, int max_NOT_MIN){
    int r, c, pmax;
    double vmax;
    // pmax: position of the maximum;   vmax: maximum value
    
    if (DIM==1){
        for (c=0; c<X_cols; c++){
            pmax = 1;
            vmax = X[0+c*X_rows];
            for (r=1; r<X_rows; r++){                
                if (X[r+c*X_rows] > vmax){
                    vmax = X[r+c*X_rows];
                    pmax = r+1;
                }
            }
            I[c] = pmax;
        }
    }else{
        if(DIM==2){
            for (r=0; r<X_rows; r++){
                pmax = 1;
                vmax = X[r+0*X_rows];
                for (c=1; c<X_cols; c++){
                    if (X[r+c*X_rows] > vmax){
                        vmax = X[r+c*X_rows];
                        pmax = c+1;
                    }
                }
                I[r] = pmax;
            }
        }
    }     
 }

  void arg_min(double *X, int X_rows, int X_cols, double *I, int DIM){
    int r, c, pmin;
    double vmin;
    // pmin: position of the minimum;   vmin: minimum value
    
    if (DIM==1){
        for (c=0; c<X_cols; c++){
            pmin = 1;
            vmin = X[0+c*X_rows];
            for (r=1; r<X_rows; r++){                
                if (X[r+c*X_rows] < vmin){
                    vmin = X[r+c*X_rows];
                    pmin = r+1;
                }
            }
            I[c] = pmin;
        }
    }else{
        if(DIM==2){
            for (r=0; r<X_rows; r++){
                pmin = 1;
                vmin = X[r+0*X_rows];
                for (c=1; c<X_cols; c++){
                    if (X[r+c*X_rows] < vmin){
                        vmin = X[r+c*X_rows];
                        pmin = c+1;
                    }
                }
                I[r] = pmin;
            }
        }
    }     
 }
 