/* Y=spktime2count(s,idx, Tw, Ncount,option)  */ 

/* transform spike trains to spike counts
   Y is neuron # x trial #  
   counts are non-overlapping 
   option=1, if neuronIdx is continuous and sorted, 0 if not 
   count from t=0; 
*/ 


#include "mex.h"
#include "math.h"
#include "time.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

int i,j,k, ID,count, temp, Nid, m1,m2,Ns, Ncount,ID1,ID2,option;
double *neuronIdx, *s;
double  Tw, T, *Pr,ts, *Y;

/******
 * Import variables from matlab
 * This is messy looking and is specific to mex.
 * Ignore if you're implementing this outside of mex.
 *******/
s = mxGetPr(prhs[0]);
m1 = mxGetM(prhs[0]);
Ns = mxGetN(prhs[0]);
if(m1!=2){
    mexErrMsgTxt("s should be 2xNs");
}

neuronIdx = mxGetPr(prhs[1]);
Nid = mxGetM(prhs[1]);
m1 = mxGetN(prhs[1]);
if(m1!=1){ 
    temp=Nid;
    Nid=m1; 
    m1=temp; 
}
if(m1!=1){
    mexErrMsgTxt("neuron id needs to be a vector or a row.");
}

Tw = mxGetScalar(prhs[2]);

Ncount = (int)mxGetScalar(prhs[3]);

option = (int)mxGetScalar(prhs[4]); /* 1, if neuronIdx is continuous and sorted, 0 if not */
if(option){  
  ID1=neuronIdx[0]; 
  ID2=neuronIdx[Nid-1];
/*  mexPrintf("ID1=%d,ID2=%d,option=%d",ID1,ID2,option); */ 
  
}

                                     
/* Allocate output vector */ 

/* mexPrintf("Ns=%d, Ncount=%d, Nid=%d",Ns,Ncount,Nid); 
 mexErrMsgTxt("stop"); */

plhs[0] = mxCreateDoubleMatrix(Nid, Ncount, mxREAL);
Y=mxGetPr(plhs[0]);

/*   main codes  */

for(i=0;i<Nid*Ncount;i++){
    Y[i]=0;
}

Pr=&s[0];
for(k=0;k<Ns;k++){
    ts=*Pr; 
    count=(int)floor(ts/Tw); 
     
    if(count<Ncount){
        Pr++;ID=(int)*Pr; 
        if(option){
            if(ID>ID1-0.1&ID<ID2+0.1){
                j=ID-ID1;
                Y[j+count*Nid]++;}
        }
        else{   
            for(j=0;j<Nid;j++){
                if(((int)neuronIdx[j])==ID){
                    
                    Y[j+count*Nid]++;}
            }
          } 
    }
    else{ 
        break;}
    Pr++;
}
    
}
    
    
    
    
