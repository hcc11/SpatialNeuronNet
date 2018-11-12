/* To use function from matlab, first compile by entering this into the Matlab command window:
   mex EIF1DRF2tausyn.c
   Then call the function like this:
   [s,Isynrecord,vr]=EIF1DRFfastslowSyn(sx, Wrf,Wrr,param);
 */

/* sx is the feedforward presynaptic spike trains, which should be 3xNsx where Nsx is the number of spikes.
      sx(1,:) is spike times, sx(2,:) is the index of the neuron (from 1 to Nx)
   Wrr is a vector of connections among the recurrent layer, containing postsynaptic cell indices, 
      sorted by the index of the presynaptic cell. The block of postsynaptic cell indices for each presynaptic
      cell is sorted as excitatory followed by inhibitory cells. I use fixed number of projections Kab to each population. 
      For example, Wrr[j*(Kee+Kie)] to Wrr[j*{Kee+Kie)+Kee-1] are connections from j to E pop and 
      Wrr[j*(Kee+Kie)+Kee] to Wrr[(j+1)*{Kee+Kie)-1] are connections from j to I pop. 
   Wrf is a vector of connections from the feedforward layer to the recurrent layer, sorted by the index of the presynaptic cell.
      The block of postsynaptic cell indices for each presynaptic cell is sorted as excitatory followed by inhibitory cells.
   
   param is a struc w/ fields: Ne, Ni, Nx, Jx, Jr, Kx, Kr, 
      gl, Cm, vlb, vth, DeltaT, vT, vl, vre, tref, tausyn, V0, T, dt,
      maxns, Irecord, Psyn
  Jx=[Jex; Jix]; Jr=[Jee, Jei; Jie, Jii];
  Kx=[Kex; Kix]; Kr=[Kee, Kei; Kie, Kii]; 
  taursyn: syn rise time const, 3x(Nsyntype), rows: X, E, I; cols: syn type
  taudsyn: syn decay time const, 3x(Nsyntype), rows: X, E, I; cols: syn type
  Psyn(i,j): percentage of synapse j for (X, E, I (i=1,2,3))
 
  
   Kab is the number of projections from each cell in pop b=e,i,x to all cells in pop a=e,i
   Iapp is the constant external input.  It should be a vector of size Nx1 or 1xN where N is the number of cells in the network. 
   Ne, Ni are the number of excitatory, inhibitory cells, N=Ne+Ni.
   Nx is the number of neurons in feedforward layer.  
   Jab is the synaptic strength of connections from b=e,i,x to a=e,i.
   
   C,gl,vl,DeltaT,VT,tref,Vth,Vre,Vlb are EIF neuron params
      They are each 2x1 vectors, for exc and inh neurons separately.
      For example, C(1) is the capacitance of exc neurons and C(2) of inh neurons 
   tausynb is the time-constant of the synapses from neurons in population b=x,e,i.
     post-synaptic currents are of the form (1/tausynb)*exp(-t/tausynb) where t>0 is 
     the time evolved since the presynaptic spike.
   V0 is vector of all membrane potential initial conditions.
      It should be Nx1 where N=Ne+Ni is the number of neurons in the recurrent network
      The first Ne elements are for exc neurons, the last Ni for inh neurons
   dt is bin size for time
   maxns is maximum number of spikes allowed for all neurons together
   Irecord is a 1x(Nrecord) matrix indicating for which neurons we should record
     the synaptic inputs and membrane potential. The first Ne elements are for exc neurons, 
     the last Ni for inh neurons
 
 
   Outputs:
   s is a 2x(maxns) matrix of spikes
      s(1,:) contains spike times.
      s(2,:) contains indices of neurons that spike
      When there are fewer than maxns spikes, extra space in s will be filled 
      with zeros.  This should be truncated in matlab by writing s=s(:,s(1,:)>0);
   Isynrecord,vr are the recorded synaptic inputs and voltage respectively.
 
 */



#include "mex.h"
#include "math.h"
#include "time.h"
#include "matrix.h"


/* A fast approximation of the exp function */
static union 
{
	  double d;
	    struct {
#ifdef LITTLE_ENDIAN
	    int j,i;
#else 
		    int i,j;
#endif
	  } n;
} _eco;
#define EXP_A (1048576/0.69314718055994530942)
#define EXP_C 60801
#define EXP(y) (_eco.n.i = EXP_A*(y) + (1072693248 - EXP_C), _eco.d)


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

int Ntref[2],Ni,Ne,Nx,Kee,Kie,Kei,Kii,Kex,Kix,k,j,i,N,Nt,m1,m2,maxns,ns,Nrecord,jj,nskiprecord,tempflag,Nsx, Nw;
double dt,*s,*v,*v0,*Isynrecord, *Psyn;
double *Isyn, *Isynprime;
int *Wrr, *Wrf, *Pr; 
double *taudsyn,*taursyn,*vr,*sx,Jex,Jix,*iXspkInd;
double Jee,Jei,Jie,Jii,T,*Irecord,*C,*Vleak,*DeltaT,*VT,*tref,*gl,*Vth,*Vre,*Vlb,xloc,yloc;
int *refstate,iXspike,jspike,*postcellX,*postcellE,*postcellI,postcell, Nsyn, isyn, Nsyntype,Ke,Ki,Kx, *syntype; 
const mxArray  *mxTmp;
double  *temp1,*temp2,  Isyntot;



/******
 * Import variables from matlab
 * This is messy looking and is specific to mex.
 * Ignore if you're implementing this outside of mex.
 *******/
sx = mxGetPr(prhs[0]);
m1 = mxGetM(prhs[0]);
Nsx = mxGetN(prhs[0]);
if(m1!=2){
    mexErrMsgTxt("sx should be Nsxx2");
}

Wrf = (int *)mxGetData(prhs[1]);
Nw = mxGetM(prhs[1]);
m1 = mxGetN(prhs[1]);
if(m1!=1){
    mexErrMsgTxt("Weight matrix Wrf must be Nwx1 where Nw is numer of connections.");
}

Wrr = (int *)mxGetData(prhs[2]);
Nw = mxGetM(prhs[2]);
m1 = mxGetN(prhs[2]);
if(m1!=1){
    mexErrMsgTxt("Weight matrix Wrr must be Nwx1 where Nw is numer of connections.");
}

/* check if prhs[3] is a struct */ 
if(mxIsStruct(prhs[3])!=1){ 
    mexErrMsgTxt("The 4th input should be the parameter struct.");
}
  
  
/* Number of exc neurons in each direction. */ 
mxTmp = mxGetField(prhs[3],0,"Ne");
Ne=(int)mxGetPr(mxTmp)[0];

/* Number of inh neurons in each direction. */ 
mxTmp = mxGetField(prhs[3],0,"Ni");
Ni=(int)mxGetPr(mxTmp)[0];

/* Number of neurons in the ffwd layer in each direction. */ 
mxTmp = mxGetField(prhs[3],0,"Nx");
Nx=(int)mxGetPr(mxTmp)[0];

mxTmp = mxGetField(prhs[3],0,"Jx");
Jex=mxGetPr(mxTmp)[0];
Jix=mxGetPr(mxTmp)[1];

mxTmp = mxGetField(prhs[3],0,"Jr");
Jee= mxGetPr(mxTmp)[0];
Jie= mxGetPr(mxTmp)[1];
Jei= mxGetPr(mxTmp)[2];
Jii= mxGetPr(mxTmp)[3];

mxTmp = mxGetField(prhs[3],0,"Kx");
Kex= (int)mxGetPr(mxTmp)[0];
Kix= (int)mxGetPr(mxTmp)[1];

mxTmp = mxGetField(prhs[3],0,"Kr");
Kee= (int)mxGetPr(mxTmp)[0];
Kie= (int)mxGetPr(mxTmp)[1];
Kei= (int)mxGetPr(mxTmp)[2];
Kii= (int)mxGetPr(mxTmp)[3];

mxTmp = mxGetField(prhs[3],0,"Cm");
C=mxGetPr(mxTmp);
m1 = mxGetN(mxTmp);
m2 = mxGetM(mxTmp);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
mxTmp = mxGetField(prhs[3],0,"gl");
gl=mxGetPr(mxTmp);
m1 = mxGetN(mxTmp);
m2 = mxGetM(mxTmp);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
mxTmp = mxGetField(prhs[3],0,"vl");
Vleak=mxGetPr(mxTmp);
m1 = mxGetN(mxTmp);
m2 = mxGetM(mxTmp);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
        
mxTmp = mxGetField(prhs[3],0,"DeltaT");        
DeltaT=mxGetPr(mxTmp);
m1 = mxGetN(mxTmp);
m2 = mxGetM(mxTmp);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");

mxTmp = mxGetField(prhs[3],0,"vT");
VT=mxGetPr(mxTmp);
m1 = mxGetN(mxTmp);
m2 = mxGetM(mxTmp);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
mxTmp = mxGetField(prhs[3],0,"tref");
tref=mxGetPr(mxTmp);
m1 = mxGetN(mxTmp);
m2 = mxGetM(mxTmp);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
mxTmp = mxGetField(prhs[3],0,"vth");
Vth=mxGetPr(mxTmp);
m1 = mxGetN(mxTmp);
m2 = mxGetM(mxTmp);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
mxTmp = mxGetField(prhs[3],0,"vre");
Vre=mxGetPr(mxTmp);
m1 = mxGetN(mxTmp);
m2 = mxGetM(mxTmp);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");
mxTmp = mxGetField(prhs[3],0,"vlb");
Vlb=mxGetPr(mxTmp);
m1 = mxGetN(mxTmp);
m2 = mxGetM(mxTmp);
if(m1*m2!=2)
    mexErrMsgTxt("All neuron parameters should be 2x1");

mxTmp = mxGetField(prhs[3],0,"taursyn");
m1=mxGetN(mxTmp);
m2=mxGetM(mxTmp);
if(m2!=3)
    mexErrMsgTxt("size(taursyn,1) should be 3");
Nsyntype=m1;
taursyn=mxGetPr(mxTmp);

mxTmp = mxGetField(prhs[3],0,"taudsyn");
m1=mxGetN(mxTmp);
m2=mxGetM(mxTmp);
if(m2!=3)
    mexErrMsgTxt("size(taudsyn,1) should be 3");
if(m1!=Nsyntype)
    mexErrMsgTxt("size(taursyn,1) should equal size(taudsyn,2)");
taudsyn=mxGetPr(mxTmp);

mxTmp = mxGetField(prhs[3],0,"Psyn");
Psyn=mxGetPr(mxTmp);
m2=mxGetM(mxTmp);
if(m2!=3)
    mexErrMsgTxt("size(Psyn,1) should be 3");
m1=mxGetN(mxTmp);
if(m1!=Nsyntype)
    mexErrMsgTxt("size(Psyn,2) should equal size(taursyn,2)");

syntype=mxMalloc(m1*m2*sizeof(int));
Nsyn = 0;
for(isyn=0;isyn<m1*m2; isyn++){
    if(Psyn[isyn]){   
        syntype[Nsyn]=isyn%3;
        Nsyn++; }  /* type 0: X, 1:E, 2:I, for updating postsyn input type */ 
}

mxTmp = mxGetField(prhs[3],0,"V0");
v0 = mxGetPr(mxTmp);
N = mxGetM(mxTmp);
m2 = mxGetN(mxTmp);
if(N==1 && m2!=1)
    N=m2;

mxTmp = mxGetField(prhs[3],0,"T");
T = mxGetPr(mxTmp)[0];
mxTmp = mxGetField(prhs[3],0,"dt");
dt =mxGetPr(mxTmp)[0];

mxTmp = mxGetField(prhs[3],0,"maxns");
maxns = (int)mxGetPr(mxTmp)[0];

mxTmp = mxGetField(prhs[3],0,"Irecord");
Irecord=mxGetPr(mxTmp);
Nrecord = mxGetN(mxTmp);
m2 = mxGetM(mxTmp);
if(m2!=1)
    mexErrMsgTxt("Irecord should be Nx1.");

/******
 * Finished importing variables.
 *******/

/* Check for consistency with total number of neurons */
if(N!=Ne+Ni)
    mexErrMsgTxt("Ne1 and/or Ni1 not consistent with size of V0");

/* Numebr of time bins */
Nt=(int)(T/dt);

/* mexPrintf("Nsyn=%d, N=%d, Nrecord=%d, Nt=%d \n",Nsyn, N,Nrecord,Nt); */

/******
 * Now allocate new variables.
 * This is also mex specific.  Use malloc in C, etc.
 *****/

/* Allocate output vector */
plhs[0] = mxCreateDoubleMatrix(2, maxns, mxREAL);
s=mxGetPr(plhs[0]);

plhs[1] = mxCreateDoubleMatrix(Nrecord*Nsyn, Nt, mxREAL);
Isynrecord=mxGetPr(plhs[1]);  

plhs[2] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
vr=mxGetPr(plhs[2]);

/* Allocate membrane potential */
v = mxMalloc(N*sizeof(double));;
refstate=mxMalloc(N*sizeof(int));

Isyn = mxMalloc(Nsyn*N*sizeof(double));  /* synp. currents, NxNsyn, Nsyn: nnz of Psyn, col i corresponds to syntype[i]   */ 
Isynprime=mxMalloc(Nsyn*N*sizeof(double));

temp1=mxMalloc(Nsyn*sizeof(double)); /* temporary constant */
temp2=mxMalloc(Nsyn*sizeof(double)); 
for (isyn=0;isyn<Nsyn;isyn++){
     temp1[isyn]=(1/taudsyn[isyn]+1/taursyn[isyn]);
     temp2[isyn]=1/(taudsyn[isyn]*taursyn[isyn]);}

Kx=Kex+Kix;
Ke=Kee+Kie;
Ki=Kei+Kii;

postcellX=mxMalloc(Kx*sizeof(int)); /* index for postsynaptic cells */ 
postcellE=mxMalloc(Ke*sizeof(int));
postcellI=mxMalloc(Ki*sizeof(int));

/*****
 * Finished allocating variables
 ****/

/* Inititalize variables */
for(j=0;j<N;j++){
    v[j]=v0[j]; 
    refstate[j]=0;}

for(jj=0;jj<Nsyn*N;jj++){
    Isyn[jj]=0;
    Isynprime[jj]=0;}

/* Record first time bin */
for(jj=0;jj<Nrecord;jj++){
  if(Irecord[jj]<1 || Irecord[jj]>N+1)
     mexErrMsgTxt("Indices in Irecord must be between 1 and N");
    for(isyn=0;isyn<Nsyn;isyn++){
        Isynrecord[isyn+jj*Nrecord]=Isyn[(int)round(Irecord[jj]-1)*Nsyn+isyn];
    }
  vr[jj]=v[(int)round(Irecord[jj]-1)]; }

/* Refractory states */
Ntref[0]=(int)round(tref[0]/dt);
Ntref[1]=(int)round(tref[1]/dt);


/* mexPrintf("Jex=%.2f, Jix=%.2f\n Jee=%.2f, Jie=%.2f\n Jei=%.2f, Jii=%.2f\n", Jex,Jix,Jee,Jie,Jei,Jii); */
             
/* Initialize number of spikes */
ns=0;

/* Time loop */
/* Exit loop and issue a warning if max number of spikes is exceeded */
iXspike=0;
iXspkInd=&sx[0]; /* even index: spk times. odd index: neuron ID */ 
for(i=1;i<Nt && ns<maxns;i++){  
    /* Update synaptic variables */ 
    for(jj=0;jj<N*Nsyn;jj++){  
          isyn=jj%Nsyn; 
          Isyn[jj]+=Isynprime[jj]*dt;   
          Isynprime[jj]-=dt*(Isynprime[jj]*temp1[isyn]+Isyn[jj]*temp2[isyn]);
         
     }
    
     /* Find all spikes in feedforwar layer at this time bin */
     /* Add to corresponding elements of JnextX */
     while(*iXspkInd<=i*dt && iXspike<Nsx){
         iXspkInd++; /* point to neuron ID */
         jspike=(int)round((*iXspkInd)-1);
         if(jspike<0 || jspike>=Nx){
             mexPrintf("\n %d %d %d %d %d\n",(int)round(*(iXspkInd-1)/dt),iXspike,i,jspike,(int)round((*iXspkInd)-1));
             mexErrMsgTxt("Out of bounds index in sx.");
         }
         Pr=&Wrf[jspike*Kx]; 
         /*
         for(k=0;k<Kx;k++){ 
             postcell=(int)((*Pr)-1);
             if(postcell<0 || postcell>=N){
                 mexPrintf("\n Wrf j=%d, postcell=%d\n",j, postcell);
                 mexErrMsgTxt("postcell out of bounds");
             }
             if (postcell<Ne){
                 for(isyn=0;isyn<Nsyn;isyn++){
                     if(syntype[isyn]==0) 
                         Isynprime[postcell*Nsyn+isyn]+=Jex*temp2[isyn];}}
             else{
                  for(isyn=0;isyn<Nsyn;isyn++){
                     if(syntype[isyn]==0) 
                         Isynprime[postcell*Nsyn+isyn]+=Jix*temp2[isyn];}
            }
             Pr++;
         } */
         
         for(k=0;k<Kx;k++){ 
             postcellX[k]=(int)((*Pr)-1);
             if(postcellX[k]<0 || postcellX[k]>=N){
                 mexPrintf("\n Wrf j=%d, postcell=%d\n",j, postcellX[k]);
                 mexErrMsgTxt("postcell out of bounds");
             }
             Pr++;
         }
         
         for(isyn=0;isyn<Nsyn;isyn++){
             if(syntype[isyn]==0){
                 for(k=0;k<Kx;k++){
                     if (postcellX[k]<Ne)
                         Isynprime[postcellX[k]*Nsyn+isyn]+=Jex*temp2[isyn];
                     else
                         Isynprime[postcellX[k]*Nsyn+isyn]+=Jix*temp2[isyn];}
             }
         }
         

         iXspike++; 
         iXspkInd++; 
     } 
     

    
    
    /* loop over neurons  */ 
    Pr=&Wrr[0];
    for(j=0;j<N;j++){      
      /* Update membrane potential */
      /* Spikes will be propagated at the END of the time bin (see below)*/
        
        if(j<Ne){

             if(refstate[j]<=0){
                Isyntot=0;
                for(isyn=0;isyn<Nsyn;isyn++){
                Isyntot+=Psyn[isyn]*Isyn[j*Nsyn+isyn];
                }
                v[j]+=fmax((Isyntot-gl[0]*(v[j]-Vleak[0])+gl[0]*DeltaT[0]*EXP((v[j]-VT[0])/DeltaT[0]))*dt/C[0],Vlb[0]-v[j]);}
             else{                 
                if(refstate[j]>1)
                   v[j]=Vth[0];
                else
                   v[j]=Vre[0];
                refstate[j]--;
             }

              /* If a spike occurs */
              if(v[j]>=Vth[0] && refstate[j]<=0 && ns<maxns){
                  
                  refstate[j]=Ntref[0];
                  v[j]=Vth[0];       /* reset membrane potential */
                  s[0+2*ns]=i*dt; /* spike time */
                  s[1+2*ns]=j+1;     /* neuron index 1 */
                  ns++;           /* update total number of spikes */
                  /* For each postsynaptic target, update synaptic inputs */
                 /* Pr=&Wrr[j*Ke]; */
                  /*
                  for(k=0;k<Ke;k++){
                      postcell=(int)((*Pr)-1);
                      if(postcell<0 || postcell>=N){
                          mexPrintf("\n Wrf j=%d, postcell=%d\n",j, postcell);
                          mexErrMsgTxt("postcell out of bounds");
                      }
                      if (postcell<Ne){
                          for(isyn=0;isyn<Nsyn;isyn++){
                              if(syntype[isyn]==1)
                                  Isynprime[postcell*Nsyn+isyn]+=Jee*temp2[isyn];}}
                          else{
                              for(isyn=0;isyn<Nsyn;isyn++){
                                  if(syntype[isyn]==1)
                                      Isynprime[postcell*Nsyn+isyn]+=Jie*temp2[isyn];}
                      }
                      Pr++;
                  } */
                  
                  for(k=0;k<Ke;k++){
                      postcellE[k]=(int)(*Pr)-1;
                      if(postcellE[k]<0 || postcellE[k]>=N){
                         mexPrintf("\n exc j=%d, postcell=%d\n",j, postcellE[k]);
                         mexErrMsgTxt("postcell out of bounds");}
                      Pr++;
                  }
                  for(isyn=0;isyn<Nsyn;isyn++){
                      if(syntype[isyn]==1){
                          for(k=0;k<Ke;k++){
                              if (postcellE[k]<Ne)
                                  Isynprime[postcellE[k]*Nsyn+isyn]+=Jee*temp2[isyn]; 
                              else
                                  Isynprime[postcellE[k]*Nsyn+isyn]+=Jie*temp2[isyn];}
                      }
                  }
                  
              }
              else{
                  Pr=Pr+Ke; }
        }
          
       else{ /* If cell is inhibitory */
            
             if(refstate[j]<=0){
                Isyntot=0;
                for(isyn=0;isyn<Nsyn;isyn++){
                Isyntot+=Psyn[isyn]*Isyn[j*Nsyn+isyn];
                }
             v[j]+=fmax((Isyntot-gl[1]*(v[j]-Vleak[1])+gl[1]*DeltaT[1]*EXP((v[j]-VT[1])/DeltaT[1]))*dt/C[1],Vlb[1]-v[j]);}
             else{                 
                if(refstate[j]>1)
                   v[j]=Vth[1];
                else
                   v[j]=Vre[1];
                refstate[j]--;
             }
             
              /* If a spike occurs */
              if(v[j]>=Vth[1] && refstate[j]<=0 && ns<maxns){                                                            
                  refstate[j]=Ntref[1];
                  v[j]=Vth[1];       /* reset membrane potential */
                  s[0+2*ns]=i*dt; /* spike time */
                  s[1+2*ns]=j+1;     /* neuron index 1 */
                  ns++;           /* update total number of spikes */
                  /* For each postsynaptic target, update synaptic inputs */
                  Pr=&Wrr[(j-Ne)*Ki+Ne*Ke];
                  /*
                  for(k=0;k<Ki;k++){
                      postcell=(int)((*Pr)-1);
                      if(postcell<0 || postcell>=N){
                          mexPrintf("\n Wrf j=%d, postcell=%d\n",j, postcell);
                          mexErrMsgTxt("postcell out of bounds");
                      }
                      if (postcell<Ne){
                          for(isyn=0;isyn<Nsyn;isyn++){
                              if(syntype[isyn]==2)
                                  Isynprime[postcell*Nsyn+isyn]+=Jei*temp2[isyn];}}
                      else{
                              for(isyn=0;isyn<Nsyn;isyn++){
                                  if(syntype[isyn]==2)
                                      Isynprime[postcell*Nsyn+isyn]+=Jii*temp2[isyn];}
                      }
                      Pr++;
                  }
                  */
                  
                  for(k=0;k<(Ki);k++){
                       postcellI[k]=(int)(*Pr)-1;
                       if(postcellI[k]<0 || postcellI[k]>=N){
                         mexPrintf("\n inh j=%d, postcell=%d\n",j, postcellI[k]);
                         mexErrMsgTxt("postcell out of bounds");
                      }
                      Pr++; 
                   }
                  for(isyn=0;isyn<Nsyn;isyn++){
                      if(syntype[isyn]==2){
                          for(k=0;k<Ki;k++){
                              if (postcellI[k]<Ne)
                                  Isynprime[postcellI[k]*Nsyn+isyn]+=Jei*temp2[isyn];
                              else
                                  Isynprime[postcellI[k]*Nsyn+isyn]+=Jii*temp2[isyn];}
                       }
                  }
                  
              }
             else{
                  Pr=Pr+Ki; }
         }
      }
              
      /* Store recorded variables */
    for(jj=0;jj<Nrecord;jj++){
      if(Irecord[jj]<1 || Irecord[jj]>N+1)
         mexErrMsgTxt("Indices in Irecord must be between 1 and N");
        for(isyn=0;isyn<Nsyn;isyn++){
            Isynrecord[isyn+jj*Nsyn+i*Nrecord*Nsyn]=Isyn[(int)round(Irecord[jj]-1)*Nsyn+isyn];
        }
      vr[jj+Nrecord*i]=v[(int)round(Irecord[jj]-1)];
    }
}

/* Issue a warning if max number of spikes reached */
if(ns>=maxns)
   mexWarnMsgTxt("Maximum number of spikes reached, simulation terminated.");
   
/* Free allocated memory */
mxFree(v);
mxFree(refstate);
mxFree(Isyn);
mxFree(Isynprime);
mxFree(temp1);
mxFree(temp2);
mxFree(postcellX);
mxFree(postcellE);
mxFree(postcellI);
mxFree(syntype); 
}







