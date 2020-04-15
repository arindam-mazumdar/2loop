#include <math.h>
#include <stdio.h>
#include "commands.h" //header file containg function "write2file", "readfile2array" which can read/write 2 coloumns of float in "filename".dat

int main()
{
///////////////////////////////////////// Read the values of k,p,q,mp,mq,np from files //////////////////////
int numrow;
int krow=52;
int sampling= 10000;
float kvalue[krow], qvalue[sampling],pvalue[sampling],mpvalue[sampling],mqvalue[sampling],npvalue[sampling];
FILE *fptr;
fptr = fopen("./input/k.dat","r"); /// values of k
for(numrow=1;numrow<=krow;numrow++)
{fscanf(fptr, "%f", &kvalue[numrow]);
}
fclose(fptr);
//printf("data at k values row %f %f\n", kvalue[12],kvalue[2]);

fptr = fopen("./input/q.dat","r"); ////////// Values of q 
for(numrow=1;numrow<=sampling;numrow++)
{fscanf(fptr, "%f", &qvalue[numrow]);
}
fclose(fptr);
//printf("data at q values row %f %f\n", qvalue[12],qvalue[992]);
fptr = fopen("./input/p.dat","r");
for(numrow=1;numrow<=sampling;numrow++)
{fscanf(fptr, "%f", &pvalue[numrow]);
}
fclose(fptr);
//printf("data at p values row %f %f\n", pvalue[12],pvalue[9345]);
fptr = fopen("./input/mp.dat","r");
for(numrow=1;numrow<=sampling;numrow++)
{fscanf(fptr, "%f", &mpvalue[numrow]);
}
fclose(fptr);
fptr = fopen("./input/mq.dat","r");
for(numrow=1;numrow<=sampling;numrow++)
{fscanf(fptr, "%f", &mqvalue[numrow]);
}
fclose(fptr);
fptr = fopen("./input/np.dat","r");
for(numrow=1;numrow<=sampling;numrow++)
{fscanf(fptr, "%f", &npvalue[numrow]);
}
fclose(fptr);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////       START   THE  MONTE CARLO INTEGRATION     /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* It will be done in 5 steps.
1.) Provide k
2.) Provide p,q.... values from the sampling
3.) Calculate the value of the kernels
4.) Provide the values of P(k-q) and P(k-p-q). 
5.) Do the multiplication and write the values of the integration in a file 
*/
for(numrow=1;numrow<=krow;numrow++)
{
float k = kvalue[numrow];//0.2;
int j;
for(j=1;j<=sampling;j++)
{
float q =qvalue[j];//0.1;
float p =pvalue[j];//0.12;
float mp=mpvalue[j];//0.25;
float mq=mqvalue[j];//0.60;
float np=npvalue[j];//0.41;
/////// others can be set to zero without loss of generality ////////
float mk=0.;
float nq=0.;
float nk=0.;
//////////// A small factor of epsilon is required to avoid divergence ///////
float eps = 0.000002;

/////////////////////////////////////////////////// value of F2sym /////////////////////
float F2symv= F2sym(k,q,mq); 
/////////////////////////////////////////////////// value of F3sym-1st type /////////////////////
float F3sym1v= F3sym1(k,q,mq);
float F3sym1vp= F3sym1(k,p,mp);
/////////////////////////////////////////////////// value of F3sym-2nd type /////////////////////
float F3sym2v= F3sym2(k,p,q,mp,mq,np);
float F3sym2vneg= F3sym2(-k,-p,-q,mp,mq,np);
//////////////////////////////////////////////// value of F4Sym kernel ///////////////////////////// F4(-q, q-k, p,-p)
float k1= -q+eps;
float k2= sqrt((1 - pow(mq,2))*pow(q,2) + pow(-k + mq*q,2));
float alpha= sqrt(1-pow(mq,2))/(mq-(k/q));
float k3= p;
float k4= -p+eps;
float m1= mq;
float m2= sqrt(1/(1+pow(alpha,2)));
float m3= mp;
float m4= mp;
float n1= mq;
float n2= 0.;
float n3= np;
float n4= np;

float F4symv = (1./24.)*(F4(k1,k2,k3,k4,m1,m2,m3,m4,n1,n2,n3,n4)+F4(k1,k2,k4,k3,m1,m2,m4,m3,n1,n2,n4,n3)+F4(k1,k3,k2,k4,m1,m3,m2,m4,n1,n3,n2,n4)+F4(k1,k3,k4,k2,m1,m3,m4,m2,n1,n3,n4,n2)+F4(k1,k4,k2,k3,m1,m4,m2,m3,n1,n4,n2,n3)+F4(k1,k4,k3,k2,m1,m4,m3,m2,n1,n4,n3,n2)+F4(k2,k1,k3,k4,m2,m1,m3,m4,n2,n1,n3,n4)+
F4(k2,k1,k4,k3,m2,m1,m4,m3,n2,n1,n4,n3)+F4(k2,k3,k1,k4,m2,m3,m1,m4,n2,n3,n1,n4)+F4(k2,k3,k4,k1,m2,m3,m4,m1,n2,n3,n4,n1)+F4(k2,k4,k1,k3,m2,m4,m1,m3,n2,n4,n1,n3)+F4(k2,k4,k3,k1,m2,m4,m3,m1,n2,n4,n3,n1)+F4(k3,k1,k2,k4,m3,m1,m2,m4,n3,n1,n2,n4)+F4(k3,k1,k4,k2,m3,m1,m4,m2,n3,n1,n4,n2)+
F4(k3,k2,k1,k4,m3,m2,m1,m4,n3,n2,n1,n4)+F4(k3,k2,k4,k1,m3,m2,m4,m1,n3,n2,n4,n1)+F4(k3,k4,k1,k2,m3,m4,m1,m2,n3,n4,n1,n2)+F4(k3,k4,k2,k1,m3,m4,m2,m1,n3,n4,n2,n1)+F4(k4,k1,k2,k3,m4,m1,m2,m3,n4,n1,n2,n3)+F4(k4,k1,k3,k2,m4,m1,m3,m2,n4,n1,n3,n2)+F4(k4,k2,k1,k3,m4,m2,m1,m3,n4,n2,n1,n3)+
F4(k4,k2,k3,k1,m4,m2,m3,m1,n4,n2,n3,n1)+F4(k4,k3,k1,k2,m4,m3,m1,m2,n4,n3,n1,n2)+F4(k4,k3,k2,k1,m4,m3,m2,m1,n4,n3,n2,n1));

//////////////////////////////////////////////// Defenition of F5Sym kernel /////////////////////////////

float F5symv = (F5(k,-p+eps,p,-q+eps,q,mk,mp,mp,-mq,mq,nk,np,np,-nq,nq) + F5(k,-p+eps,p,q,-q+eps,mk,mp,mp,mq,-mq,nk,np,np,nq,-nq) + F5(k,-p+eps,-q+eps,p,q,mk,mp,-mq,mp,mq,nk,np,-nq,np,nq) + F5(k,-p+eps,-q+eps,q,p,mk,mp,-mq,mq,mp,nk,np,-nq,nq,np) + F5(k,-p+eps,q,p,-q+eps,mk,mp,mq,mp,-mq,nk,np,nq,np,-nq) + 
     F5(k,-p+eps,q,-q+eps,p,mk,mp,mq,-mq,mp,nk,np,nq,-nq,np) + F5(k,p,-p+eps,-q+eps,q,mk,mp,mp,-mq,mq,nk,np,np,-nq,nq) + F5(k,p,-p+eps,q,-q+eps,mk,mp,mp,mq,-mq,nk,np,np,nq,-nq) + F5(k,p,-q+eps,-p+eps,q,mk,mp,-mq,mp,mq,nk,np,-nq,np,nq) + 
     F5(k,p,-q+eps,q,-p+eps,mk,mp,-mq,mq,mp,nk,np,-nq,nq,np) + F5(k,p,q,-p+eps,-q+eps,mk,mp,mq,mp,-mq,nk,np,nq,np,-nq) + F5(k,p,q,-q+eps,-p+eps,mk,mp,mq,-mq,mp,nk,np,nq,-nq,np) + F5(k,-q+eps,-p+eps,p,q,mk,-mq,mp,mp,mq,nk,-nq,np,np,nq) + 
     F5(k,-q+eps,-p+eps,q,p,mk,-mq,mp,mq,mp,nk,-nq,np,nq,np) + F5(k,-q+eps,p,-p+eps,q,mk,-mq,mp,mp,mq,nk,-nq,np,np,nq) + F5(k,-q+eps,p,q,-p+eps,mk,-mq,mp,mq,mp,nk,-nq,np,nq,np) + F5(k,-q+eps,q,-p+eps,p,mk,-mq,mq,mp,mp,nk,-nq,nq,np,np) + 
     F5(k,-q+eps,q,p,-p+eps,mk,-mq,mq,mp,mp,nk,-nq,nq,np,np) + F5(k,q,-p+eps,p,-q+eps,mk,mq,mp,mp,-mq,nk,nq,np,np,-nq) + F5(k,q,-p+eps,-q+eps,p,mk,mq,mp,-mq,mp,nk,nq,np,-nq,np) + F5(k,q,p,-p+eps,-q+eps,mk,mq,mp,mp,-mq,nk,nq,np,np,-nq) + 
     F5(k,q,p,-q+eps,-p+eps,mk,mq,mp,-mq,mp,nk,nq,np,-nq,np) + F5(k,q,-q+eps,-p+eps,p,mk,mq,-mq,mp,mp,nk,nq,-nq,np,np) + F5(k,q,-q+eps,p,-p+eps,mk,mq,-mq,mp,mp,nk,nq,-nq,np,np) + F5(-p+eps,k,p,-q+eps,q,mp,mk,mp,-mq,mq,np,nk,np,-nq,nq) + 
     F5(-p+eps,k,p,q,-q+eps,mp,mk,mp,mq,-mq,np,nk,np,nq,-nq) + F5(-p+eps,k,-q+eps,p,q,mp,mk,-mq,mp,mq,np,nk,-nq,np,nq) + F5(-p+eps,k,-q+eps,q,p,mp,mk,-mq,mq,mp,np,nk,-nq,nq,np) + F5(-p+eps,k,q,p,-q+eps,mp,mk,mq,mp,-mq,np,nk,nq,np,-nq) + 
     F5(-p+eps,k,q,-q+eps,p,mp,mk,mq,-mq,mp,np,nk,nq,-nq,np) + F5(-p+eps,p,k,-q+eps,q,mp,mp,mk,-mq,mq,np,np,nk,-nq,nq) + F5(-p+eps,p,k,q,-q+eps,mp,mp,mk,mq,-mq,np,np,nk,nq,-nq) + F5(-p+eps,p,-q+eps,k,q,mp,mp,-mq,mk,mq,np,np,-nq,nk,nq) + 
     F5(-p+eps,p,-q+eps,q,k,mp,mp,-mq,mq,mk,np,np,-nq,nq,nk) + F5(-p+eps,p,q,k,-q+eps,mp,mp,mq,mk,-mq,np,np,nq,nk,-nq) + F5(-p+eps,p,q,-q+eps,k,mp,mp,mq,-mq,mk,np,np,nq,-nq,nk) + F5(-p+eps,-q+eps,k,p,q,mp,-mq,mk,mp,mq,np,-nq,nk,np,nq) + 
     F5(-p+eps,-q+eps,k,q,p,mp,-mq,mk,mq,mp,np,-nq,nk,nq,np) + F5(-p+eps,-q+eps,p,k,q,mp,-mq,mp,mk,mq,np,-nq,np,nk,nq) + F5(-p+eps,-q+eps,p,q,k,mp,-mq,mp,mq,mk,np,-nq,np,nq,nk) + F5(-p+eps,-q+eps,q,k,p,mp,-mq,mq,mk,mp,np,-nq,nq,nk,np) + 
     F5(-p+eps,-q+eps,q,p,k,mp,-mq,mq,mp,mk,np,-nq,nq,np,nk) + F5(-p+eps,q,k,p,-q+eps,mp,mq,mk,mp,-mq,np,nq,nk,np,-nq) + F5(-p+eps,q,k,-q+eps,p,mp,mq,mk,-mq,mp,np,nq,nk,-nq,np) + F5(-p+eps,q,p,k,-q+eps,mp,mq,mp,mk,-mq,np,nq,np,nk,-nq) + 
     F5(-p+eps,q,p,-q+eps,k,mp,mq,mp,-mq,mk,np,nq,np,-nq,nk) + F5(-p+eps,q,-q+eps,k,p,mp,mq,-mq,mk,mp,np,nq,-nq,nk,np) + F5(-p+eps,q,-q+eps,p,k,mp,mq,-mq,mp,mk,np,nq,-nq,np,nk) + F5(p,k,-p+eps,-q+eps,q,mp,mk,mp,-mq,mq,np,nk,np,-nq,nq) + 
     F5(p,k,-p+eps,q,-q+eps,mp,mk,mp,mq,-mq,np,nk,np,nq,-nq) + F5(p,k,-q+eps,-p+eps,q,mp,mk,-mq,mp,mq,np,nk,-nq,np,nq) + F5(p,k,-q+eps,q,-p+eps,mp,mk,-mq,mq,mp,np,nk,-nq,nq,np) + F5(p,k,q,-p+eps,-q+eps,mp,mk,mq,mp,-mq,np,nk,nq,np,-nq) + 
     F5(p,k,q,-q+eps,-p+eps,mp,mk,mq,-mq,mp,np,nk,nq,-nq,np) + F5(p,-p+eps,k,-q+eps,q,mp,mp,mk,-mq,mq,np,np,nk,-nq,nq) + F5(p,-p+eps,k,q,-q+eps,mp,mp,mk,mq,-mq,np,np,nk,nq,-nq) + F5(p,-p+eps,-q+eps,k,q,mp,mp,-mq,mk,mq,np,np,-nq,nk,nq) + 
     F5(p,-p+eps,-q+eps,q,k,mp,mp,-mq,mq,mk,np,np,-nq,nq,nk) + F5(p,-p+eps,q,k,-q+eps,mp,mp,mq,mk,-mq,np,np,nq,nk,-nq) + F5(p,-p+eps,q,-q+eps,k,mp,mp,mq,-mq,mk,np,np,nq,-nq,nk) + F5(p,-q+eps,k,-p+eps,q,mp,-mq,mk,mp,mq,np,-nq,nk,np,nq) + 
     F5(p,-q+eps,k,q,-p+eps,mp,-mq,mk,mq,mp,np,-nq,nk,nq,np) + F5(p,-q+eps,-p+eps,k,q,mp,-mq,mp,mk,mq,np,-nq,np,nk,nq) + F5(p,-q+eps,-p+eps,q,k,mp,-mq,mp,mq,mk,np,-nq,np,nq,nk) + F5(p,-q+eps,q,k,-p+eps,mp,-mq,mq,mk,mp,np,-nq,nq,nk,np) + 
     F5(p,-q+eps,q,-p+eps,k,mp,-mq,mq,mp,mk,np,-nq,nq,np,nk) + F5(p,q,k,-p+eps,-q+eps,mp,mq,mk,mp,-mq,np,nq,nk,np,-nq) + F5(p,q,k,-q+eps,-p+eps,mp,mq,mk,-mq,mp,np,nq,nk,-nq,np) + F5(p,q,-p+eps,k,-q+eps,mp,mq,mp,mk,-mq,np,nq,np,nk,-nq) + 
     F5(p,q,-p+eps,-q+eps,k,mp,mq,mp,-mq,mk,np,nq,np,-nq,nk) + F5(p,q,-q+eps,k,-p+eps,mp,mq,-mq,mk,mp,np,nq,-nq,nk,np) + F5(p,q,-q+eps,-p+eps,k,mp,mq,-mq,mp,mk,np,nq,-nq,np,nk) + F5(-q+eps,k,-p+eps,p,q,-mq,mk,mp,mp,mq,-nq,nk,np,np,nq) + 
     F5(-q+eps,k,-p+eps,q,p,-mq,mk,mp,mq,mp,-nq,nk,np,nq,np) + F5(-q+eps,k,p,-p+eps,q,-mq,mk,mp,mp,mq,-nq,nk,np,np,nq) + F5(-q+eps,k,p,q,-p+eps,-mq,mk,mp,mq,mp,-nq,nk,np,nq,np) + F5(-q+eps,k,q,-p+eps,p,-mq,mk,mq,mp,mp,-nq,nk,nq,np,np) + 
     F5(-q+eps,k,q,p,-p+eps,-mq,mk,mq,mp,mp,-nq,nk,nq,np,np) + F5(-q+eps,-p+eps,k,p,q,-mq,mp,mk,mp,mq,-nq,np,nk,np,nq) + F5(-q+eps,-p+eps,k,q,p,-mq,mp,mk,mq,mp,-nq,np,nk,nq,np) + F5(-q+eps,-p+eps,p,k,q,-mq,mp,mp,mk,mq,-nq,np,np,nk,nq) + 
     F5(-q+eps,-p+eps,p,q,k,-mq,mp,mp,mq,mk,-nq,np,np,nq,nk) + F5(-q+eps,-p+eps,q,k,p,-mq,mp,mq,mk,mp,-nq,np,nq,nk,np) + F5(-q+eps,-p+eps,q,p,k,-mq,mp,mq,mp,mk,-nq,np,nq,np,nk) + F5(-q+eps,p,k,-p+eps,q,-mq,mp,mk,mp,mq,-nq,np,nk,np,nq) + 
     F5(-q+eps,p,k,q,-p+eps,-mq,mp,mk,mq,mp,-nq,np,nk,nq,np) + F5(-q+eps,p,-p+eps,k,q,-mq,mp,mp,mk,mq,-nq,np,np,nk,nq) + F5(-q+eps,p,-p+eps,q,k,-mq,mp,mp,mq,mk,-nq,np,np,nq,nk) + F5(-q+eps,p,q,k,-p+eps,-mq,mp,mq,mk,mp,-nq,np,nq,nk,np) + 
     F5(-q+eps,p,q,-p+eps,k,-mq,mp,mq,mp,mk,-nq,np,nq,np,nk) + F5(-q+eps,q,k,-p+eps,p,-mq,mq,mk,mp,mp,-nq,nq,nk,np,np) + F5(-q+eps,q,k,p,-p+eps,-mq,mq,mk,mp,mp,-nq,nq,nk,np,np) + F5(-q+eps,q,-p+eps,k,p,-mq,mq,mp,mk,mp,-nq,nq,np,nk,np) + 
     F5(-q+eps,q,-p+eps,p,k,-mq,mq,mp,mp,mk,-nq,nq,np,np,nk) + F5(-q+eps,q,p,k,-p+eps,-mq,mq,mp,mk,mp,-nq,nq,np,nk,np) + F5(-q+eps,q,p,-p+eps,k,-mq,mq,mp,mp,mk,-nq,nq,np,np,nk) + F5(q,k,-p+eps,p,-q+eps,mq,mk,mp,mp,-mq,nq,nk,np,np,-nq) + 
     F5(q,k,-p+eps,-q+eps,p,mq,mk,mp,-mq,mp,nq,nk,np,-nq,np) + F5(q,k,p,-p+eps,-q+eps,mq,mk,mp,mp,-mq,nq,nk,np,np,-nq) + F5(q,k,p,-q+eps,-p+eps,mq,mk,mp,-mq,mp,nq,nk,np,-nq,np) + F5(q,k,-q+eps,-p+eps,p,mq,mk,-mq,mp,mp,nq,nk,-nq,np,np) + 
     F5(q,k,-q+eps,p,-p+eps,mq,mk,-mq,mp,mp,nq,nk,-nq,np,np) + F5(q,-p+eps,k,p,-q+eps,mq,mp,mk,mp,-mq,nq,np,nk,np,-nq) + F5(q,-p+eps,k,-q+eps,p,mq,mp,mk,-mq,np,nq,np,nk,-nq,np) + F5(q,-p+eps,p,k,-q+eps,mq,mp,mp,mk,-mq,nq,np,np,nk,-nq) + 
     F5(q,-p+eps,p,-q+eps,k,mq,mp,mp,-mq,mk,nq,np,np,-nq,nk) + F5(q,-p+eps,-q+eps,k,p,mq,mp,-mq,mk,mp,nq,np,-nq,nk,np) + F5(q,-p+eps,-q+eps,p,k,mq,mp,-mq,mp,mk,nq,np,-nq,np,nk) + F5(q,p,k,-p+eps,-q+eps,mq,mp,mk,mp,-mq,nq,np,nk,np,-nq) + 
     F5(q,p,k,-q+eps,-p+eps,mq,mp,mk,-mq,mp,nq,np,nk,-nq,np) + F5(q,p,-p+eps,k,-q+eps,mq,mp,mp,mk,-mq,nq,np,np,nk,-nq) + F5(q,p,-p+eps,-q+eps,k,mq,mp,mp,-mq,mk,nq,np,np,-nq,nk) + F5(q,p,-q+eps,k,-p+eps,mq,mp,-mq,mk,mp,nq,np,-nq,nk,np) + 
     F5(q,p,-q+eps,-p+eps,k,mq,mp,-mq,mp,mk,nq,np,-nq,np,nk) + F5(q,p,-q+eps,-p+eps,k,mq,mp,-mq,mp,mk,nq,np,-nq,np,nk) + F5(q,-q+eps,k,p,-p+eps,mq,-mq,mk,mp,mp,nq,-nq,nk,np,np) + F5(q,-q+eps,-p+eps,k,p,mq,-mq,mp,mk,mp,nq,-nq,np,nk,np) + 
     F5(q,-q+eps,-p+eps,p,k,mq,-mq,mp,mp,mk,nq,-nq,np,np,nk) + F5(q,-q+eps,p,k,-p+eps,mq,-mq,mp,mk,mp,nq,-nq,np,nk,np) + F5(q,-q+eps,p,-p+eps,k,mq,-mq,mp,mp,mk,nq,-nq,np,np,nk))/120. ;

///////////////////////////////////////////////////////// Write the values of the kernels ///////////
char file[32];
sprintf(file, "./input/%d/F2.dat", numrow);
fptr = fopen(file,"a");
fprintf(fptr,"%f\n", F2symv);
fclose(fptr);
sprintf(file, "./input/%d/F3Iq.dat", numrow);
fptr = fopen(file,"a");
fprintf(fptr,"%f\n", F3sym1v);
fclose(fptr);
sprintf(file, "./input/%d/F3Ip.dat", numrow);
fptr = fopen(file,"a");
fprintf(fptr,"%f\n", F3sym1vp);
fclose(fptr);
sprintf(file, "./input/%d/F3II.dat", numrow);
fptr = fopen(file,"a");
fprintf(fptr,"%f\n", F3sym2v);
fclose(fptr);
sprintf(file, "./input/%d/F3IIneg.dat", numrow);
fptr = fopen(file,"a");
fprintf(fptr,"%f\n", F3sym2vneg);
fclose(fptr);
sprintf(file, "./input/%d/F4.dat", numrow);
fptr = fopen(file,"a");
fprintf(fptr,"%f\n", F4symv);
fclose(fptr);
sprintf(file, "./input/%d/F5.dat", numrow);
fptr = fopen(file,"a");
fprintf(fptr,"%f\n", F5symv);
fclose(fptr);


//write2file(file, kvalue[numrow], F5symv);
//printf("Value of F2 is %f\n",b);
//printf("Values of F5 F4 F3I F3II F2 are \n %.6f \n %.6f \n %.6f \n %.6f \n %.6f \n", F5symv,F4symv,F3sym1v,F3sym2v,F2symv); //

}
printf("Finished the writting the kernels for k= %f\n",kvalue[numrow]);

}
return 0;
} 

/*float P_k_q[sampling], P_k_q_p[sampling]
char file[32];
int x=44;
sprintf(file, "./input/%d/P11_k-q.dat", x);
fptr = fopen(file,"r");
for(numrow=1;numrow<=sampling;numrow++)
{fscanf(fptr, "%f", &kvalue[numrow]);
}
fclose(fptr); */

