#include <math.h>
//#include <stdio.h>

// This file contains

//One loop kernels F2 and F3
float F2sym(float k, float q, float mq)
{
float temp;
temp = (pow(k,2)*(7*k*mq + (3 - 10*pow(mq,2))*q))/(14.*q*(pow(k,2) - 2*k*mq*q + pow(q,2)));
return (temp);
}
///////////////////////////////////////////////// F3 functions ///////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
float F3sym1(float k, float q, float mq) // F3sym (k,q,-q)
{float temp;
temp = -(k*(21*pow(k,5)*pow(mq,2) + pow(k,3)*(-10 + 44*pow(mq,2) - 76*pow(mq,4))*pow(q,2) - k*(10 - 59*pow(mq,2) + 28*pow(mq,4))*pow(q,4)))/
   (126.*pow(q,2)*(pow(k,2) - 2*k*mq*q + pow(q,2))*(pow(k,2) + 2*k*mq*q + pow(q,2)));
return(temp);
}
////////////////////////////////////////////////////////////////////////////
float F3sym2(float k, float p, float q, float mp, float mq, float np)// F3sym(q,p,k-q-p)
{ float temp1;
temp1 = (pow(k,2)*(21*pow(k,6)*mp*mq*(pow(p,2) + 2*(mp*mq + sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np)*p*q + pow(q,2)) + 
       pow(p,3)*pow(q,3)*((8 + 56*pow(mp,3)*sqrt(1 - pow(mp,2))*mq*sqrt(1 - pow(mq,2))*np + 
             2*mp*sqrt(1 - pow(mp,2))*mq*sqrt(1 - pow(mq,2))*(-15 + 14*pow(mq,2))*np - 8*pow(np,2) - 14*pow(mq,4)*pow(np,2) + 
             pow(mq,2)*(-7 + 22*pow(np,2)) + 28*pow(mp,4)*(-pow(np,2) + pow(mq,2)*(1 + pow(np,2))) + 
             pow(mp,2)*(-21 + 36*pow(np,2) + 14*pow(mq,4)*(1 + pow(np,2)) - 2*pow(mq,2)*(11 + 25*pow(np,2))))*pow(p,2) + 
          2*(-(mp*mq*(6 + 24*pow(np,2) + 42*pow(mq,4)*pow(np,2) + pow(mq,2)*(7 - 66*pow(np,2)))) + 
             14*pow(mp,4)*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*(-pow(np,2) + pow(mq,2)*(3 + pow(np,2))) + 
             14*pow(mp,5)*(-3*mq*pow(np,2) + pow(mq,3)*(1 + 3*pow(np,2))) + 
             pow(mp,2)*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*
              (-7 + 22*pow(np,2) + 14*pow(mq,4)*(3 + pow(np,2)) - 12*pow(mq,2)*(2 + 3*pow(np,2))) + 
             sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*
              (8 - 8*pow(np,2) - 14*pow(mq,4)*pow(np,2) + pow(mq,2)*(-7 + 22*pow(np,2))) + 
             pow(mp,3)*mq*(-7 + 66*pow(np,2) + 14*pow(mq,4)*(1 + 3*pow(np,2)) - 4*pow(mq,2)*(2 + 27*pow(np,2))))*p*q + 
          (8 + 28*pow(mp,3)*sqrt(1 - pow(mp,2))*mq*sqrt(1 - pow(mq,2))*np + 
             2*mp*sqrt(1 - pow(mp,2))*mq*sqrt(1 - pow(mq,2))*(-15 + 28*pow(mq,2))*np - 8*pow(np,2) - 28*pow(mq,4)*pow(np,2) + 
             3*pow(mq,2)*(-7 + 12*pow(np,2)) + 14*pow(mp,4)*(-pow(np,2) + pow(mq,2)*(1 + pow(np,2))) + 
             pow(mp,2)*(-7 + 22*pow(np,2) + 28*pow(mq,4)*(1 + pow(np,2)) - 2*pow(mq,2)*(11 + 25*pow(np,2))))*pow(q,2)) + 
       pow(k,2)*p*q*((8 - 2*mp*sqrt(1 - pow(mp,2))*mq*sqrt(1 - pow(mq,2))*np + 
             56*pow(mp,3)*sqrt(1 - pow(mp,2))*mq*sqrt(1 - pow(mq,2))*np - 3*pow(np,2) + pow(mq,2)*(-25 + 3*pow(np,2)) + 
             pow(mp,2)*(-26 + 31*pow(np,2) + pow(mq,2)*(57 - 31*pow(np,2))) + 
             28*pow(mp,4)*(-pow(np,2) + pow(mq,2)*(1 + pow(np,2))))*pow(p,4) + 
          (pow(mp,3)*mq*(-86 + 60*pow(np,2) + pow(mq,2)*(260 - 60*pow(np,2))) + 
             2*pow(mp,2)*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*(-22 + 17*pow(np,2) + pow(mq,2)*(173 - 17*pow(np,2))) + 
             28*pow(mp,4)*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*(-pow(np,2) + pow(mq,2)*(7 + pow(np,2))) + 
             2*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*(8 - 3*pow(np,2) + pow(mq,2)*(-35 + 3*pow(np,2))) + 
             28*pow(mp,5)*(-5*mq*pow(np,2) + pow(mq,3)*(3 + 5*pow(np,2))) - 
             mp*mq*(13 - 80*pow(np,2) + 8*pow(mq,2)*(7 + 10*pow(np,2))))*pow(p,3)*q + 
          (16 - 6*pow(np,2) - 6*pow(mq,4)*pow(np,2) + 3*pow(mq,2)*(-23 + 4*pow(np,2)) + 
             2*mp*sqrt(1 - pow(mp,2))*mq*sqrt(1 - pow(mq,2))*np*(15 + 44*pow(np,2) + pow(mq,2)*(20 - 44*pow(np,2))) + 
             8*pow(mp,3)*sqrt(1 - pow(mp,2))*mq*sqrt(1 - pow(mq,2))*np*(5 - 11*pow(np,2) + pow(mq,2)*(89 + 11*pow(np,2))) + 
             pow(mp,4)*(-6*pow(np,2) + pow(mq,2)*(34 - 482*pow(np,2)) + 8*pow(mq,4)*(39 + 61*pow(np,2))) + 
             pow(mp,2)*(-69 + 12*pow(np,2) + pow(mq,4)*(34 - 482*pow(np,2)) + pow(mq,2)*(36 + 470*pow(np,2))))*pow(p,2)*
           pow(q,2) + (-2*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*
              (-8 + 3*pow(np,2) + 14*pow(mq,4)*pow(np,2) + pow(mq,2)*(22 - 17*pow(np,2))) + 
             2*pow(mp,2)*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*
              (-35 + 3*pow(np,2) + pow(mq,2)*(173 - 17*pow(np,2)) + 14*pow(mq,4)*(7 + pow(np,2))) + 
             4*pow(mp,3)*mq*(pow(mq,2)*(65 - 15*pow(np,2)) + 7*pow(mq,4)*(3 + 5*pow(np,2)) - 2*(7 + 10*pow(np,2))) + 
             mp*mq*(-13 + 80*pow(np,2) - 140*pow(mq,4)*pow(np,2) + pow(mq,2)*(-86 + 60*pow(np,2))))*p*pow(q,3) + 
          (8 + 2*mp*sqrt(1 - pow(mp,2))*mq*sqrt(1 - pow(mq,2))*(-1 + 28*pow(mq,2))*np - 3*pow(np,2) - 28*pow(mq,4)*pow(np,2) + 
             pow(mq,2)*(-26 + 31*pow(np,2)) + pow(mp,2)*
              (-25 + 3*pow(np,2) + pow(mq,2)*(57 - 31*pow(np,2)) + 28*pow(mq,4)*(1 + pow(np,2))))*pow(q,4)) - 
       9*pow(k,5)*(pow(mq,3)*(-pow(np,2) + pow(mp,2)*(17 + pow(np,2)))*p*pow(q,2) - 
          mp*q*((2 + (-1 + pow(mp,2))*pow(np,2))*pow(p,2) + 2*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*p*q + pow(q,2)) + 
          mp*pow(mq,2)*q*((6 - pow(np,2) + pow(mp,2)*(17 + pow(np,2)))*pow(p,2) + 
             18*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*p*q + 8*pow(q,2)) + 
          mq*p*((-1 + 8*pow(mp,2))*pow(p,2) + 2*sqrt(1 - pow(mp,2))*(-1 + 9*pow(mp,2))*sqrt(1 - pow(mq,2))*np*p*q + 
             (-2 + pow(np,2) - pow(mp,2)*(-6 + pow(np,2)))*pow(q,2))) + 
       pow(k,4)*(8*pow(mp,4)*(-4*pow(np,2) + pow(mq,2)*(19 + 4*pow(np,2)))*pow(p,3)*q + 
          4*pow(mp,3)*mq*pow(p,2)*(15*pow(p,2) + 46*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*p*q + 
             6*(-5*pow(np,2) + pow(mq,2)*(21 + 5*pow(np,2)))*pow(q,2)) - 
          2*pow(mp,2)*p*q*((31 - 15*pow(np,2) + pow(mq,2)*(-113 + 15*pow(np,2)))*pow(p,2) + 
             2*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*(15 + pow(np,2) - pow(mq,2)*(155 + pow(np,2)))*p*q + 
             (15 + pow(np,2) - 4*pow(mq,4)*(19 + 4*pow(np,2)) + pow(mq,2)*(-113 + 15*pow(np,2)))*pow(q,2)) + 
          2*p*q*((4 + pow(np,2) - pow(mq,2)*(15 + pow(np,2)))*pow(p,2) + 
             2*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*(4 + pow(np,2) - pow(mq,2)*(15 + pow(np,2)))*p*q + 
             (4 + pow(np,2) - 16*pow(mq,4)*pow(np,2) + pow(mq,2)*(-31 + 15*pow(np,2)))*pow(q,2)) + 
          mp*mq*(3*pow(p,4) + 32*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*pow(p,3)*q - 
             6*(7 + 20*(-1 + pow(mq,2))*pow(np,2))*pow(p,2)*pow(q,2) + 
             8*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*(4 + 23*pow(mq,2))*np*p*pow(q,3) + 3*(1 + 20*pow(mq,2))*pow(q,4))) + 
       pow(k,3)*(-28*pow(mq,5)*(-pow(np,2) + pow(mp,2)*(1 + pow(np,2)))*p*pow(q,4) - 
          2*mp*pow(mq,4)*p*pow(q,3)*(-3*(8 + 45*pow(np,2))*p + pow(mp,2)*(221 + 135*pow(np,2))*p + 
             28*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*q) + 
          2*pow(mq,3)*p*pow(q,2)*((12 + pow(np,2) + 2*pow(mp,2)*(-44 + 67*pow(np,2)) - pow(mp,4)*(221 + 135*pow(np,2)))*
              pow(p,2) + sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*(24 + 13*pow(np,2) - pow(mp,2)*(343 + 13*pow(np,2)))*p*q + 
             (26 - 3*pow(np,2) + pow(mp,2)*(-115 + 3*pow(np,2)))*pow(q,2)) + 
          mp*q*((2 - 22*pow(np,2) + 28*pow(mp,4)*pow(np,2) + pow(mp,2)*(52 - 6*pow(np,2)))*pow(p,4) + 
             2*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*(-2 - 13*pow(np,2) + pow(mp,2)*(24 + 13*pow(np,2)))*pow(p,3)*q + 
             (11 - 2*pow(np,2) + 2*pow(mp,2)*(12 + pow(np,2)))*pow(p,2)*pow(q,2) + 
             28*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*p*pow(q,3) + 9*pow(q,4)) - 
          mp*pow(mq,2)*q*((-15 - 22*pow(np,2) + pow(mp,2)*(230 - 6*pow(np,2)) + 28*pow(mp,4)*(1 + pow(np,2)))*pow(p,4) + 
             2*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*(33 - 13*pow(np,2) + pow(mp,2)*(343 + 13*pow(np,2)))*pow(p,3)*q + 
             (-73 + 268*pow(np,2) - 4*pow(mp,2)*(-44 + 67*pow(np,2)))*pow(p,2)*pow(q,2) + 
             140*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*p*pow(q,3) + 30*pow(q,4)) + 
          mq*p*((9 - 30*pow(mp,2))*pow(p,4) - 28*sqrt(1 - pow(mp,2))*(-1 + 5*pow(mp,2) + 2*pow(mp,4))*sqrt(1 - pow(mq,2))*np*
              pow(p,3)*q + (11 - 2*pow(np,2) + pow(mp,2)*(73 - 268*pow(np,2)) + 6*pow(mp,4)*(8 + 45*pow(np,2)))*pow(p,2)*
              pow(q,2) + 2*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*(-2 - 13*pow(np,2) + pow(mp,2)*(-33 + 13*pow(np,2)))*p*
              pow(q,3) + (2 - 22*pow(np,2) + pow(mp,2)*(15 + 22*pow(np,2)))*pow(q,4))) - 
       k*pow(p,2)*pow(q,2)*(28*pow(mq,5)*(-pow(np,2) + pow(mp,2)*(1 - 5*pow(np,2)) + pow(mp,4)*(2 + 6*pow(np,2)))*p*
           pow(q,2) + 2*mp*pow(mq,4)*q*((-14 - 51*pow(np,2) + pow(mp,2)*(45 - 33*pow(np,2)) + 28*pow(mp,4)*(1 + 3*pow(np,2)))*
              pow(p,2) + 28*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*(1 - pow(np,2) + pow(mp,2)*(3 + pow(np,2)))*p*q + 
             28*(-pow(np,2) + pow(mp,2)*(1 + pow(np,2)))*pow(q,2)) + 
          pow(mq,3)*((-14 - 3*pow(np,2) + pow(mp,2)*(3 - 53*pow(np,2)) + 56*pow(mp,4)*(1 + pow(np,2)))*pow(p,3) + 
             2*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*
              (-14 - 3*pow(np,2) + pow(mp,2)*(93 - 25*pow(np,2)) + 28*pow(mp,4)*(3 + pow(np,2)))*pow(p,2)*q + 
             6*(-7 + 6*pow(np,2) + pow(mp,4)*(15 - 11*pow(np,2)) + pow(mp,2)*(1 + 5*pow(np,2)))*p*pow(q,2) + 
             112*pow(mp,2)*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*pow(q,3)) + 
          mq*((7 + 3*pow(np,2) - 56*pow(mp,4)*pow(np,2) + pow(mp,2)*(-31 + 53*pow(np,2)))*pow(p,3) + 
             2*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*
              (2 + 3*pow(np,2) - 28*pow(mp,4)*(-1 + pow(np,2)) + pow(mp,2)*(-1 + 25*pow(np,2)))*pow(p,2)*q - 
             (2 + 8*pow(np,2) - 5*pow(mp,2)*(-9 + 22*pow(np,2)) + 2*pow(mp,4)*(14 + 51*pow(np,2)))*p*pow(q,2) + 
             2*sqrt(1 - pow(mp,2))*(-5 + 3*pow(mp,2))*sqrt(1 - pow(mq,2))*np*pow(q,3)) - 
          mp*((2 + 42*pow(mp,2))*pow(p,2)*q - 6*pow(1 - pow(mp,2),1.5)*sqrt(1 - pow(mq,2))*pow(np,3)*p*pow(q,2) + 
             7*(-1 + 2*pow(mp,2))*pow(q,3) + (-1 + pow(mp,2))*pow(np,2)*q*(4*(-2 + 7*pow(mp,2))*pow(p,2) + 3*pow(q,2)) + 
             2*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*p*(5*pow(p,2) + 2*(-1 + 7*pow(mp,2))*pow(q,2))) + 
          mp*pow(mq,2)*((-45 + 6*pow(mp,2) + 28*pow(mp,4))*pow(p,2)*q + 
             50*pow(1 - pow(mp,2),1.5)*sqrt(1 - pow(mq,2))*pow(np,3)*p*pow(q,2) + (-31 + 3*pow(mp,2))*pow(q,3) - 
             (-1 + pow(mp,2))*pow(np,2)*q*(10*(11 + 14*pow(mp,2))*pow(p,2) + 53*pow(q,2)) + 
             2*sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np*p*((3 + 56*pow(mp,2))*pow(p,2) + (-1 + 93*pow(mp,2))*pow(q,2))))))/
   (126.*p*(pow(k,2) - 2*k*mp*p + pow(p,2))*q*(pow(k,2) - 2*k*mq*q + pow(q,2))*
     (pow(p,2) + 2*(mp*mq + sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np)*p*q + pow(q,2))*
     (pow(k,2) + pow(p,2) + 2*(mp*mq + sqrt(1 - pow(mp,2))*sqrt(1 - pow(mq,2))*np)*p*q + pow(q,2) - 2*k*(mp*p + mq*q)));
return (temp1);
}
////////////////////////////////////////////////////// F4 functions ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
float F4(float k1, float k2, float k3, float k4, float m1, float m2, float m3, float m4, float n1, float n2, float n3, float n4)
{
float temp3;
temp3 = (((((3*(k1*m1*(k1*m1 + k2*m2) + 
                 k1*sqrt(1 - pow(m1,2))*n1*
                  (k1*sqrt(1 - pow(m1,2))*n1 + 
                    k2*sqrt(1 - pow(m2,2))*n2) + 
                 k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*
                  (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
                    k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))))/
             (pow(k1,2)*pow(m1,2) + 
               pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + 
               pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2))) + 
            (2*(k1*k2*m1*m2 + 
                 k1*k2*sqrt(1 - pow(m1,2))*sqrt(1 - pow(m2,2))*
                  n1*n2 + k1*k2*sqrt(1 - pow(m1,2))*
                  sqrt(1 - pow(m2,2))*sqrt(1 - pow(n1,2))*
                  sqrt(1 - pow(n2,2)))*
               (pow(k1*m1 + k2*m2,2) + 
                 pow(k1*sqrt(1 - pow(m1,2))*n1 + 
                   k2*sqrt(1 - pow(m2,2))*n2,2) + 
                 pow(k1*sqrt(1 - pow(m1,2))*
                    sqrt(1 - pow(n1,2)) + 
                   k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2)))/
             ((pow(k1,2)*pow(m1,2) + 
                 pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + 
                 pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2)))*
               (pow(k2,2)*pow(m2,2) + 
                 pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + 
                 pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))))*
          ((3*((k1*m1 + k2*m2)*(k1*m1 + k2*m2 + k3*m3) + 
                 (k1*sqrt(1 - pow(m1,2))*n1 + 
                    k2*sqrt(1 - pow(m2,2))*n2)*
                  (k1*sqrt(1 - pow(m1,2))*n1 + 
                    k2*sqrt(1 - pow(m2,2))*n2 + 
                    k3*sqrt(1 - pow(m3,2))*n3) + 
                 (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
                    k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))*
                  (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
                    k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                    k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))))/
             (pow(k1*m1 + k2*m2,2) + 
               pow(k1*sqrt(1 - pow(m1,2))*n1 + 
                 k2*sqrt(1 - pow(m2,2))*n2,2) + 
               pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
                 k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2)) + 
            (3*(k3*(k1*m1 + k2*m2)*m3 + 
                 k3*sqrt(1 - pow(m3,2))*
                  (k1*sqrt(1 - pow(m1,2))*n1 + 
                    k2*sqrt(1 - pow(m2,2))*n2)*n3 + 
                 k3*sqrt(1 - pow(m3,2))*
                  (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
                    k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))*
                  sqrt(1 - pow(n3,2)))*
               (pow(k1*m1 + k2*m2 + k3*m3,2) + 
                 pow(k1*sqrt(1 - pow(m1,2))*n1 + 
                   k2*sqrt(1 - pow(m2,2))*n2 + 
                   k3*sqrt(1 - pow(m3,2))*n3,2) + 
                 pow(k1*sqrt(1 - pow(m1,2))*
                    sqrt(1 - pow(n1,2)) + 
                   k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                   k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)))/
             ((pow(k1*m1 + k2*m2,2) + 
                 pow(k1*sqrt(1 - pow(m1,2))*n1 + 
                   k2*sqrt(1 - pow(m2,2))*n2,2) + 
                 pow(k1*sqrt(1 - pow(m1,2))*
                    sqrt(1 - pow(n1,2)) + 
                   k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2))*
               (pow(k3,2)*pow(m3,2) + 
                 pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + 
                 pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))))))/126. \
+ ((3*(k1*m1*(k1*m1 + k2*m2 + k3*m3) + 
               k1*sqrt(1 - pow(m1,2))*n1*
                (k1*sqrt(1 - pow(m1,2))*n1 + 
                  k2*sqrt(1 - pow(m2,2))*n2 + 
                  k3*sqrt(1 - pow(m3,2))*n3) + 
               k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*
                (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
                  k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                  k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))))*
             ((5*(k2*m2*(k2*m2 + k3*m3) + 
                    k2*sqrt(1 - pow(m2,2))*n2*
                     (k2*sqrt(1 - pow(m2,2))*n2 + 
                       k3*sqrt(1 - pow(m3,2))*n3) + 
                    k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                     (k2*sqrt(1 - pow(m2,2))*
                       sqrt(1 - pow(n2,2)) + 
                       k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))\
)))/(pow(k2,2)*pow(m2,2) + 
                  pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + 
                  pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
               ((k2*k3*m2*m3 + 
                    k2*k3*sqrt(1 - pow(m2,2))*
                     sqrt(1 - pow(m3,2))*n2*n3 + 
                    k2*k3*sqrt(1 - pow(m2,2))*
                     sqrt(1 - pow(m3,2))*sqrt(1 - pow(n2,2))*
                     sqrt(1 - pow(n3,2)))*
                  (pow(k2*m2 + k3*m3,2) + 
                    pow(k2*sqrt(1 - pow(m2,2))*n2 + 
                      k3*sqrt(1 - pow(m3,2))*n3,2) + 
                    pow(k2*sqrt(1 - pow(m2,2))*
                       sqrt(1 - pow(n2,2)) + 
                      k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),
                     2)))/
                ((pow(k2,2)*pow(m2,2) + 
                    pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + 
                    pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                  (pow(k3,2)*pow(m3,2) + 
                    pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + 
                    pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))))))/
           (7.*(pow(k1,2)*pow(m1,2) + 
               pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + 
               pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2)))) + 
          (3*(k1*m1*(k2*m2 + k3*m3) + 
               k1*sqrt(1 - pow(m1,2))*n1*
                (k2*sqrt(1 - pow(m2,2))*n2 + 
                  k3*sqrt(1 - pow(m3,2))*n3) + 
               k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*
                (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                  k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))))*
             (pow(k1*m1 + k2*m2 + k3*m3,2) + 
               pow(k1*sqrt(1 - pow(m1,2))*n1 + 
                 k2*sqrt(1 - pow(m2,2))*n2 + 
                 k3*sqrt(1 - pow(m3,2))*n3,2) + 
               pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
                 k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                 k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))*
             ((3*(k2*m2*(k2*m2 + k3*m3) + 
                    k2*sqrt(1 - pow(m2,2))*n2*
                     (k2*sqrt(1 - pow(m2,2))*n2 + 
                       k3*sqrt(1 - pow(m3,2))*n3) + 
                    k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                     (k2*sqrt(1 - pow(m2,2))*
                       sqrt(1 - pow(n2,2)) + 
                       k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))\
))/(pow(k2,2)*pow(m2,2) + 
                  pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + 
                  pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
               (2*(k2*k3*m2*m3 + 
                    k2*k3*sqrt(1 - pow(m2,2))*
                     sqrt(1 - pow(m3,2))*n2*n3 + 
                    k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*
                     sqrt(1 - pow(n2,2))*sqrt(1 - pow(n3,2)))*
                  (pow(k2*m2 + k3*m3,2) + 
                    pow(k2*sqrt(1 - pow(m2,2))*n2 + 
                      k3*sqrt(1 - pow(m3,2))*n3,2) + 
                    pow(k2*sqrt(1 - pow(m2,2))*
                       sqrt(1 - pow(n2,2)) + 
                      k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)\
))/((pow(k2,2)*pow(m2,2) + 
                    pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + 
                    pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                  (pow(k3,2)*pow(m3,2) + 
                    pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + 
                    pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))))))/
           (7.*(pow(k1,2)*pow(m1,2) + 
               pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + 
               pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2)))*
             (pow(k2*m2 + k3*m3,2) + 
               pow(k2*sqrt(1 - pow(m2,2))*n2 + 
                 k3*sqrt(1 - pow(m3,2))*n3,2) + 
               pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                 k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))))/18.)*
     ((9*((k1*m1 + k2*m2 + k3*m3)*(k1*m1 + k2*m2 + k3*m3 + k4*m4) + 
            (k1*sqrt(1 - pow(m1,2))*n1 + 
               k2*sqrt(1 - pow(m2,2))*n2 + 
               k3*sqrt(1 - pow(m3,2))*n3)*
             (k1*sqrt(1 - pow(m1,2))*n1 + 
               k2*sqrt(1 - pow(m2,2))*n2 + 
               k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4\
) + (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
               k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
               k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
             (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
               k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
               k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
               k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
        (pow(k1*m1 + k2*m2 + k3*m3,2) + 
          pow(k1*sqrt(1 - pow(m1,2))*n1 + 
            k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) \
+ pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
            k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
            k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)) + 
       ((k4*(k1*m1 + k2*m2 + k3*m3)*m4 + 
            k4*sqrt(1 - pow(m4,2))*
             (k1*sqrt(1 - pow(m1,2))*n1 + 
               k2*sqrt(1 - pow(m2,2))*n2 + 
               k3*sqrt(1 - pow(m3,2))*n3)*n4 + 
            k4*sqrt(1 - pow(m4,2))*
             (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
               k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
               k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
             sqrt(1 - pow(n4,2)))*
          (pow(k1*m1 + k2*m2 + k3*m3 + k4*m4,2) + 
            pow(k1*sqrt(1 - pow(m1,2))*n1 + 
              k2*sqrt(1 - pow(m2,2))*n2 + 
              k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,
             2) + pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
              k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
              k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
              k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
        ((pow(k1*m1 + k2*m2 + k3*m3,2) + 
            pow(k1*sqrt(1 - pow(m1,2))*n1 + 
              k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,
             2) + pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
              k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
              k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))*
          (pow(k4,2)*pow(m4,2) + 
            pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + 
            pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/33. + 
  (((3*(k1*m1*(k1*m1 + k2*m2) + 
            k1*sqrt(1 - pow(m1,2))*n1*
             (k1*sqrt(1 - pow(m1,2))*n1 + 
               k2*sqrt(1 - pow(m2,2))*n2) + 
            k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*
             (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
               k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))))/
        (pow(k1,2)*pow(m1,2) + 
          pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + 
          pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2))) + 
       (2*(k1*k2*m1*m2 + k1*k2*sqrt(1 - pow(m1,2))*
             sqrt(1 - pow(m2,2))*n1*n2 + 
            k1*k2*sqrt(1 - pow(m1,2))*sqrt(1 - pow(m2,2))*
             sqrt(1 - pow(n1,2))*sqrt(1 - pow(n2,2)))*
          (pow(k1*m1 + k2*m2,2) + 
            pow(k1*sqrt(1 - pow(m1,2))*n1 + 
              k2*sqrt(1 - pow(m2,2))*n2,2) + 
            pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
              k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2)))/
        ((pow(k1,2)*pow(m1,2) + 
            pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + 
            pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2)))*
          (pow(k2,2)*pow(m2,2) + 
            pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + 
            pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))))*
     ((9*((k1*m1 + k2*m2)*(k1*m1 + k2*m2 + k3*m3 + k4*m4) + 
            (k1*sqrt(1 - pow(m1,2))*n1 + 
               k2*sqrt(1 - pow(m2,2))*n2)*
             (k1*sqrt(1 - pow(m1,2))*n1 + 
               k2*sqrt(1 - pow(m2,2))*n2 + 
               k3*sqrt(1 - pow(m3,2))*n3 + 
               k4*sqrt(1 - pow(m4,2))*n4) + 
            (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
               k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))*
             (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
               k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
               k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
               k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
          ((5*(k3*m3*(k3*m3 + k4*m4) + 
                 k3*sqrt(1 - pow(m3,2))*n3*
                  (k3*sqrt(1 - pow(m3,2))*n3 + 
                    k4*sqrt(1 - pow(m4,2))*n4) + 
                 k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                  (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                    k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
             (pow(k3,2)*pow(m3,2) + 
               pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + 
               pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
            ((k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*
                  sqrt(1 - pow(m4,2))*n3*n4 + 
                 k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*
                  sqrt(1 - pow(n3,2))*sqrt(1 - pow(n4,2)))*
               (pow(k3*m3 + k4*m4,2) + 
                 pow(k3*sqrt(1 - pow(m3,2))*n3 + 
                   k4*sqrt(1 - pow(m4,2))*n4,2) + 
                 pow(k3*sqrt(1 - pow(m3,2))*
                    sqrt(1 - pow(n3,2)) + 
                   k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
             ((pow(k3,2)*pow(m3,2) + 
                 pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + 
                 pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
               (pow(k4,2)*pow(m4,2) + 
                 pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + 
                 pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/
        (7.*(pow(k1*m1 + k2*m2,2) + 
            pow(k1*sqrt(1 - pow(m1,2))*n1 + 
              k2*sqrt(1 - pow(m2,2))*n2,2) + 
            pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
              k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2))) + 
       (((k1*m1 + k2*m2)*(k3*m3 + k4*m4) + 
            (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2)*
             (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) \
+ (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
               k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))*
             (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
               k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
          (pow(k1*m1 + k2*m2 + k3*m3 + k4*m4,2) + 
            pow(k1*sqrt(1 - pow(m1,2))*n1 + 
              k2*sqrt(1 - pow(m2,2))*n2 + 
              k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,
             2) + pow(k1*sqrt(1 - pow(m1,2))*
               sqrt(1 - pow(n1,2)) + 
              k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
              k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
              k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
          ((3*(k3*m3*(k3*m3 + k4*m4) + 
                 k3*sqrt(1 - pow(m3,2))*n3*
                  (k3*sqrt(1 - pow(m3,2))*n3 + 
                    k4*sqrt(1 - pow(m4,2))*n4) + 
                 k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                  (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                    k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
             (pow(k3,2)*pow(m3,2) + 
               pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + 
               pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
            (2*(k3*k4*m3*m4 + 
                 k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*n3*
                  n4 + k3*k4*sqrt(1 - pow(m3,2))*
                  sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*
                  sqrt(1 - pow(n4,2)))*
               (pow(k3*m3 + k4*m4,2) + 
                 pow(k3*sqrt(1 - pow(m3,2))*n3 + 
                   k4*sqrt(1 - pow(m4,2))*n4,2) + 
                 pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                   k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
             ((pow(k3,2)*pow(m3,2) + 
                 pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + 
                 pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
               (pow(k4,2)*pow(m4,2) + 
                 pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + 
                 pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/
        (7.*(pow(k1*m1 + k2*m2,2) + 
            pow(k1*sqrt(1 - pow(m1,2))*n1 + 
              k2*sqrt(1 - pow(m2,2))*n2,2) + 
            pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
              k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2))*
          (pow(k3*m3 + k4*m4,2) + 
            pow(k3*sqrt(1 - pow(m3,2))*n3 + 
              k4*sqrt(1 - pow(m4,2))*n4,2) + 
            pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
              k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))))/231. + 
  ((9*(k1*m1*(k1*m1 + k2*m2 + k3*m3 + k4*m4) + 
          k1*sqrt(1 - pow(m1,2))*n1*
           (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + 
             k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) \
+ k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*
           (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
             k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
             k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
             k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
        ((((3*(k2*m2*(k2*m2 + k3*m3) + 
                    k2*sqrt(1 - pow(m2,2))*n2*
                     (k2*sqrt(1 - pow(m2,2))*n2 + 
                       k3*sqrt(1 - pow(m3,2))*n3) + 
                    k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                     (k2*sqrt(1 - pow(m2,2))*
                       sqrt(1 - pow(n2,2)) + 
                       k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))\
)))/(pow(k2,2)*pow(m2,2) + 
                  pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + 
                  pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
               (2*(k2*k3*m2*m3 + 
                    k2*k3*sqrt(1 - pow(m2,2))*
                     sqrt(1 - pow(m3,2))*n2*n3 + 
                    k2*k3*sqrt(1 - pow(m2,2))*
                     sqrt(1 - pow(m3,2))*sqrt(1 - pow(n2,2))*
                     sqrt(1 - pow(n3,2)))*
                  (pow(k2*m2 + k3*m3,2) + 
                    pow(k2*sqrt(1 - pow(m2,2))*n2 + 
                      k3*sqrt(1 - pow(m3,2))*n3,2) + 
                    pow(k2*sqrt(1 - pow(m2,2))*
                       sqrt(1 - pow(n2,2)) + 
                      k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),
                     2)))/
                ((pow(k2,2)*pow(m2,2) + 
                    pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + 
                    pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                  (pow(k3,2)*pow(m3,2) + 
                    pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + 
                    pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))))*
             ((7*((k2*m2 + k3*m3)*(k2*m2 + k3*m3 + k4*m4) + 
                    (k2*sqrt(1 - pow(m2,2))*n2 + 
                       k3*sqrt(1 - pow(m3,2))*n3)*
                     (k2*sqrt(1 - pow(m2,2))*n2 + 
                       k3*sqrt(1 - pow(m3,2))*n3 + 
                       k4*sqrt(1 - pow(m4,2))*n4) + 
                    (k2*sqrt(1 - pow(m2,2))*
                       sqrt(1 - pow(n2,2)) + 
                       k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))\
)*(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                       k3*sqrt(1 - pow(m3,2))*
                       sqrt(1 - pow(n3,2)) + 
                       k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))\
))/(pow(k2*m2 + k3*m3,2) + 
                  pow(k2*sqrt(1 - pow(m2,2))*n2 + 
                    k3*sqrt(1 - pow(m3,2))*n3,2) + 
                  pow(k2*sqrt(1 - pow(m2,2))*
                     sqrt(1 - pow(n2,2)) + 
                    k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)) \
+ ((k4*(k2*m2 + k3*m3)*m4 + k4*sqrt(1 - pow(m4,2))*
                     (k2*sqrt(1 - pow(m2,2))*n2 + 
                       k3*sqrt(1 - pow(m3,2))*n3)*n4 + 
                    k4*sqrt(1 - pow(m4,2))*
                     (k2*sqrt(1 - pow(m2,2))*
                       sqrt(1 - pow(n2,2)) + 
                       k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))\
)*sqrt(1 - pow(n4,2)))*(pow(k2*m2 + k3*m3 + k4*m4,2) + 
                    pow(k2*sqrt(1 - pow(m2,2))*n2 + 
                      k3*sqrt(1 - pow(m3,2))*n3 + 
                      k4*sqrt(1 - pow(m4,2))*n4,2) + 
                    pow(k2*sqrt(1 - pow(m2,2))*
                       sqrt(1 - pow(n2,2)) + 
                      k3*sqrt(1 - pow(m3,2))*
                       sqrt(1 - pow(n3,2)) + 
                      k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)\
))/((pow(k2*m2 + k3*m3,2) + 
                    pow(k2*sqrt(1 - pow(m2,2))*n2 + 
                      k3*sqrt(1 - pow(m3,2))*n3,2) + 
                    pow(k2*sqrt(1 - pow(m2,2))*
                       sqrt(1 - pow(n2,2)) + 
                      k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)\
)*(pow(k4,2)*pow(m4,2) + 
                    pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + 
                    pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/
           126. + (((k2*m2*(k2*m2 + k3*m3 + k4*m4) + 
                  k2*sqrt(1 - pow(m2,2))*n2*
                   (k2*sqrt(1 - pow(m2,2))*n2 + 
                     k3*sqrt(1 - pow(m3,2))*n3 + 
                     k4*sqrt(1 - pow(m4,2))*n4) + 
                  k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                   (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                     k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                     k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
                ((5*(k3*m3*(k3*m3 + k4*m4) + 
                       k3*sqrt(1 - pow(m3,2))*n3*
                       (k3*sqrt(1 - pow(m3,2))*n3 + 
                       k4*sqrt(1 - pow(m4,2))*n4) + 
                       k3*sqrt(1 - pow(m3,2))*
                        sqrt(1 - pow(n3,2))*
                        (k3*sqrt(1 - pow(m3,2))*
                       sqrt(1 - pow(n3,2)) + 
                        k4*sqrt(1 - pow(m4,2))*
                        sqrt(1 - pow(n4,2)))))/
                   (pow(k3,2)*pow(m3,2) + 
                     pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + 
                     pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) \
+ ((k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*
                       n3*n4 + 
                       k3*k4*sqrt(1 - pow(m3,2))*
                        sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*
                        sqrt(1 - pow(n4,2)))*
                     (pow(k3*m3 + k4*m4,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*n3 + 
                       k4*sqrt(1 - pow(m4,2))*n4,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*
                       sqrt(1 - pow(n3,2)) + 
                        k4*sqrt(1 - pow(m4,2))*
                        sqrt(1 - pow(n4,2)),2)))/
                   ((pow(k3,2)*pow(m3,2) + 
                       pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + 
                       pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))\
)*(pow(k4,2)*pow(m4,2) + 
                       pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + 
                       pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))\
)))/(pow(k2,2)*pow(m2,2) + 
                pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + 
                pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
             ((k2*m2*(k3*m3 + k4*m4) + 
                  k2*sqrt(1 - pow(m2,2))*n2*
                   (k3*sqrt(1 - pow(m3,2))*n3 + 
                     k4*sqrt(1 - pow(m4,2))*n4) + 
                  k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                   (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                     k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
                (pow(k2*m2 + k3*m3 + k4*m4,2) + 
                  pow(k2*sqrt(1 - pow(m2,2))*n2 + 
                    k3*sqrt(1 - pow(m3,2))*n3 + 
                    k4*sqrt(1 - pow(m4,2))*n4,2) + 
                  pow(k2*sqrt(1 - pow(m2,2))*
                     sqrt(1 - pow(n2,2)) + 
                    k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                    k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
                ((3*(k3*m3*(k3*m3 + k4*m4) + 
                       k3*sqrt(1 - pow(m3,2))*n3*
                        (k3*sqrt(1 - pow(m3,2))*n3 + 
                        k4*sqrt(1 - pow(m4,2))*n4) + 
                       k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                        (k3*sqrt(1 - pow(m3,2))*
                        sqrt(1 - pow(n3,2)) + 
                        k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))\
)))/(pow(k3,2)*pow(m3,2) + 
                     pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + 
                     pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                  (2*(k3*k4*m3*m4 + 
                       k3*k4*sqrt(1 - pow(m3,2))*
                        sqrt(1 - pow(m4,2))*n3*n4 + 
                       k3*k4*sqrt(1 - pow(m3,2))*
                        sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*
                        sqrt(1 - pow(n4,2)))*
                     (pow(k3*m3 + k4*m4,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*n3 + 
                        k4*sqrt(1 - pow(m4,2))*n4,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*
                        sqrt(1 - pow(n3,2)) + 
                        k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),
                        2)))/
                   ((pow(k3,2)*pow(m3,2) + 
                       pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + 
                       pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                     (pow(k4,2)*pow(m4,2) + 
                       pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + 
                       pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))\
))/(7.*(pow(k2,2)*pow(m2,2) + 
                  pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + 
                  pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                (pow(k3*m3 + k4*m4,2) + 
                  pow(k3*sqrt(1 - pow(m3,2))*n3 + 
                    k4*sqrt(1 - pow(m4,2))*n4,2) + 
                  pow(k3*sqrt(1 - pow(m3,2))*
                     sqrt(1 - pow(n3,2)) + 
                    k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))))/
           18.))/
      (pow(k1,2)*pow(m1,2) + 
        pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + 
        pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2))) + 
     ((k1*m1*(k2*m2 + k3*m3 + k4*m4) + 
          k1*sqrt(1 - pow(m1,2))*n1*
           (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + 
             k4*sqrt(1 - pow(m4,2))*n4) + 
          k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*
           (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
             k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
             k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
        (pow(k1*m1 + k2*m2 + k3*m3 + k4*m4,2) + 
          pow(k1*sqrt(1 - pow(m1,2))*n1 + 
            k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + 
            k4*sqrt(1 - pow(m4,2))*n4,2) + 
          pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
            k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
            k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
            k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
        ((((3*(k2*m2*(k2*m2 + k3*m3) + 
                    k2*sqrt(1 - pow(m2,2))*n2*
                     (k2*sqrt(1 - pow(m2,2))*n2 + 
                       k3*sqrt(1 - pow(m3,2))*n3) + 
                    k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                     (k2*sqrt(1 - pow(m2,2))*
                       sqrt(1 - pow(n2,2)) + 
                       k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))\
))/(pow(k2,2)*pow(m2,2) + 
                  pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + 
                  pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
               (2*(k2*k3*m2*m3 + 
                    k2*k3*sqrt(1 - pow(m2,2))*
                     sqrt(1 - pow(m3,2))*n2*n3 + 
                    k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*
                     sqrt(1 - pow(n2,2))*sqrt(1 - pow(n3,2)))*
                  (pow(k2*m2 + k3*m3,2) + 
                    pow(k2*sqrt(1 - pow(m2,2))*n2 + 
                      k3*sqrt(1 - pow(m3,2))*n3,2) + 
                    pow(k2*sqrt(1 - pow(m2,2))*
                       sqrt(1 - pow(n2,2)) + 
                      k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)\
))/((pow(k2,2)*pow(m2,2) + 
                    pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + 
                    pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                  (pow(k3,2)*pow(m3,2) + 
                    pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + 
                    pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))))*
             ((3*((k2*m2 + k3*m3)*(k2*m2 + k3*m3 + k4*m4) + 
                    (k2*sqrt(1 - pow(m2,2))*n2 + 
                       k3*sqrt(1 - pow(m3,2))*n3)*
                     (k2*sqrt(1 - pow(m2,2))*n2 + 
                       k3*sqrt(1 - pow(m3,2))*n3 + 
                       k4*sqrt(1 - pow(m4,2))*n4) + 
                    (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                       k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
                     (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                       k3*sqrt(1 - pow(m3,2))*
                       sqrt(1 - pow(n3,2)) + 
                       k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))\
)/(pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + 
                    k3*sqrt(1 - pow(m3,2))*n3,2) + 
                  pow(k2*sqrt(1 - pow(m2,2))*
                     sqrt(1 - pow(n2,2)) + 
                    k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)) + 
               (3*(k4*(k2*m2 + k3*m3)*m4 + 
                    k4*sqrt(1 - pow(m4,2))*
                     (k2*sqrt(1 - pow(m2,2))*n2 + 
                       k3*sqrt(1 - pow(m3,2))*n3)*n4 + 
                    k4*sqrt(1 - pow(m4,2))*
                     (k2*sqrt(1 - pow(m2,2))*
                       sqrt(1 - pow(n2,2)) + 
                       k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
                     sqrt(1 - pow(n4,2)))*
                  (pow(k2*m2 + k3*m3 + k4*m4,2) + 
                    pow(k2*sqrt(1 - pow(m2,2))*n2 + 
                      k3*sqrt(1 - pow(m3,2))*n3 + 
                      k4*sqrt(1 - pow(m4,2))*n4,2) + 
                    pow(k2*sqrt(1 - pow(m2,2))*
                       sqrt(1 - pow(n2,2)) + 
                      k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                      k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))\
)/((pow(k2*m2 + k3*m3,2) + 
                    pow(k2*sqrt(1 - pow(m2,2))*n2 + 
                      k3*sqrt(1 - pow(m3,2))*n3,2) + 
                    pow(k2*sqrt(1 - pow(m2,2))*
                       sqrt(1 - pow(n2,2)) + 
                      k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))*
                  (pow(k4,2)*pow(m4,2) + 
                    pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + 
                    pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/
           126. + ((3*(k2*m2*(k2*m2 + k3*m3 + k4*m4) + 
                  k2*sqrt(1 - pow(m2,2))*n2*
                   (k2*sqrt(1 - pow(m2,2))*n2 + 
                     k3*sqrt(1 - pow(m3,2))*n3 + 
                     k4*sqrt(1 - pow(m4,2))*n4) + 
                  k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                   (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                     k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                     k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
                ((5*(k3*m3*(k3*m3 + k4*m4) + 
                       k3*sqrt(1 - pow(m3,2))*n3*
                        (k3*sqrt(1 - pow(m3,2))*n3 + 
                        k4*sqrt(1 - pow(m4,2))*n4) + 
                       k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                        (k3*sqrt(1 - pow(m3,2))*
                        sqrt(1 - pow(n3,2)) + 
                        k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))\
)))/(pow(k3,2)*pow(m3,2) + 
                     pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + 
                     pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                  ((k3*k4*m3*m4 + 
                       k3*k4*sqrt(1 - pow(m3,2))*
                        sqrt(1 - pow(m4,2))*n3*n4 + 
                       k3*k4*sqrt(1 - pow(m3,2))*
                        sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*
                        sqrt(1 - pow(n4,2)))*
                     (pow(k3*m3 + k4*m4,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*n3 + 
                        k4*sqrt(1 - pow(m4,2))*n4,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*
                        sqrt(1 - pow(n3,2)) + 
                        k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),
                        2)))/
                   ((pow(k3,2)*pow(m3,2) + 
                       pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + 
                       pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                     (pow(k4,2)*pow(m4,2) + 
                       pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + 
                       pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))\
))/(7.*(pow(k2,2)*pow(m2,2) + 
                  pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + 
                  pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))) + 
             (3*(k2*m2*(k3*m3 + k4*m4) + 
                  k2*sqrt(1 - pow(m2,2))*n2*
                   (k3*sqrt(1 - pow(m3,2))*n3 + 
                     k4*sqrt(1 - pow(m4,2))*n4) + 
                  k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                   (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                     k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
                (pow(k2*m2 + k3*m3 + k4*m4,2) + 
                  pow(k2*sqrt(1 - pow(m2,2))*n2 + 
                    k3*sqrt(1 - pow(m3,2))*n3 + 
                    k4*sqrt(1 - pow(m4,2))*n4,2) + 
                  pow(k2*sqrt(1 - pow(m2,2))*
                     sqrt(1 - pow(n2,2)) + 
                    k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                    k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
                ((3*(k3*m3*(k3*m3 + k4*m4) + 
                       k3*sqrt(1 - pow(m3,2))*n3*
                        (k3*sqrt(1 - pow(m3,2))*n3 + 
                        k4*sqrt(1 - pow(m4,2))*n4) + 
                       k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                        (k3*sqrt(1 - pow(m3,2))*
                        sqrt(1 - pow(n3,2)) + 
                        k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))\
))/(pow(k3,2)*pow(m3,2) + 
                     pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + 
                     pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                  (2*(k3*k4*m3*m4 + 
                       k3*k4*sqrt(1 - pow(m3,2))*
                        sqrt(1 - pow(m4,2))*n3*n4 + 
                       k3*k4*sqrt(1 - pow(m3,2))*
                        sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*
                        sqrt(1 - pow(n4,2)))*
                     (pow(k3*m3 + k4*m4,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*n3 + 
                        k4*sqrt(1 - pow(m4,2))*n4,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*
                        sqrt(1 - pow(n3,2)) + 
                        k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),
                        2)))/
                   ((pow(k3,2)*pow(m3,2) + 
                       pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + 
                       pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                     (pow(k4,2)*pow(m4,2) + 
                       pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + 
                       pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))))\
)/(7.*(pow(k2,2)*pow(m2,2) + 
                  pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + 
                  pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                (pow(k3*m3 + k4*m4,2) + 
                  pow(k3*sqrt(1 - pow(m3,2))*n3 + 
                    k4*sqrt(1 - pow(m4,2))*n4,2) + 
                  pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                    k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))))/
           18.))/
      ((pow(k1,2)*pow(m1,2) + 
          pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + 
          pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2)))*
        (pow(k2*m2 + k3*m3 + k4*m4,2) + 
          pow(k2*sqrt(1 - pow(m2,2))*n2 + 
            k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
          pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
            k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
            k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))))/33. ;
return (temp3); 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// F5 functions ///////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Postion of k at fisrt argument, total 6 terms //////////////////////////
float F5(float k1, float k2, float k3, float k4, float k5, float m1, float m2, float m3, float m4, float m5, float n1, float n2, float n3, float n4, float n5) /* F5(k,q,-q,p,-p)*/
{ float temp2;
temp2 = (((((((3*(k1*m1*(k1*m1 + k2*m2) + k1*sqrt(1 - pow(m1,2))*n1*(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2) + 
                       k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*
                        (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))))/
                   (pow(k1,2)*pow(m1,2) + pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2))) + 
                  (2*(k1*k2*m1*m2 + k1*k2*sqrt(1 - pow(m1,2))*sqrt(1 - pow(m2,2))*n1*n2 + 
                       k1*k2*sqrt(1 - pow(m1,2))*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n1,2))*sqrt(1 - pow(n2,2)))*
                     (pow(k1*m1 + k2*m2,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2,2) + 
                       pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2)))/
                   ((pow(k1,2)*pow(m1,2) + pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2)))*
                     (pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))))*
                ((3*((k1*m1 + k2*m2)*(k1*m1 + k2*m2 + k3*m3) + 
                       (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2)*
                        (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3) + 
                       (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))*
                        (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                          k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))))/
                   (pow(k1*m1 + k2*m2,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2,2) + 
                     pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2)) + 
                  (3*(k3*(k1*m1 + k2*m2)*m3 + k3*sqrt(1 - pow(m3,2))*(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2)*n3 + 
                       k3*sqrt(1 - pow(m3,2))*(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))*
                        sqrt(1 - pow(n3,2)))*(pow(k1*m1 + k2*m2 + k3*m3,2) + 
                       pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                       pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                         k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)))/
                   ((pow(k1*m1 + k2*m2,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2,2) + 
                       pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2))*
                     (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))))))/126. + 
             ((3*(k1*m1*(k1*m1 + k2*m2 + k3*m3) + k1*sqrt(1 - pow(m1,2))*n1*
                      (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3) + 
                     k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*
                      (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                        k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))))*
                   ((5*(k2*m2*(k2*m2 + k3*m3) + k2*sqrt(1 - pow(m2,2))*n2*(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3) + 
                          k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                           (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))))/
                      (pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
                     ((k2*k3*m2*m3 + k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*n2*n3 + 
                          k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n2,2))*sqrt(1 - pow(n3,2)))*
                        (pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)))/
                      ((pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                        (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))))))/
                 (7.*(pow(k1,2)*pow(m1,2) + pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2)))) + 
                (3*(k1*m1*(k2*m2 + k3*m3) + k1*sqrt(1 - pow(m1,2))*n1*(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3) + 
                     k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*
                      (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))))*
                   (pow(k1*m1 + k2*m2 + k3*m3,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,
                      2) + pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                       k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))*
                   ((3*(k2*m2*(k2*m2 + k3*m3) + k2*sqrt(1 - pow(m2,2))*n2*(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3) + 
                          k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                           (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))))/
                      (pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
                     (2*(k2*k3*m2*m3 + k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*n2*n3 + 
                          k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n2,2))*sqrt(1 - pow(n3,2)))*
                        (pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)))/
                      ((pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                        (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))))))/
                 (7.*(pow(k1,2)*pow(m1,2) + pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2)))*
                   (pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))))/18.)*
           ((3*((k1*m1 + k2*m2 + k3*m3)*(k1*m1 + k2*m2 + k3*m3 + k4*m4) + 
                  (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3)*
                   (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                  (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                     k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
                   (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                     k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
              (pow(k1*m1 + k2*m2 + k3*m3,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                  k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)) + 
             (4*(k4*(k1*m1 + k2*m2 + k3*m3)*m4 + k4*sqrt(1 - pow(m4,2))*
                   (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3)*n4 + 
                  k4*sqrt(1 - pow(m4,2))*(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                     k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*sqrt(1 - pow(n4,2)))*
                (pow(k1*m1 + k2*m2 + k3*m3 + k4*m4,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + 
                    k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                  pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                    k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
              ((pow(k1*m1 + k2*m2 + k3*m3,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                  pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                    k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))*
                (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/33. + 
        (((3*(k1*m1*(k1*m1 + k2*m2) + k1*sqrt(1 - pow(m1,2))*n1*(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2) + 
                  k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*
                   (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))))/
              (pow(k1,2)*pow(m1,2) + pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2))) + 
             (2*(k1*k2*m1*m2 + k1*k2*sqrt(1 - pow(m1,2))*sqrt(1 - pow(m2,2))*n1*n2 + 
                  k1*k2*sqrt(1 - pow(m1,2))*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n1,2))*sqrt(1 - pow(n2,2)))*
                (pow(k1*m1 + k2*m2,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2,2) + 
                  pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2)))/
              ((pow(k1,2)*pow(m1,2) + pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2)))*
                (pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))))*
           ((3*((k1*m1 + k2*m2)*(k1*m1 + k2*m2 + k3*m3 + k4*m4) + 
                  (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2)*
                   (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                  (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))*
                   (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                     k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
                ((5*(k3*m3*(k3*m3 + k4*m4) + k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                       k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                        (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                   (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                  ((k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*n3*n4 + 
                       k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*sqrt(1 - pow(n4,2)))*
                     (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                   ((pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                     (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/
              (7.*(pow(k1*m1 + k2*m2,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2,2) + 
                  pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2))) + 
             (4*((k1*m1 + k2*m2)*(k3*m3 + k4*m4) + (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2)*
                   (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                  (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))*
                   (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
                (pow(k1*m1 + k2*m2 + k3*m3 + k4*m4,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + 
                    k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                  pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                    k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
                ((3*(k3*m3*(k3*m3 + k4*m4) + k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                       k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                        (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                   (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                  (2*(k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*n3*n4 + 
                       k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*sqrt(1 - pow(n4,2)))*
                     (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                   ((pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                     (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/
              (7.*(pow(k1*m1 + k2*m2,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2,2) + 
                  pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2))*
                (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                  pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))))/231. + 
        ((3*(k1*m1*(k1*m1 + k2*m2 + k3*m3 + k4*m4) + k1*sqrt(1 - pow(m1,2))*n1*
                 (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*
                 (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                   k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
              ((((3*(k2*m2*(k2*m2 + k3*m3) + k2*sqrt(1 - pow(m2,2))*n2*(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3) + 
                          k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                           (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))))/
                      (pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
                     (2*(k2*k3*m2*m3 + k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*n2*n3 + 
                          k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n2,2))*sqrt(1 - pow(n3,2)))*
                        (pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)))/
                      ((pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                        (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))))*
                   ((7*((k2*m2 + k3*m3)*(k2*m2 + k3*m3 + k4*m4) + 
                          (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3)*
                           (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                          (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
                           (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                             k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                      (pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                        pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)) + 
                     ((k4*(k2*m2 + k3*m3)*m4 + k4*sqrt(1 - pow(m4,2))*(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3)*n4 + 
                          k4*sqrt(1 - pow(m4,2))*(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
                           sqrt(1 - pow(n4,2)))*(pow(k2*m2 + k3*m3 + k4*m4,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                            k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                      ((pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))*
                        (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/126.
                  + (((k2*m2*(k2*m2 + k3*m3 + k4*m4) + k2*sqrt(1 - pow(m2,2))*n2*
                         (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                        k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                         (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                           k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
                      ((5*(k3*m3*(k3*m3 + k4*m4) + k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                             k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                              (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                         (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                        ((k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*n3*n4 + 
                             k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*sqrt(1 - pow(n4,2)))*
                           (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                         ((pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                           (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/
                    (pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
                   ((k2*m2*(k3*m3 + k4*m4) + k2*sqrt(1 - pow(m2,2))*n2*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                        k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                         (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
                      (pow(k2*m2 + k3*m3 + k4*m4,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,
                         2) + pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                          k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
                      ((3*(k3*m3*(k3*m3 + k4*m4) + k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                             k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                              (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                         (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                        (2*(k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*n3*n4 + 
                             k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*sqrt(1 - pow(n4,2)))*
                           (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                         ((pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                           (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/
                    (7.*(pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                      (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                        pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))))/18.))/
            (pow(k1,2)*pow(m1,2) + pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2))) + 
           (4*(k1*m1*(k2*m2 + k3*m3 + k4*m4) + k1*sqrt(1 - pow(m1,2))*n1*
                 (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*
                 (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                   k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
              (pow(k1*m1 + k2*m2 + k3*m3 + k4*m4,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + 
                  k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                  k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
              ((((3*(k2*m2*(k2*m2 + k3*m3) + k2*sqrt(1 - pow(m2,2))*n2*(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3) + 
                          k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                           (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))))/
                      (pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
                     (2*(k2*k3*m2*m3 + k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*n2*n3 + 
                          k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n2,2))*sqrt(1 - pow(n3,2)))*
                        (pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)))/
                      ((pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                        (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))))*
                   ((3*((k2*m2 + k3*m3)*(k2*m2 + k3*m3 + k4*m4) + 
                          (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3)*
                           (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                          (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
                           (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                             k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                      (pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                        pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)) + 
                     (3*(k4*(k2*m2 + k3*m3)*m4 + k4*sqrt(1 - pow(m4,2))*(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3)*n4 + 
                          k4*sqrt(1 - pow(m4,2))*(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
                           sqrt(1 - pow(n4,2)))*(pow(k2*m2 + k3*m3 + k4*m4,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                            k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                      ((pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))*
                        (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/126.
                  + ((3*(k2*m2*(k2*m2 + k3*m3 + k4*m4) + k2*sqrt(1 - pow(m2,2))*n2*
                         (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                        k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                         (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                           k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
                      ((5*(k3*m3*(k3*m3 + k4*m4) + k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                             k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                              (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                         (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                        ((k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*n3*n4 + 
                             k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*sqrt(1 - pow(n4,2)))*
                           (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                         ((pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                           (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/
                    (7.*(pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))) + 
                   (3*(k2*m2*(k3*m3 + k4*m4) + k2*sqrt(1 - pow(m2,2))*n2*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                        k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                         (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
                      (pow(k2*m2 + k3*m3 + k4*m4,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,
                         2) + pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                          k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
                      ((3*(k3*m3*(k3*m3 + k4*m4) + k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                             k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                              (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                         (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                        (2*(k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*n3*n4 + 
                             k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*sqrt(1 - pow(n4,2)))*
                           (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                         ((pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                           (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/
                    (7.*(pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                      (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                        pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))))/18.))/
            ((pow(k1,2)*pow(m1,2) + pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2)))*
              (pow(k2*m2 + k3*m3 + k4*m4,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                  k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))))/33.)*
      ((11*((k1*m1 + k2*m2 + k3*m3 + k4*m4)*(k1*m1 + k2*m2 + k3*m3 + k4*m4 + k5*m5) + 
             (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*
              (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + 
                k5*sqrt(1 - pow(m5,2))*n5) + (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*
              (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
         (pow(k1*m1 + k2*m2 + k3*m3 + k4*m4,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + 
             k4*sqrt(1 - pow(m4,2))*n4,2) + pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
             k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)) + 
        ((k5*(k1*m1 + k2*m2 + k3*m3 + k4*m4)*m5 + k5*sqrt(1 - pow(m5,2))*
              (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*n5 + 
             k5*sqrt(1 - pow(m5,2))*(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*sqrt(1 - pow(n5,2)))*
           (pow(k1*m1 + k2*m2 + k3*m3 + k4*m4 + k5*m5,2) + 
             pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + 
               k5*sqrt(1 - pow(m5,2))*n5,2) + pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
               k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
               k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
         ((pow(k1*m1 + k2*m2 + k3*m3 + k4*m4,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + 
               k4*sqrt(1 - pow(m4,2))*n4,2) + pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
               k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
               k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
           (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/52. + 
   (((((3*(k1*m1*(k1*m1 + k2*m2) + k1*sqrt(1 - pow(m1,2))*n1*(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2) + 
                  k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*
                   (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))))/
              (pow(k1,2)*pow(m1,2) + pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2))) + 
             (2*(k1*k2*m1*m2 + k1*k2*sqrt(1 - pow(m1,2))*sqrt(1 - pow(m2,2))*n1*n2 + 
                  k1*k2*sqrt(1 - pow(m1,2))*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n1,2))*sqrt(1 - pow(n2,2)))*
                (pow(k1*m1 + k2*m2,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2,2) + 
                  pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2)))/
              ((pow(k1,2)*pow(m1,2) + pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2)))*
                (pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))))*
           ((3*((k1*m1 + k2*m2)*(k1*m1 + k2*m2 + k3*m3) + (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2)*
                   (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3) + 
                  (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))*
                   (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                     k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))))/
              (pow(k1*m1 + k2*m2,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2,2) + 
                pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2)) + 
             (3*(k3*(k1*m1 + k2*m2)*m3 + k3*sqrt(1 - pow(m3,2))*(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2)*n3 + 
                  k3*sqrt(1 - pow(m3,2))*(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))*
                   sqrt(1 - pow(n3,2)))*(pow(k1*m1 + k2*m2 + k3*m3,2) + 
                  pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                  pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                    k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)))/
              ((pow(k1*m1 + k2*m2,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2,2) + 
                  pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2))*
                (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))))))/126. + 
        ((3*(k1*m1*(k1*m1 + k2*m2 + k3*m3) + k1*sqrt(1 - pow(m1,2))*n1*
                 (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3) + 
                k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*
                 (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                   k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))))*
              ((5*(k2*m2*(k2*m2 + k3*m3) + k2*sqrt(1 - pow(m2,2))*n2*(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3) + 
                     k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                      (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))))/
                 (pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
                ((k2*k3*m2*m3 + k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*n2*n3 + 
                     k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n2,2))*sqrt(1 - pow(n3,2)))*
                   (pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)))/
                 ((pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                   (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))))))/
            (7.*(pow(k1,2)*pow(m1,2) + pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2)))) + 
           (3*(k1*m1*(k2*m2 + k3*m3) + k1*sqrt(1 - pow(m1,2))*n1*(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3) + 
                k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*
                 (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))))*
              (pow(k1*m1 + k2*m2 + k3*m3,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                  k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))*
              ((3*(k2*m2*(k2*m2 + k3*m3) + k2*sqrt(1 - pow(m2,2))*n2*(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3) + 
                     k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                      (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))))/
                 (pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
                (2*(k2*k3*m2*m3 + k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*n2*n3 + 
                     k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n2,2))*sqrt(1 - pow(n3,2)))*
                   (pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)))/
                 ((pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                   (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))))))/
            (7.*(pow(k1,2)*pow(m1,2) + pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2)))*
              (pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))))/18.)*
      ((11*((k1*m1 + k2*m2 + k3*m3)*(k1*m1 + k2*m2 + k3*m3 + k4*m4 + k5*m5) + 
             (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3)*
              (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + 
                k5*sqrt(1 - pow(m5,2))*n5) + (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
              (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
           ((5*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                  k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                   (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
              (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
             ((k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                  k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                  pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
              ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/
         (7.*(pow(k1*m1 + k2*m2 + k3*m3,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
             pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
               k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))) + 
        (((k1*m1 + k2*m2 + k3*m3)*(k4*m4 + k5*m5) + (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3)*
              (k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
             (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
              (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
           (pow(k1*m1 + k2*m2 + k3*m3 + k4*m4 + k5*m5,2) + 
             pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + 
               k5*sqrt(1 - pow(m5,2))*n5,2) + pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
               k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
               k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))*
           ((3*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                  k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                   (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
              (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
             (2*(k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                  k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                  pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
              ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/
         (7.*(pow(k1*m1 + k2*m2 + k3*m3,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
             pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
               k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))*
           (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
             pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))))/52. + 
   (((3*(k1*m1*(k1*m1 + k2*m2) + k1*sqrt(1 - pow(m1,2))*n1*(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2) + 
             k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
                k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))))/
         (pow(k1,2)*pow(m1,2) + pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2))) + 
        (2*(k1*k2*m1*m2 + k1*k2*sqrt(1 - pow(m1,2))*sqrt(1 - pow(m2,2))*n1*n2 + 
             k1*k2*sqrt(1 - pow(m1,2))*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n1,2))*sqrt(1 - pow(n2,2)))*
           (pow(k1*m1 + k2*m2,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2,2) + 
             pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2)))/
         ((pow(k1,2)*pow(m1,2) + pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2)))*
           (pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))))*
      ((11*((k1*m1 + k2*m2)*(k1*m1 + k2*m2 + k3*m3 + k4*m4 + k5*m5) + 
             (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2)*
              (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + 
                k5*sqrt(1 - pow(m5,2))*n5) + (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))*
              (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
                k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
           ((((3*(k3*m3*(k3*m3 + k4*m4) + k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                       k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                        (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                   (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                  (2*(k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*n3*n4 + 
                       k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*sqrt(1 - pow(n4,2)))*
                     (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                   ((pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                     (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))))*
                ((7*((k3*m3 + k4*m4)*(k3*m3 + k4*m4 + k5*m5) + 
                       (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*
                        (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                       (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*
                        (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                          k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                   (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                     pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)) + 
                  ((k5*(k3*m3 + k4*m4)*m5 + k5*sqrt(1 - pow(m5,2))*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*n5 + 
                       k5*sqrt(1 - pow(m5,2))*(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*
                        sqrt(1 - pow(n5,2)))*(pow(k3*m3 + k4*m4 + k5*m5,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                         k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                   ((pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
                     (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/126. + 
             (((k3*m3*(k3*m3 + k4*m4 + k5*m5) + k3*sqrt(1 - pow(m3,2))*n3*
                      (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                     k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                      (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                        k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                   ((5*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                          k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                           (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                      (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
                     ((k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                          k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                        (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                          pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                      ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                        (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/
                 (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                ((k3*m3*(k4*m4 + k5*m5) + k3*sqrt(1 - pow(m3,2))*n3*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                     k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                      (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                   (pow(k3*m3 + k4*m4 + k5*m5,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,
                      2) + pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                       k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))*
                   ((3*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                          k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                           (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                      (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
                     (2*(k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                          k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                        (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                          pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                      ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                        (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/
                 (7.*(pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                   (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                     pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))))/18.))/
         (pow(k1*m1 + k2*m2,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2,2) + 
           pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2)) + 
        (((k1*m1 + k2*m2)*(k3*m3 + k4*m4 + k5*m5) + (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2)*
              (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
             (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)))*
              (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
           (pow(k1*m1 + k2*m2 + k3*m3 + k4*m4 + k5*m5,2) + 
             pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + 
               k5*sqrt(1 - pow(m5,2))*n5,2) + pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + 
               k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
               k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))*
           ((((3*(k3*m3*(k3*m3 + k4*m4) + k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                       k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                        (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                   (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                  (2*(k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*n3*n4 + 
                       k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*sqrt(1 - pow(n4,2)))*
                     (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                   ((pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                     (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))))*
                ((3*((k3*m3 + k4*m4)*(k3*m3 + k4*m4 + k5*m5) + 
                       (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*
                        (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                       (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*
                        (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                          k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                   (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                     pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)) + 
                  (3*(k5*(k3*m3 + k4*m4)*m5 + k5*sqrt(1 - pow(m5,2))*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*n5 + 
                       k5*sqrt(1 - pow(m5,2))*(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*
                        sqrt(1 - pow(n5,2)))*(pow(k3*m3 + k4*m4 + k5*m5,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                         k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                   ((pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                       pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
                     (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/126. + 
             ((3*(k3*m3*(k3*m3 + k4*m4 + k5*m5) + k3*sqrt(1 - pow(m3,2))*n3*
                      (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                     k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                      (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                        k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                   ((5*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                          k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                           (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                      (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
                     ((k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                          k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                        (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                          pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                      ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                        (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/
                 (7.*(pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))) + 
                (3*(k3*m3*(k4*m4 + k5*m5) + k3*sqrt(1 - pow(m3,2))*n3*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                     k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                      (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                   (pow(k3*m3 + k4*m4 + k5*m5,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,
                      2) + pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                       k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))*
                   ((3*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                          k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                           (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                      (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
                     (2*(k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                          k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                        (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                          pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                      ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                        (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/
                 (7.*(pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                   (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                     pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))))/18.))/
         ((pow(k1*m1 + k2*m2,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2,2) + 
             pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)),2))*
           (pow(k3*m3 + k4*m4 + k5*m5,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
               k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))))/364. + 
   ((11*(k1*m1*(k1*m1 + k2*m2 + k3*m3 + k4*m4 + k5*m5) + k1*sqrt(1 - pow(m1,2))*n1*
            (k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + 
              k5*sqrt(1 - pow(m5,2))*n5) + k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*
            (k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
              k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
              k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
         ((((((3*(k2*m2*(k2*m2 + k3*m3) + k2*sqrt(1 - pow(m2,2))*n2*(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3) + 
                          k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                           (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))))/
                      (pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
                     (2*(k2*k3*m2*m3 + k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*n2*n3 + 
                          k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n2,2))*sqrt(1 - pow(n3,2)))*
                        (pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)))/
                      ((pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                        (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))))*
                   ((3*((k2*m2 + k3*m3)*(k2*m2 + k3*m3 + k4*m4) + 
                          (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3)*
                           (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                          (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
                           (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                             k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                      (pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                        pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)) + 
                     (3*(k4*(k2*m2 + k3*m3)*m4 + k4*sqrt(1 - pow(m4,2))*(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3)*n4 + 
                          k4*sqrt(1 - pow(m4,2))*(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
                           sqrt(1 - pow(n4,2)))*(pow(k2*m2 + k3*m3 + k4*m4,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                            k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                      ((pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))*
                        (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/126.
                  + ((3*(k2*m2*(k2*m2 + k3*m3 + k4*m4) + k2*sqrt(1 - pow(m2,2))*n2*
                         (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                        k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                         (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                           k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
                      ((5*(k3*m3*(k3*m3 + k4*m4) + k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                             k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                              (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                         (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                        ((k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*n3*n4 + 
                             k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*sqrt(1 - pow(n4,2)))*
                           (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                         ((pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                           (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/
                    (7.*(pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))) + 
                   (3*(k2*m2*(k3*m3 + k4*m4) + k2*sqrt(1 - pow(m2,2))*n2*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                        k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                         (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
                      (pow(k2*m2 + k3*m3 + k4*m4,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,
                         2) + pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                          k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
                      ((3*(k3*m3*(k3*m3 + k4*m4) + k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                             k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                              (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                         (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                        (2*(k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*n3*n4 + 
                             k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*sqrt(1 - pow(n4,2)))*
                           (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                         ((pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                           (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/
                    (7.*(pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                      (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                        pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))))/18.)*
              ((9*((k2*m2 + k3*m3 + k4*m4)*(k2*m2 + k3*m3 + k4*m4 + k5*m5) + 
                     (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*
                      (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                     (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                        k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*
                      (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                        k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                 (pow(k2*m2 + k3*m3 + k4*m4,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                   pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                     k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)) + 
                ((k5*(k2*m2 + k3*m3 + k4*m4)*m5 + k5*sqrt(1 - pow(m5,2))*
                      (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*n5 + 
                     k5*sqrt(1 - pow(m5,2))*(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                        k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*sqrt(1 - pow(n5,2)))*
                   (pow(k2*m2 + k3*m3 + k4*m4 + k5*m5,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                       k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                 ((pow(k2*m2 + k3*m3 + k4*m4,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                       k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
                   (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/33. + 
           (((3*(k2*m2*(k2*m2 + k3*m3) + k2*sqrt(1 - pow(m2,2))*n2*(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3) + 
                     k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                      (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))))/
                 (pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
                (2*(k2*k3*m2*m3 + k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*n2*n3 + 
                     k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n2,2))*sqrt(1 - pow(n3,2)))*
                   (pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)))/
                 ((pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                   (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))))*
              ((9*((k2*m2 + k3*m3)*(k2*m2 + k3*m3 + k4*m4 + k5*m5) + 
                     (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3)*
                      (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                     (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
                      (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                        k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                   ((5*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                          k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                           (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                      (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
                     ((k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                          k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                        (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                          pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                      ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                        (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/
                 (7.*(pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))) + 
                (((k2*m2 + k3*m3)*(k4*m4 + k5*m5) + (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3)*
                      (k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                     (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
                      (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                   (pow(k2*m2 + k3*m3 + k4*m4 + k5*m5,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                       k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))*
                   ((3*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                          k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                           (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                      (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
                     (2*(k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                          k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                        (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                          pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                      ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                        (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/
                 (7.*(pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))*
                   (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                     pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))))/231. + 
           ((9*(k2*m2*(k2*m2 + k3*m3 + k4*m4 + k5*m5) + k2*sqrt(1 - pow(m2,2))*n2*
                    (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                   k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                    (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                      k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                 ((((3*(k3*m3*(k3*m3 + k4*m4) + k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                             k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                              (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                         (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                        (2*(k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*n3*n4 + 
                             k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*sqrt(1 - pow(n4,2)))*
                           (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                         ((pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                           (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))))*
                      ((7*((k3*m3 + k4*m4)*(k3*m3 + k4*m4 + k5*m5) + 
                             (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*
                              (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                             (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*
                              (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                                k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                         (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                           pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)) + 
                        ((k5*(k3*m3 + k4*m4)*m5 + k5*sqrt(1 - pow(m5,2))*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*n5 + 
                             k5*sqrt(1 - pow(m5,2))*(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                                k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*sqrt(1 - pow(n5,2)))*
                           (pow(k3*m3 + k4*m4 + k5*m5,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                               k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                         ((pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
                           (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/
                    126. + (((k3*m3*(k3*m3 + k4*m4 + k5*m5) + 
                           k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                           k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                            (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                              k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                         ((5*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                                k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                                 (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                            (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
                           ((k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                                k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                              (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                                pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                            ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                              (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2)))))
                         )/(pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                      ((k3*m3*(k4*m4 + k5*m5) + k3*sqrt(1 - pow(m3,2))*n3*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                           k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                            (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                         (pow(k3*m3 + k4*m4 + k5*m5,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + 
                             k5*sqrt(1 - pow(m5,2))*n5,2) + 
                           pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                             k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))*
                         ((3*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                                k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                                 (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                            (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
                           (2*(k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                                k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                              (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                                pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                            ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                              (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2)))))
                         )/(7.*(pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                         (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                           pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))))/18.))/
               (pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
              ((k2*m2*(k3*m3 + k4*m4 + k5*m5) + k2*sqrt(1 - pow(m2,2))*n2*
                    (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                   k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                    (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                      k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                 (pow(k2*m2 + k3*m3 + k4*m4 + k5*m5,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + 
                     k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                   pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                     k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))*
                 ((((3*(k3*m3*(k3*m3 + k4*m4) + k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                             k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                              (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                         (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                        (2*(k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*n3*n4 + 
                             k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*sqrt(1 - pow(n4,2)))*
                           (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                         ((pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                           (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))))*
                      ((3*((k3*m3 + k4*m4)*(k3*m3 + k4*m4 + k5*m5) + 
                             (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*
                              (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                             (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*
                              (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                                k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                         (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                           pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)) + 
                        (3*(k5*(k3*m3 + k4*m4)*m5 + k5*sqrt(1 - pow(m5,2))*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*n5 + 
                             k5*sqrt(1 - pow(m5,2))*(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                                k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*sqrt(1 - pow(n5,2)))*
                           (pow(k3*m3 + k4*m4 + k5*m5,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                               k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                         ((pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
                           (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/
                    126. + ((3*(k3*m3*(k3*m3 + k4*m4 + k5*m5) + 
                           k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                           k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                            (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                              k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                         ((5*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                                k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                                 (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                            (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
                           ((k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                                k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                              (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                                pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                            ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                              (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2)))))
                         )/(7.*(pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))))
                        + (3*(k3*m3*(k4*m4 + k5*m5) + k3*sqrt(1 - pow(m3,2))*n3*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                           k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                            (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                         (pow(k3*m3 + k4*m4 + k5*m5,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + 
                             k5*sqrt(1 - pow(m5,2))*n5,2) + 
                           pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                             k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))*
                         ((3*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                                k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                                 (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                            (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
                           (2*(k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                                k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                              (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                                pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                            ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                              (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2)))))
                         )/(7.*(pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                         (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                           pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))))/18.))/
               ((pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                 (pow(k3*m3 + k4*m4 + k5*m5,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                   pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                     k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))))/33.))/
       (pow(k1,2)*pow(m1,2) + pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2))) + 
      ((k1*m1*(k2*m2 + k3*m3 + k4*m4 + k5*m5) + k1*sqrt(1 - pow(m1,2))*n1*
            (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
           k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2))*(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
              k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
              k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
         (pow(k1*m1 + k2*m2 + k3*m3 + k4*m4 + k5*m5,2) + pow(k1*sqrt(1 - pow(m1,2))*n1 + k2*sqrt(1 - pow(m2,2))*n2 + 
             k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
           pow(k1*sqrt(1 - pow(m1,2))*sqrt(1 - pow(n1,2)) + k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + 
             k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
             k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))*
         ((((((3*(k2*m2*(k2*m2 + k3*m3) + k2*sqrt(1 - pow(m2,2))*n2*(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3) + 
                          k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                           (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))))/
                      (pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
                     (2*(k2*k3*m2*m3 + k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*n2*n3 + 
                          k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n2,2))*sqrt(1 - pow(n3,2)))*
                        (pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)))/
                      ((pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                        (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))))*
                   ((3*((k2*m2 + k3*m3)*(k2*m2 + k3*m3 + k4*m4) + 
                          (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3)*
                           (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                          (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
                           (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                             k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                      (pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                        pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)) + 
                     (3*(k4*(k2*m2 + k3*m3)*m4 + k4*sqrt(1 - pow(m4,2))*(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3)*n4 + 
                          k4*sqrt(1 - pow(m4,2))*(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
                           sqrt(1 - pow(n4,2)))*(pow(k2*m2 + k3*m3 + k4*m4,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                            k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                      ((pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                          pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))*
                        (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/126.
                  + ((3*(k2*m2*(k2*m2 + k3*m3 + k4*m4) + k2*sqrt(1 - pow(m2,2))*n2*
                         (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                        k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                         (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                           k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
                      ((5*(k3*m3*(k3*m3 + k4*m4) + k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                             k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                              (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                         (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                        ((k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*n3*n4 + 
                             k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*sqrt(1 - pow(n4,2)))*
                           (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                         ((pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                           (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/
                    (7.*(pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))) + 
                   (3*(k2*m2*(k3*m3 + k4*m4) + k2*sqrt(1 - pow(m2,2))*n2*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                        k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                         (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))))*
                      (pow(k2*m2 + k3*m3 + k4*m4,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,
                         2) + pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                          k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
                      ((3*(k3*m3*(k3*m3 + k4*m4) + k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                             k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                              (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                         (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                        (2*(k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*n3*n4 + 
                             k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*sqrt(1 - pow(n4,2)))*
                           (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                         ((pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                           (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))))))/
                    (7.*(pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                      (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                        pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))))/18.)*
              ((3*((k2*m2 + k3*m3 + k4*m4)*(k2*m2 + k3*m3 + k4*m4 + k5*m5) + 
                     (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*
                      (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                     (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                        k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*
                      (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                        k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                 (pow(k2*m2 + k3*m3 + k4*m4,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                   pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                     k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)) + 
                (4*(k5*(k2*m2 + k3*m3 + k4*m4)*m5 + k5*sqrt(1 - pow(m5,2))*
                      (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*n5 + 
                     k5*sqrt(1 - pow(m5,2))*(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                        k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*sqrt(1 - pow(n5,2)))*
                   (pow(k2*m2 + k3*m3 + k4*m4 + k5*m5,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                       k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                 ((pow(k2*m2 + k3*m3 + k4*m4,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                       k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
                   (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/33. + 
           (((3*(k2*m2*(k2*m2 + k3*m3) + k2*sqrt(1 - pow(m2,2))*n2*(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3) + 
                     k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                      (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))))/
                 (pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
                (2*(k2*k3*m2*m3 + k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*n2*n3 + 
                     k2*k3*sqrt(1 - pow(m2,2))*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n2,2))*sqrt(1 - pow(n3,2)))*
                   (pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2)))/
                 ((pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                   (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))))*
              ((3*((k2*m2 + k3*m3)*(k2*m2 + k3*m3 + k4*m4 + k5*m5) + 
                     (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3)*
                      (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                     (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
                      (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                        k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                   ((5*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                          k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                           (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                      (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
                     ((k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                          k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                        (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                          pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                      ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                        (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/
                 (7.*(pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))) + 
                (4*((k2*m2 + k3*m3)*(k4*m4 + k5*m5) + (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3)*
                      (k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                     (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)))*
                      (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                   (pow(k2*m2 + k3*m3 + k4*m4 + k5*m5,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                       k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))*
                   ((3*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                          k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                           (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                      (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
                     (2*(k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                          k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                        (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                          pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                      ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                        (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/
                 (7.*(pow(k2*m2 + k3*m3,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3,2) + 
                     pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)),2))*
                   (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                     pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))))/231. + 
           ((3*(k2*m2*(k2*m2 + k3*m3 + k4*m4 + k5*m5) + k2*sqrt(1 - pow(m2,2))*n2*
                    (k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                   k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                    (k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                      k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                 ((((3*(k3*m3*(k3*m3 + k4*m4) + k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                             k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                              (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                         (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                        (2*(k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*n3*n4 + 
                             k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*sqrt(1 - pow(n4,2)))*
                           (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                         ((pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                           (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))))*
                      ((7*((k3*m3 + k4*m4)*(k3*m3 + k4*m4 + k5*m5) + 
                             (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*
                              (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                             (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*
                              (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                                k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                         (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                           pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)) + 
                        ((k5*(k3*m3 + k4*m4)*m5 + k5*sqrt(1 - pow(m5,2))*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*n5 + 
                             k5*sqrt(1 - pow(m5,2))*(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                                k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*sqrt(1 - pow(n5,2)))*
                           (pow(k3*m3 + k4*m4 + k5*m5,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                               k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                         ((pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
                           (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/
                    126. + (((k3*m3*(k3*m3 + k4*m4 + k5*m5) + 
                           k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                           k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                            (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                              k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                         ((5*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                                k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                                 (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                            (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
                           ((k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                                k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                              (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                                pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                            ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                              (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2)))))
                         )/(pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                      ((k3*m3*(k4*m4 + k5*m5) + k3*sqrt(1 - pow(m3,2))*n3*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                           k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                            (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                         (pow(k3*m3 + k4*m4 + k5*m5,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + 
                             k5*sqrt(1 - pow(m5,2))*n5,2) + 
                           pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                             k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))*
                         ((3*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                                k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                                 (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                            (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
                           (2*(k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                                k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                              (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                                pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                            ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                              (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2)))))
                         )/(7.*(pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                         (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                           pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))))/18.))/
               (pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2))) + 
              (4*(k2*m2*(k3*m3 + k4*m4 + k5*m5) + k2*sqrt(1 - pow(m2,2))*n2*
                    (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                   k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2))*
                    (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                      k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                 (pow(k2*m2 + k3*m3 + k4*m4 + k5*m5,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + 
                     k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                   pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                     k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))*
                 ((((3*(k3*m3*(k3*m3 + k4*m4) + k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4) + 
                             k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                              (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))))/
                         (pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))) + 
                        (2*(k3*k4*m3*m4 + k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*n3*n4 + 
                             k3*k4*sqrt(1 - pow(m3,2))*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n3,2))*sqrt(1 - pow(n4,2)))*
                           (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)))/
                         ((pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                           (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))))*
                      ((3*((k3*m3 + k4*m4)*(k3*m3 + k4*m4 + k5*m5) + 
                             (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*
                              (k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                             (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*
                              (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                                k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                         (pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                           pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2)) + 
                        (3*(k5*(k3*m3 + k4*m4)*m5 + k5*sqrt(1 - pow(m5,2))*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4)*n5 + 
                             k5*sqrt(1 - pow(m5,2))*(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
                                k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)))*sqrt(1 - pow(n5,2)))*
                           (pow(k3*m3 + k4*m4 + k5*m5,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                               k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                         ((pow(k3*m3 + k4*m4,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4,2) + 
                             pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)),2))*
                           (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2))))))/
                    126. + ((3*(k3*m3*(k3*m3 + k4*m4 + k5*m5) + 
                           k3*sqrt(1 - pow(m3,2))*n3*(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                           k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                            (k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                              k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                         ((5*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                                k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                                 (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                            (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
                           ((k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                                k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                              (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                                pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                            ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                              (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2)))))
                         )/(7.*(pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2))))
                        + (3*(k3*m3*(k4*m4 + k5*m5) + k3*sqrt(1 - pow(m3,2))*n3*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                           k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2))*
                            (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2))))*
                         (pow(k3*m3 + k4*m4 + k5*m5,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + 
                             k5*sqrt(1 - pow(m5,2))*n5,2) + 
                           pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                             k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))*
                         ((3*(k4*m4*(k4*m4 + k5*m5) + k4*sqrt(1 - pow(m4,2))*n4*(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5) + 
                                k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2))*
                                 (k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)))))/
                            (pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2))) + 
                           (2*(k4*k5*m4*m5 + k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*n4*n5 + 
                                k4*k5*sqrt(1 - pow(m4,2))*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n4,2))*sqrt(1 - pow(n5,2)))*
                              (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                                pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2)))/
                            ((pow(k4,2)*pow(m4,2) + pow(k4,2)*(1 - pow(m4,2))*pow(n4,2) + pow(k4,2)*(1 - pow(m4,2))*(1 - pow(n4,2)))*
                              (pow(k5,2)*pow(m5,2) + pow(k5,2)*(1 - pow(m5,2))*pow(n5,2) + pow(k5,2)*(1 - pow(m5,2))*(1 - pow(n5,2)))))
                         )/(7.*(pow(k3,2)*pow(m3,2) + pow(k3,2)*(1 - pow(m3,2))*pow(n3,2) + pow(k3,2)*(1 - pow(m3,2))*(1 - pow(n3,2)))*
                         (pow(k4*m4 + k5*m5,2) + pow(k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                           pow(k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))))/18.))/
               ((pow(k2,2)*pow(m2,2) + pow(k2,2)*(1 - pow(m2,2))*pow(n2,2) + pow(k2,2)*(1 - pow(m2,2))*(1 - pow(n2,2)))*
                 (pow(k3*m3 + k4*m4 + k5*m5,2) + pow(k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + k5*sqrt(1 - pow(m5,2))*n5,2) + 
                   pow(k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + 
                     k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))))/33.))/
       ((pow(k1,2)*pow(m1,2) + pow(k1,2)*(1 - pow(m1,2))*pow(n1,2) + pow(k1,2)*(1 - pow(m1,2))*(1 - pow(n1,2)))*
         (pow(k2*m2 + k3*m3 + k4*m4 + k5*m5,2) + pow(k2*sqrt(1 - pow(m2,2))*n2 + k3*sqrt(1 - pow(m3,2))*n3 + k4*sqrt(1 - pow(m4,2))*n4 + 
             k5*sqrt(1 - pow(m5,2))*n5,2) + pow(k2*sqrt(1 - pow(m2,2))*sqrt(1 - pow(n2,2)) + k3*sqrt(1 - pow(m3,2))*sqrt(1 - pow(n3,2)) + 
             k4*sqrt(1 - pow(m4,2))*sqrt(1 - pow(n4,2)) + k5*sqrt(1 - pow(m5,2))*sqrt(1 - pow(n5,2)),2))))/52. ;

return (temp2);
}
