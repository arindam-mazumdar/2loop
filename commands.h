void write2file(char x[], float y, float z)
{

char file[32];
FILE *fptr;
sprintf(file, "./output/%s.dat", x);
fptr = fopen(file,"a");
fprintf(fptr,"%f %f\n", y, z);
fclose(fptr);


}

/*void readfile2array(char x[], int maxrow, float* arrayname[maxrow][2])//arguments 1. FILENAME, 2. variable1, 3. variable2, 4.array-name
{
char file[32];
int numcol,numrow;
FILE *fptr;
sprintf(file, "./input/%s.dat", x);
fptr = fopen(file,"r");
for(numrow=1;numrow<maxrow;numrow++)
{for(numcol=1;numcol<3;numcol++)
{fscanf(fptr, "%f",&arrayname[numrow][numcol]);
}}
fclose(fptr);


}*/

float F2sym(float k, float q, float mq);
float F3sym1(float k, float q, float mq);
float F3sym2(float k, float q, float p, float mp, float mq, float np);
float F4(float k1, float k2, float k3, float k4, float m1, float m2, float m3, float m4, float n1, float n2, float n3, float n4);
float F5(float k1, float k2, float k3, float k4, float k5, float m1, float m2, float m3, float m4, float m5, float n1, float n2, float n3, float n4, float n5);
//float F511(float k, float q, float p, float m2, float m1, float n1);
//float F512(float k, float q, float p, float m2, float m1, float n1);
//float F513(float k, float q, float p, float m2, float m1, float n1);
//float F521(float k, float q, float p, float m2, float m1, float n1);
//float F531(float k, float q, float p, float m2, float m1, float n1);
