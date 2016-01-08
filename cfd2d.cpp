#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include "riemann.h"

using namespace std;

#define gamma 1.4

double **Q1, **Q2, **Q3, **Q4, **dF1, **dF2, **dF3, **dF4, **dG1, **dG2, **dG3, **dG4;
double dx, dy, dt;
double **rho, **u, **v, **e, **p, **c;
int nmax;
double imax, jmax, cfl;


//*******************************

double **Array2D(unsigned int row, unsigned int col)
{
double **x;
x = new double *[row]; 
for (int count = 0; count < row; count++)
{
 x[count] = new double[col];
}
return x;
}

//*******************************

void cons_2_prim()
{
 for(int i = 0; i < imax; i++)
 { 
  for(int j = 0; j < jmax; j++)
  {
   rho[i][j] = Q1[i][j];
   u[i][j] = Q2[i][j]/Q1[i][j];
   v[i][j] = Q3[i][j]/Q1[i][j];
   double ke = u[i][j]*u[i][j]/2.0 + v[i][j]*v[i][j]/2.0;
   e[i][j] = Q4[i][j]/Q1[i][j] - ke;
   p[i][j] = (gamma-1.0)*rho[i][j]*e[i][j];
   c[i][j] = sqrt(gamma*p[i][j]/rho[i][j]);
  }
 }
}

//*******************************

void Initialize(void)
{

  // read input file
  string s;  
  ifstream f;
  f.open ("input_file", ios_base::in);
  f>>s>>nmax;
  f>>s>>dx;
  f>>s>>dy;
  f>>s>>imax;
  f>>s>>jmax;
  f>>s>>cfl;
  f.close();
  cout<<"nmax = "<<nmax<<endl;
  cout<<"dx = "<<dx<<endl;
  cout<<"dy = "<<dy<<endl;
  cout<<"imax = "<<imax<<endl;
  cout<<"jmax = "<<jmax<<endl;
  cout<<"CFL = "<<cfl<<endl;
  

  // conservative
  Q1 = Array2D(imax,jmax);
  Q2 = Array2D(imax,jmax);
  Q3 = Array2D(imax,jmax);
  Q4 = Array2D(imax,jmax);
  
  // delta flux
  dF1 = Array2D(imax,jmax);
  dF2 = Array2D(imax,jmax);
  dF3 = Array2D(imax,jmax);
  dF4 = Array2D(imax,jmax);
  dG1 = Array2D(imax,jmax);
  dG2 = Array2D(imax,jmax);
  dG3 = Array2D(imax,jmax);
  dG4 = Array2D(imax,jmax);
  
  // primitive
  rho = Array2D(imax,jmax);
  u = Array2D(imax,jmax);
  v = Array2D(imax,jmax);
  e = Array2D(imax,jmax);
  p = Array2D(imax,jmax);
  c = Array2D(imax,jmax);

  //driver 
  double rho2 = 2.0;
  double u2 = 0.0;
  double v2 = 0.0;
  double p2 = 2.0;
  double e2 = p2/(gamma-1.0)/rho2;
  double qq2[4];
  qq2[0] = rho2;
  qq2[1] = rho2*u2;
  qq2[2] = rho2*v2;
  qq2[3] = rho2*(e2+u2*u2/2.0+v2*v2/2.0);

  //driven
  double rho1 = 1.0;
  double u1 = 0.0;
  double v1 = 0.0;
  double p1 = 1.0;
  double e1 = p1/(gamma-1.0)/rho1;
  double qq1[4];
  qq1[0] = rho1;
  qq1[1] = rho1*u1;
  qq1[2] = rho1*v1;
  qq1[3] = rho1*(e1+u1*u1/2.0+v1*v1/2.0);

  double x, y, r;

  for(int i=0; i<imax; i++)
  {
   for(int j=0; j<jmax; j++)
   {
    x = (double(i)+0.5)*dx;   
    y = (double(j)+0.5)*dy;   
    r = sqrt(x*x + y*y);

    if(r < 0.2)
    {
     Q1[i][j] = qq2[0];
     Q2[i][j] = qq2[1];
     Q3[i][j] = qq2[2];
     Q4[i][j] = qq2[3];
    }
    else
    {
     Q1[i][j] = qq1[0];
     Q2[i][j] = qq1[1];
     Q3[i][j] = qq1[2];
     Q4[i][j] = qq1[3];
    } 
    dF1[i][j] = 0.0;
    dF2[i][j] = 0.0;
    dF3[i][j] = 0.0;
    dF4[i][j] = 0.0;
    dG1[i][j] = 0.0;
    dG2[i][j] = 0.0;
    dG3[i][j] = 0.0;
    dG4[i][j] = 0.0;
   }
  }

}


//*******************************

void calc_dt()
{
 double vel = 0.0;
 double upc;
 for(int i = 0; i < imax; i++)
 {
  for(int j = 0; j < jmax; j++)
  { 
   upc = sqrt(u[i][j]*u[i][j] + v[i][j]*v[i][j]) + c[i][j];
   vel = (vel>upc)? vel : upc;
  }
 }
 dt = cfl*sqrt((dx*dx+dy*dy)/2.0)/vel; 
}

//*******************************

void calc_min_max()
{
 double min_rho = 1.0e5;
 double max_rho = 0.0;
 double min_p = 1.0e5;
 double max_p = 0.0;
 
 for(int i = 0; i < imax; i++)
 {
  for(int j = 0; j < jmax; j++)
  {
   min_rho = (min_rho<rho[i][j])? min_rho : rho[i][j];
   min_p = (min_p<p[i][j])? min_p : p[i][j];
   max_rho = (max_rho>rho[i][j])? max_rho : rho[i][j];
   max_p = (max_p>p[i][j])? max_p : p[i][j];
  }
 }
 cout<<"min/max rho = "<<min_rho<<" "<<max_rho<<endl;
 cout<<"min/max p = "<<min_p<<" "<<max_p<<endl;
}


//*******************************

void calc_flux()
{
 double left[6], right[6];
 double fimh[4], fiph[4]; 
 double gjmh[4], gjph[4]; 

 for(int i = 0; i < imax; i++)
 {
  for(int j = 0; j < jmax; j++)
  {


  // i-1/2 
  if(i>0)
  {
   left[0] = rho[i-1][j];
   left[1] = u[i-1][j];
   left[2] = v[i-1][j];
   left[3] = e[i-1][j];
   left[4] = p[i-1][j];
   left[5] = c[i-1][j];
  }
  else
  {
   // ghost cell : wall bc
   left[0] = rho[i][j];
   left[1] = -u[i][j];
   left[2] = v[i][j];
   left[3] = e[i][j];
   left[4] = p[i][j];
   left[5] = c[i][j];
  }
  right[0] = rho[i][j];
  right[1] = u[i][j];
  right[2] = v[i][j];
  right[3] = e[i][j];
  right[4] = p[i][j];
  right[5] = c[i][j];    
  riemann_solver(gamma,1,left,right,fimh); 

  // i+1/2
  left[0] = rho[i][j];
  left[1] = u[i][j];
  left[2] = v[i][j];
  left[3] = e[i][j];
  left[4] = p[i][j];
  left[5] = c[i][j];
  if(i<imax-1) 
  {  
   right[0] = rho[i+1][j];
   right[1] = u[i+1][j];
   right[2] = v[i+1][j];
   right[3] = e[i+1][j];
   right[4] = p[i+1][j];
   right[5] = c[i+1][j]; 
  }
  else
  {
   // ghost cell : outflow bc
   right[0] = rho[i][j];
   right[1] = u[i][j];
   right[2] = v[i][j];
   right[3] = e[i][j];
   right[4] = p[i][j];
   right[5] = c[i][j];
  }
  riemann_solver(gamma,1,left,right,fiph);

  // delta flux 
  dF1[i][j] = (fiph[0] - fimh[0])/dx;
  dF2[i][j] = (fiph[1] - fimh[1])/dx;
  dF3[i][j] = (fiph[2] - fimh[2])/dx;
  dF4[i][j] = (fiph[3] - fimh[3])/dx;


//---------------


  // j-1/2 
  if(j>0)
  {
   left[0] = rho[i][j-1];
   left[1] = u[i][j-1];
   left[2] = v[i][j-1];
   left[3] = e[i][j-1];
   left[4] = p[i][j-1];
   left[5] = c[i][j-1];
  }
  else
  {
   // ghost cell : wall bc
   left[0] = rho[i][j];
   left[1] = u[i][j];
   left[2] = v[i][j];
   left[3] = e[i][j];
   left[4] = p[i][j];
   left[5] = c[i][j];
  }
  right[0] = rho[i][j];
  right[1] = u[i][j];
  right[2] = v[i][j];
  right[3] = e[i][j];
  right[4] = p[i][j];
  right[5] = c[i][j];    
  riemann_solver(gamma,2,left,right,gjmh); 

  // j+1/2
  left[0] = rho[i][j];
  left[1] = u[i][j];
  left[2] = v[i][j];
  left[3] = e[i][j];
  left[4] = p[i][j];
  left[5] = c[i][j];
  if(j<jmax-1) 
  {  
   right[0] = rho[i][j+1];
   right[1] = u[i][j+1];
   right[2] = v[i][j+1];
   right[3] = e[i][j+1];
   right[4] = p[i][j+1];
   right[5] = c[i][j+1]; 
  }
  else
  {
   // ghost cell : outflow bc
   right[0] = rho[i][j];
   right[1] = u[i][j];
   right[2] = -v[i][j];
   right[3] = e[i][j];
   right[4] = p[i][j];
   right[5] = c[i][j];
  }
  riemann_solver(gamma,2,left,right,gjph);

  // delta flux 
  dG1[i][j] = (gjph[0] - gjmh[0])/dy;
  dG2[i][j] = (gjph[1] - gjmh[1])/dy;
  dG3[i][j] = (gjph[2] - gjmh[2])/dy;
  dG4[i][j] = (gjph[3] - gjmh[3])/dy;

  }
 }
 
}


//*******************************

void Update(void)
{
 for(int i = 0; i < imax; i++)
 {
  for(int j = 0; j < jmax; j++)
  {
   Q1[i][j] -= dt*(dF1[i][j]+dG1[i][j]);
   Q2[i][j] -= dt*(dF2[i][j]+dG2[i][j]);
   Q3[i][j] -= dt*(dF3[i][j]+dG3[i][j]);
   Q4[i][j] -= dt*(dF4[i][j]+dG4[i][j]);
  }
 }
}

//*******************************

void Evolve(void)
{
 for(int n = 0; n < nmax; n++)
 {
  cons_2_prim();
  calc_dt();
  cout<<"*****************"<<endl;
  cout<<"n = "<<n<<"; dt = "<<dt<<" \n";
  calc_min_max();
  calc_flux();
  Update();
 } 
  cout<<"*****************"<<endl;
}

//*******************************

void Output_results()
{
 double x, y, vel;
 ofstream f, fx, fy;

 // (x,y)
 fx.open ("x_output.csv", ios_base::out);
 fy.open ("y_output.csv", ios_base::out);
 for (int i = 0; i < imax; i++)
 {
  for (int j = 0; j < jmax; j++)
  {
   x = (double(i)+0.5)*dx;
   y = (double(j)+0.5)*dy;
   if(j<jmax-1)
   {
    fx<<x<<",";
    fy<<y<<",";
   }
   else
   {
    fx<<x;
    fy<<y;
   }
  }
  fx<<endl;
  fy<<endl;
 }
 fx.close();
 fy.close();

 // rho 
 f.open ("rho_output.csv", ios_base::out);
 for (int i = 0; i < imax; i++)
 {
  for (int j = 0; j < jmax; j++)
  {
   if(j<jmax-1)
   {
    f<<rho[i][j]<<",";
   }
   else
   {
    f<<rho[i][j];
   }
  }
  f<<endl;
 }
 f.close(); 

 // p 
 f.open ("p_output.csv", ios_base::out);
 for (int i = 0; i < imax; i++)
 {
  for (int j = 0; j < jmax; j++)
  {
   if(j<jmax-1)
   {
    f<<p[i][j]<<",";
   }
   else
   {
    f<<rho[i][j];
   }
  }
  f<<endl;
 }
 f.close(); 

 // vel 
 f.open ("vel_output.csv", ios_base::out);
 for (int i = 0; i < imax; i++)
 {
  for (int j = 0; j < jmax; j++)
  {
   vel = sqrt(pow(u[i][j],2.0) + pow(v[i][j],2.0));
   if(j<jmax-1)
   {
    f<<vel<<",";
   }
   else
   {
    f<<vel;
   }
  }
  f<<endl;
 }
 f.close(); 

}

//*******************************

int main()
{
 Initialize();
 Evolve();
 Output_results();
 return 0;
}

//*******************************

