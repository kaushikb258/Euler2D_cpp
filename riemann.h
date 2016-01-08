#include <stdlib.h>

void riemann_solver(double gam, int dir, double left[], double right[], double flux[])
{
 double rhol = left[0];
 double ul = left[1];
 double vl = left[2];
 double el = left[3];
 double pl = left[4];
 double cl = left[5];

 double rhor = right[0];
 double ur = right[1];
 double vr = right[2];
 double er = right[3];
 double pr = right[4];
 double cr = right[5];

 double ql, qr, nx, ny;
 if(dir==1)
 {
  ql = ul;
  qr = ur;
  nx = 1.0;
  ny = 0.0;
 }
 else if(dir==2)
 {
  ql = vl;
  qr = vr;
  nx = 0.0;
  ny = 1.0;
 }

 double cwave, wl, wr;
 wl = abs(ql)+cl;
 wr = abs(qr)+cr;
 cwave = (wl>wr)? wl : wr;

 double fl[4], fr[4];
 
 // left flux
 fl[0] = rhol*ql;
 fl[1] = rhol*ql*ul + pl*nx;
 fl[2] = rhol*ql*vl + pl*ny;
 fl[3] = rhol*ql*(el+ul*ul/2.0+vl*vl/2.0) + pl*ql;

 // right flux
 fr[0] = rhor*qr;
 fr[1] = rhor*qr*ur + pr*nx;
 fr[2] = rhor*qr*vr + pr*ny;
 fr[3] = rhor*qr*(er+ur*ur/2.0+vr*vr/2.0) + pr*qr;

 double Ql[4], Qr[4];

 Ql[0] = rhol;
 Ql[1] = rhol*ul;
 Ql[2] = rhol*vl;
 Ql[3] = rhol*(el+ul*ul/2.0+vl*vl/2.0);
 
 Qr[0] = rhor;
 Qr[1] = rhor*ur;
 Qr[2] = rhor*vr;
 Qr[3] = rhor*(er+ur*ur/2.0+vr*vr/2.0); 

 for(int i = 0; i<4; i++)
 {
  flux[i] = 0.5*(fl[i]+fr[i]) - 0.5*cwave*(Qr[i]-Ql[i]);
 }

}

//*******************************
