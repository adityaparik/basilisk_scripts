#include "utils.h"
#include "fractions.h"

double interfacearea_axi (scalar f) {
  double sb = 0.;
  double area = 0.;

  foreach( reduction(+:sb) reduction(+:area) ) {
    sb += f[]*2*pi*y*sq(Delta);
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = interface_normal (point, f), p;
      double alpha = plane_alpha (f[], n);
      area += 2*pi*y*pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
    }    
  }
  return area;
}

vector vsx[];        // Viscous Stress Tensor: Row 1
vector vsy[];        // Viscous Stress Tensor: Row 2
vector vsz[];        // Viscous Stress Tensor: Row 3


void viscous_tensor (vector u) {
  foreach() {
      vsx.x[] = mu2( (u.x[1,0,0] - u.x[-1,0,0])/(2.*Delta) + (u.x[1,0,0] - u.x[-1,0,0])/(2.*Delta) );
      vsx.y[] = mu2( (u.x[0,1,0] - u.x[0,-1,0])/(2.*Delta) + (u.y[1,0,0] - u.y[-1,0,0])/(2.*Delta) );
      vsx.z[] = mu2( (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta) + (u.z[1,0,0] - u.z[-1,0,0])/(2.*Delta) );
  
      vsy.x[] = mu2( (u.y[1,0,0] - u.y[-1,0,0])/(2.*Delta) + (u.x[0,1,0] - u.x[0,-1,0])/(2.*Delta) );
      vsy.y[] = mu2( (u.y[0,1,0] - u.y[0,-1,0])/(2.*Delta) + (u.y[0,1,0] - u.y[0,-1,0])/(2.*Delta) );
      vsy.z[] = mu2( (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta) + (u.z[0,1,0] - u.z[0,-1,0])/(2.*Delta) );
  
      vsz.x[] = mu2( (u.z[1,0,0] - u.z[-1,0,0])/(2.*Delta) + (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta) );
      vsz.y[] = mu2( (u.z[0,1,0] - u.z[0,-1,0])/(2.*Delta) + (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta) );
      vsz.z[] = mu2( (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta) + (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta) );   
  }
}
  
// Total Stress Tensor Calculation : ts = -p + vs

vector tsx[];        // Total Stress Tensor: Row 1
vector tsy[];        // Total Stress Tensor: Row 2
vector tsz[];        // Total Stress Tensor: Row 3

void total_stress (vector u) {
  viscous_tensor(u);
  foreach() {
      tsx.x[] = -p[] + vsx.x[];
      tsx.y[] = vsx.y[];
      tsx.z[] = vsx.z[];
  
      tsy.x[] = vsy.x[];
      tsy.y[] = -p[] + vsy.y[];
      tsy.z[] = vsy.z[];
  
      tsz.x[] = vsz.x[];
      tsz.y[] = vsz.y[];
      tsz.z[] = -p[] + vsz.z[];
  }
}


// Function to obtain surface integral of surface stresses
vector Tr[];        // Traction Vector to store stress vector at each interface point

double Fx=0;        // x-component of total force on an object
double Fy=0;        // y-component of total force on an object
double Fz=0;        // z-component of total force on an object
double Fmag=0;      // magnitude of total force on an object

void surfaceforces_axi (scalar f, vector u) {
  total_stress(u);        // Viscous and total stress tensors are created

  foreach(reduction(+:Fx) reduction(+:Fy) reduction(+:Fz)) {
      if (f[] > 1e-4 && f[] < 1. - 1e-4) {
          // VOF interface segment outward normal, alpha calculation
          coord m = mycs (point, f);
          double alpha = plane_alpha (f[], m);
          coord pp;       // Variable to store interface segment centroid
          
          Tr.x[] = ( tsx.x[]*m.x + tsx.y[]*m.y + tsx.z[]*m.z );
          Tr.y[] = ( tsy.x[]*m.x + tsy.y[]*m.y + tsy.z[]*m.z );
          Tr.z[] = ( tsz.x[]*m.x + tsz.y[]*m.y + tsz.z[]*m.z );
          
          // Area of VOF interface segment in the current cell
          double segmentarea = 2*pi*y*Delta*plane_area_center(m, alpha, &pp);
  
          Fx += Tr.x[] * segmentarea;
          Fy += Tr.y[] * segmentarea;
          Fz += Tr.z[] * segmentarea;
      }
  }
  Fmag = sqrt(sq(Fx) + sq(Fy) + sq(Fz));
}
  
