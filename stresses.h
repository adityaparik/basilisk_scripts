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
      area += 2*pi*y*pow(Delta, 2 - 1)*plane_area_center (n, alpha, &p);
    }    
  }
  return area;
}

// Function to calculate Rate of Deformation Tensor
void deformation_tensor (const vector u, tensor D) {
  foreach() {
      D.x.x[] = (u.x[1,0,0] - u.x[-1,0,0])/(2.*Delta) + (u.x[1,0,0] - u.x[-1,0,0])/(2.*Delta);
      D.x.y[] = (u.x[0,1,0] - u.x[0,-1,0])/(2.*Delta) + (u.y[1,0,0] - u.y[-1,0,0])/(2.*Delta);
#if dimension==3
      D.x.z[] = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta) + (u.z[1,0,0] - u.z[-1,0,0])/(2.*Delta);
#endif
  
      D.y.x[] = (u.y[1,0,0] - u.y[-1,0,0])/(2.*Delta) + (u.x[0,1,0] - u.x[0,-1,0])/(2.*Delta);
      D.y.y[] = (u.y[0,1,0] - u.y[0,-1,0])/(2.*Delta) + (u.y[0,1,0] - u.y[0,-1,0])/(2.*Delta);
#if dimension==3
      D.y.z[] = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta) + (u.z[0,1,0] - u.z[0,-1,0])/(2.*Delta);
#endif
  
      D.z.x[] = (u.z[1,0,0] - u.z[-1,0,0])/(2.*Delta) + (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
      D.z.y[] = (u.z[0,1,0] - u.z[0,-1,0])/(2.*Delta) + (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
#if dimension==3
      D.z.z[] = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta) + (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);   
#endif
  }
}
  
// Newtonian Stress Tensor Calculation : ns = -p + mu*D 
void newtonian_stress (const scalar p, const tensor D, tensor ns) {
    deformation_tensor (u);
    foreach() {
        ns.x.x[] = -p[] + mu2*D.x.x[];
        ns.x.y[] = mu2*D.x.y[];
#if dimension==3
        ns.x.z[] = mu2*D.x.z[];
#endif
    
        ns.y.x[] = mu2*D.y.x[];
        ns.y.y[] = -p[] + mu2*D.y.y[];
#if dimension==3
        ns.y.z[] = mu2*D.y.z[];
#endif
   
        ns.z.x[] = mu2*D.z.x[];
        ns.z.y[] = mu2*D.z.y[];
#if dimension==3
        ns.z.z[] = -p[] + mu2*D.z.z[];
#endif
    }
}


// Function to obtain surface integral of surface stresses

coord surfaceforces (const scalar f, const tensor ts, vector Tr) {
    foreach() {
#if dimension!=3
        coord F_drop = {0.,0.};
#else
        coord F_drop = {0.,0.,0.};
#endif
        if (f[] > 1e-6 && f[] < 1. - 1e-6) {
            // VOF interface segment outward normal, alpha calculation
            coord m = mycs (point, f);
            double alpha = plane_alpha (f[], m);
            coord pp;       // Variable to store interface segment centroid
            
            Tr.x[] = ( tsx.x[]*m.x + tsx.y[]*m.y + tsx.z[]*m.z );
            Tr.y[] = ( tsy.x[]*m.x + tsy.y[]*m.y + tsy.z[]*m.z );
#if dimension==3
            Tr.z[] = ( tsz.x[]*m.x + tsz.y[]*m.y + tsz.z[]*m.z );
            // Area of VOF interface segment for 3D domain
            double segmentarea = Delta*Delta*plane_area_center(m, alpha, &pp);
#elif AXI
            // Area of VOF interface segment for Axisymmetric domain
            double segmentarea = 2*pi*y*Delta*plane_area_center(m, alpha, &pp);
#endif
            F_drop.x += Tr.x[] * segmentarea;
            F_drop.y += Tr.y[] * segmentarea;
#if dimension==3
            F_drop.z += Tr.z[] * segmentarea;
#endif

            return F_drop;
        }
    }
}
