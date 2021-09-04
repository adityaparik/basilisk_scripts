// This header file stores many custom functions
// If using Stress Tensor calculation functions, you must define if 
// you wish to use a specific type of Constitutive relationship.

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

int n_drops=-1;

typedef struct {
  double vol;
  double ke;        // Droplet absolute KE
//  double rel_ke=0.;
  double s_en;           // Droplet Surface Energy
#if dimension==3
  double cmz;
  double uz;
#endif
  double cmx;
  double cmy;
  double ux;
  double uy;
} DropStats;

DropStats * atoms ( const scalar f, const vector u, scalar m ) {
  foreach()
    m[] = f[] > 1e-3;
  n_drops = tag(m);
  
  if (n_drops==-1)
    static DropStats atom[n_drops];
  
// Initialize the stat variables to zero  
  for (int j = 0; j < n_drops; j++) {
    atom[j].vol;
    atom[j].ke=0.;        // Droplet absolute KE
  //  double rel_ke=0.;
    atom[j].se=0.;           // Droplet Surface Energy
    atom[j].cmx=0.;
    atom[j].cmy=0.;
    atom[j].ux=0.;
    atom[j].uy=0.;
#if dimension==3
    atom[j].cmz=0.;
    atom[j].uz=0.;
#endif
  }

  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      atom[j].vol += dv()*f[];
      coord p = {x,y,z};
      atom[j].cmx += dv()*f[]*p.x;
      atom[j].cmy += dv()*f[]*p.y;
      atom[j].ux += u.x[]*dv()*f[];
      atom[j].uy += u.y[]*dv()*f[];
#if dimension==3
      atom[j].cmz += dv()*f[]*p.z;
      atom[j].uz += u.z[]*dv()*f[];
#endif
      }
    }

  MPI_Allreduce (MPI_IN_PLACE, atom, n_drops, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  

  for (int j = 0; j < n_drops; j++) {
    atom[j].cm.x = atom[j].cm.x/atom[j].vol;
    atom[j].cm.y = atom[j].cm.y/atom[j].vol;
    atom[j].u.x = atom[j].cm.x/atom[j].vol;
    atom[j].u.y = atom[j].cm.y/atom[j].vol;
#if dimension==3
    atom[j].cm.z = atom[j].cm.z/atom[j].vol;
    atom[j].u.z = atom[j].cm.z/atom[j].vol;
#endif
  }
  foreach( reduction(min:xmin) reduction(max:xmax) reduction(min:ymin) reduction(max:ymax) ) {
    if (m[] > 0) {
      int j = m[] - 1;
      if (x>atom[j].xmax)  atom[j].xmax = x;
      if (x<atom[j].xmin)  atom[j].xmin = x;
      if (y>atom[j].ymax)  atom[j].ymax = y;
      if (y<atom[j].ymin)  atom[j].ymin = y;
#if dimension==3
      if (z>atom[j].zmax)  atom[j].zmax = z;
      if (z<atom[j].zmin)  atom[j].zmin = z;
#endif
    }
  }
}

    
        







// Function to calculate Rate of Deformation Tensor
void stress_tensors (const vector u, const scalar p, tensor D, tensor ts) {
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
  
      D.z.x[] = (u.z[1,0,0] - u.z[-1,0,0])/(2.*Delta) + (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
      D.z.y[] = (u.z[0,1,0] - u.z[0,-1,0])/(2.*Delta) + (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
      D.z.z[] = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta) + (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);   
#endif
  }

#ifdef NEWTONIAN
// Newtonian Stress Tensor Calculation : ns = -p + mu*D 
  foreach() {
      ts.x.x[] = -p[] + mu2*D.x.x[];
      ts.x.y[] = mu2*D.x.y[];
#if dimension==3
      ts.x.z[] = mu2*D.x.z[];
#endif
      ts.y.x[] = mu2*D.y.x[];
      ts.y.y[] = -p[] + mu2*D.y.y[];
#if dimension==3
      ts.y.z[] = mu2*D.y.z[];
 
      ts.z.x[] = mu2*D.z.x[];
      ts.z.y[] = mu2*D.z.y[];
      ts.z.z[] = -p[] + mu2*D.z.z[];
#endif
    }
#endif

#ifdef GNF

#endif
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
