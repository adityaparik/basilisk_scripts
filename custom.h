// This header file stores many custom functions
// If using Stress Tensor calculation functions, you must define if 
// you wish to use a specific type of Constitutive relationship.

#include "utils.h"
#include "fractions.h"
#include "tag.h"

double interfacearea (scalar f) {
  double area = 0.;
  foreach( reduction(+:area) ) {
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = interface_normal (point, f), p;
      double alpha = plane_alpha (f[], n);
#if dimension==3
      area += pow(Delta,2) * plane_area_center (n, alpha, &p);
#else
      area += cm[]*Delta * plane_area_center (n, alpha, &p);
#endif
    }    
  }
  return area;
}

// You must use tag function to obtain n_drops and m[]
int n_drops;
scalar m[];

typedef struct {
  double vol;
  double KE;        // Droplet absolute KE
  double SE;      // Droplet Surface Energy
  coord com;
  coord u;
  coord l_max;      // Maximum extent of the droplet
  coord l_min;      // Minimum extent of the droplet
} DropStats;

void atoms ( const scalar f, const vector u, DropStats *atom ) {
  // Initialize the all atom array components to zero  
  for (int j = 0; j < n_drops+1; j++) {
    atom[j].vol = atom[j].KE = atom[j].SE=0.;
    atom[j].com.x = 0.; atom[j].com.y = 0.;
    atom[j].u.x = 0.; atom[j].u.y = 0.;
    atom[j].l_max.x=-WIDTH; atom[j].l_max.y=-WIDTH;
    atom[j].l_min.x=WIDTH; atom[j].l_min.y=WIDTH;
#if dimension==3
    atom[j].com.z = 0.;
    atom[j].u.z = 0.;
    atom[j].l_max.z=-WIDTH;
    atom[j].l_min.z=WIDTH;
#endif
  }

  foreach_leaf() {

    if (m[] > 0) {
      int j = m[];
      atom[j].vol += dv()*f[];      // dv() is missing 2pi factor 
      coord p = {x,y,z};

// Surface Energy Calculation for the current cell
      if (f[] > 1e-6 && f[] < 1. - 1e-6) {
        coord n = interface_normal (point, f), p;
        double alpha = plane_alpha (f[], n);
#if dimension==3
        atom[j].SE += pow(Delta,2)*plane_area_center (n, alpha, &p);
#else
// Axisymmetric domain has a 2pi factor missing
        atom[j].SE += cm[]*Delta*plane_area_center (n, alpha, &p);
#endif
      }

      foreach_dimension() {
        atom[j].KE += dv()*f[] * sq(u.x[]); // dv() is missing 2pi factor for AXI
        atom[j].com.x += dv()*f[] * p.x;
        atom[j].u.x += dv()*f[] * u.x[];
      }
    }
  }
//  MPI_Allreduce (MPI_IN_PLACE, atom, n_drops+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  for (int j = 1; j < n_drops+1; j++) {
    atom[0].vol += atom[j].vol;

// Assigning Bulk stats before atom fragment stats are modified 
    atom[0].com.x += atom[j].com.x;
    atom[0].com.y += atom[j].com.y;
    atom[0].u.x += atom[j].u.x;
    atom[0].u.y += atom[j].u.y;

    atom[j].com.x = atom[j].com.x/atom[j].vol;
    atom[j].com.y = atom[j].com.y/atom[j].vol;
    atom[j].u.x = atom[j].u.x/atom[j].vol;
    atom[j].u.y = atom[j].u.y/atom[j].vol;
#if dimension==3
    atom[0].com.z += atom[j].com.z;
    atom[0].u.z += atom[j].u.z;

    atom[j].com.z = atom[j].com.z/atom[j].vol;
    atom[j].u.z = atom[j].u.z/atom[j].vol;
#endif
    atom[j].SE = sigma*atom[j].SE;   // Missing 2pi factor for AXI
    atom[j].KE = rho1*atom[j].KE;      // Missing 2pi factor for AXI

    atom[0].SE += atom[j].SE;          // Missing 2pi factor for AXI
    atom[0].KE += atom[j].KE;          // Missing 2pi factor for AXI
  }

// Assigning Bulk Stats
  atom[0].com.x = atom[0].com.x/atom[0].vol;
  atom[0].com.y = atom[0].com.y/atom[0].vol;
  atom[0].u.x = atom[0].u.x/atom[0].vol;
  atom[0].u.y = atom[0].u.y/atom[0].vol;
#if dimension==3
  atom[0].com.z = atom[0].com.z/atom[0].vol;
  atom[0].u.z = atom[0].u.z/atom[0].vol;
#endif

  foreach_leaf() {
    if (m[] > 0) {
      int j = m[];
      if (x<atom[j].l_min.x)  atom[j].l_min.x = x;
      if (y<atom[j].l_min.y)  atom[j].l_min.y = y;
#if dimension==3
      if (z<atom[j].l_min.z)  atom[j].l_min.z = z;
#endif
    }
  }
//  MPI_Allreduce (MPI_IN_PLACE, atom, n_drops+1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  foreach_leaf() {
    if (m[] > 0) {
      int j = m[];
      if (x>atom[j].l_max.x)  atom[j].l_max.x = x;
      if (y>atom[j].l_max.y)  atom[j].l_max.y = y;
#if dimension==3
      if (z>atom[j].l_max.z)  atom[j].l_max.z = z;
#endif
    }
  }
//  MPI_Allreduce (MPI_IN_PLACE, atom, n_drops+1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  for (int j = 1; j < n_drops+1; j++) {
    if (atom[j].l_max.x>atom[0].l_max.x)  atom[0].l_max.x = atom[j].l_max.x;
    if (atom[j].l_max.y>atom[0].l_max.y)  atom[0].l_max.y = atom[j].l_max.y;
    if (atom[j].l_max.z>atom[0].l_max.z)  atom[0].l_max.z = atom[j].l_max.z;

    if (atom[j].l_min.x<atom[0].l_min.x)  atom[0].l_min.x = atom[j].l_min.x;
    if (atom[j].l_min.y<atom[0].l_min.y)  atom[0].l_min.y = atom[j].l_min.y;
    if (atom[j].l_min.z<atom[0].l_min.z)  atom[0].l_min.z = atom[j].l_min.z;
  }

} // End of atoms function


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
  
      D.z.x[] = (u.z[1,0,0] - u.z[-1,0,0])/(2.*Delta) + (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
      D.z.y[] = (u.z[0,1,0] - u.z[0,-1,0])/(2.*Delta) + (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
      D.z.z[] = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta) + (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);   
#endif
  }
}

// Function to calculate viscous stress Tensor
void viscous_stress (const tensor D, tensor vs) {
#ifdef NEWTONIAN
// Newtonian viscous stress Tensor Calculation : vs = mu*D 
  foreach() {
      vs.x.x[] = mu2*D.x.x[];
      vs.x.y[] = mu2*D.x.y[];
      vs.y.x[] = mu2*D.y.x[];
      vs.y.y[] = mu2*D.y.y[];
#if dimension==3
      vs.x.z[] = mu2*D.x.z[];
      vs.y.z[] = mu2*D.y.z[];
      vs.z.x[] = mu2*D.z.x[];
      vs.z.y[] = mu2*D.z.y[];
      vs.z.z[] = mu2*D.z.z[];
#endif
  } 
#endif
}

// Function to calculate total stress Tensor
void total_stress (const scalar p, const tensor vs, tensor ts) {
  foreach() {
      ts.x.x[] = -p[] + vs.x.x[];
      ts.x.y[] = vs.x.y[];
#if dimension==3
      ts.x.z[] = vs.x.z[];
#endif
      ts.y.x[] = vs.y.x[];
      ts.y.y[] = -p[] + vs.y.y[];
#if dimension==3
      ts.y.z[] = vs.y.z[];
 
      ts.z.x[] = vs.z.x[];
      ts.z.y[] = vs.z.y[];
      ts.z.z[] = -p[] + vs.z.z[];
#endif
  }
}

// Function to obtain surface stresses for a given stress tensor

coord surfaceforces (const scalar f, const tensor ts, vector Tr) {
#if dimension==3
  coord F_drop = {0.,0.,0.};
#else
  coord F_drop = {0.,0.};
#endif
  foreach() {
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      // VOF interface segment outward normal, alpha calculation
      coord m = mycs (point, f);
      double alpha = plane_alpha (f[], m);
      coord pp;       // Variable to store interface segment centroid
      
      Tr.x[] = ( ts.x.x[]*m.x + ts.x.y[]*m.y + ts.x.z[]*m.z );
      Tr.y[] = ( ts.y.x[]*m.x + ts.y.y[]*m.y + ts.y.z[]*m.z );
#if dimension==3
      Tr.z[] = ( ts.z.x[]*m.x + ts.z.y[]*m.y + ts.z.z[]*m.z );
      // Area of VOF interface segment for 3D domain
      double segmentarea = Delta*Delta*plane_area_center(m, alpha, &pp);
#else
      // Area of VOF interface segment for AXI or 2D
      double segmentarea = cm[]*Delta*plane_area_center(m, alpha, &pp);
#endif
     F_drop.x += Tr.x[] * segmentarea;
     F_drop.y += Tr.y[] * segmentarea;
#if dimension==3
     F_drop.z += Tr.z[] * segmentarea;
#endif
    }
  } 

  return F_drop;
}

