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

typedef struct {
  double vol;
  double ke;        // Droplet absolute KE
//  double rel_ke=0.;
  double s_en;           // Droplet Surface Energy
  double xmax;
  double xmin;
  double ymax;
  double ymin;
#if dimension==3
  double zmax;
  double zmin;
#endif
  coord cm;
  coord vel;
} DropStats;

typedef struct {
  Dropstats *ptr = NULL;
  Dropstats bulk;
  int n_drops;
} Fragments;

DropStats dropstats ( const scalar f, const vector u ) {
  Fragments drops;
  scalar m[];
  foreach()
    m[] = f[] > 1e-3;
  drops.n_drops = tag(m);

  static DropStats atom[n_drops];
  drops.ptr = atom;

// Initialize the stat variables to zero  
  for (int j = 0; j < n_drops; j++) {
    atom[j].vol;
    atom[j].ke=0.;        // Droplet absolute KE
  //  double rel_ke=0.;
    atom[j].se=0.;           // Droplet Surface Energy
    atom[j].xmax=0;
    atom[j].xmin=WIDTH;
    atom[j].ymax=0;
    atom[j].ymin=WIDTH;
  #if dimension==3
    atom[j].zmax=0;
    atom[j].zmin=WIDTH;
    atom[j].cm={0.,0.,0.};
    atom[j].vel={0.,0.,0.};
  #else
    atom[j].cm={0.,0.};
    atom[j].vel={0.,0.};
  #endif
  }

  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      atom[j].vol += dv()*f[];
      coord p = {x,y,z};
      foreach_dimension()
        atom[j].cm.x += dv()*f[]*p.x;
      
      atom[j].vel.x += u.x[]*dv()*f[];
      atom[j].vel.y += u.y[]*dv()*f[];
    }

  MPI_Allreduce (MPI_IN_PLACE, atom, n_drops, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (int j = 0; j < n_drops; j++) {
    drops.bulk.cm.x += b[j].x;
    bulk.cm.y += b[j].y;
    bulk.cm.z += b[j].z;
    bulk.vol += vol[j];
  }

  bulk_com[0] = bulk_com[0]/bulk_vol;    // Bulk COM for camera location
  bulk_com[1] = bulk_com[1]/bulk_vol;
  bulk_com[2] = bulk_com[2]/bulk_vol;
  
  double xmin = WIDTH, xmax = 0.;
  double ymin = WIDTH, ymax = 0.;

  foreach( reduction(min:xmin) reduction(max:xmax) reduction(min:ymin) reduction(max:ymax) ) {
    if( f[]>0.01 ) {
      if ( fabs(y) < 0.02 ) {
        if (x>xmax)  xmax = x;
        if (x<xmin)  xmin = x;
      }
      if (y>ymax)  ymax = y;
      if (y<ymin)  ymin = y;
    }
  }

  if (i==0)
    fprintf ( fout, "iteration_nos\ttime\tdroplet_index\tvolume\tx_loc\ty_loc\tz_loc\tx_vel\ty_vel\txy_AR\n" );
  
  int j = 0;
  
  fprintf ( fout, "%.8d %.8f %.8d %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n",
            i, t, j, vol[j], b[j].x/vol[j], b[j].y/vol[j], b[j].z/vol[j], velx[j]/vol[j], vely[j]/vol[j], 0.5*(xmax-xmin)/(ymax-ymin) );

  for (j = 1; j < n_drops; j++)
       fprintf ( fout, "%.8d %.8f %.8d %.8f %.8f %.8f %.8f %.8f %.8f\n",
              i, t, j, vol[j], b[j].x/vol[j], b[j].y/vol[j], b[j].z/vol[j], velx[j]/vol[j], vely[j]/vol[j] );
  
  fflush (fout);
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
