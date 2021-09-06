double sigma = 1./93.0;
double rho1 = 10;
double rho2 = 1;
double mu1 = 3.2791e-2;
double mu2 = 1.037e-2;
#define WIDTH 32
#define dimension 2

#include "view.h"
#include "axi.h"
#include "custom.h"
#include <sys/stat.h>
#include <sys/types.h>

scalar f[];
scalar p[];
scalar omega[];
vector g[];
vector u[];
// scalar umag[];

char dumpname[80];

int main (int argc, char * argv[]) {
  if (argc>1) sprintf (dumpname, "%s", argv[1]);

  FILE * fp = fopen ( "dropstats.txt", "a");

  if (restore (file = dumpname)) {
    fprintf (fp, "Dumpfile restored: %s\n", dumpname);

    if (is_constant(cm)) {
      scalar * l = list_copy (all);
      cm = new scalar;
      free (all);
      all = list_concat ({cm}, l);
      free (l);
    }
    scalar cmv = cm;
    foreach()
      cmv[] = y;
    cm[top] = dirichlet(y);
    cm[bottom] = dirichlet(y);
    
    foreach()
      m[] = f[]>1e-3;
    n_drops = tag(m);
    printf ("%d\n", n_drops);
    DropStats atom[n_drops+1];
    atoms (f, u, atom);

//    tensor D[], vs[];
//    vector vTr[];
//
//    deformation_tensor (u, D);
//    viscous_stress (D, vs);
//    coord F_drop = surfaceforces ( f, vs, vTr);
//
//    scalar vTr_mag[];
//    foreach()
//      vTr_mag[] = sqrt ( sq(vTr.x[]) + sq(vTr.y[]) );

    for (int j = 0; j < n_drops+1; j++)
      if (pid()==0)
        fprintf ( fp, "%.8d %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n",
                  j, 2*pi*atom[j].vol, atom[j].com.x, atom[j].u.x, pi*atom[j].KE,
                  2*pi*atom[j].SE, atom[j].l_max.x, atom[j].l_min.x, atom[j].l_max.y, atom[j].l_min.y );

//    fprintf (fout, "%.8f %.8f\n", F_drop.x, F_drop.y);
//
//    char dirname[80];
//    sprintf (dirname, "./images");
//    mkdir( dirname, 0777);
//    
//    char label[80];
//  
//    clear();
//    view(fov=3,camera="front",tx=-atom[0].com.x/WIDTH, ty=-0.05 );
//    box();
//    draw_vof("f");
//    sprintf (label, "%s/f_%s.ppm", dirname, dumpname);
//    save (label);
//
//    clear();
//    view(fov=3,camera="front",tx=-atom[0].com.x/WIDTH, ty=-0.05 );
//    box();
//    squares( "omega", min=-30, max=30, linear=true, map=cool_warm );
//    draw_vof("f");
//    sprintf (label, "%s/omega_%s.ppm", dirname, dumpname);
//    save (label);
//
//    clear();
//    view(fov=3,camera="front",tx=-atom[0].com.x/WIDTH, ty=-0.05 );
//    box();
//    squares( "vTr_mag", min=0, max=1, linear=true, map=cool_warm );
//    draw_vof("f");
//    sprintf (label, "%s/vTr_%s.ppm", dirname, dumpname);
//    save (label);
//
  }

  else fprintf (ferr, "Dumpfile \"%s\" does not exist!!\n", dumpname );
  fclose(fp);
}
