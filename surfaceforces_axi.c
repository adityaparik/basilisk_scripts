vector vsx[];        // Viscous Stress Tensor: Row 1
vector vsy[];        // Viscous Stress Tensor: Row 2
vector vsz[];        // Viscous Stress Tensor: Row 3

vector tsx[];        // Total Stress Tensor: Row 1
vector tsy[];        // Total Stress Tensor: Row 2
vector tsz[];        // Total Stress Tensor: Row 3

// Viscous Stress Tensor calculation : $\mu_o {\Delta u} + {\Delta u}^T$
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

// Total Stress Tensor Calculation : ts = -p + vs
foreach(){
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

vector Tr[];        // Traction Vector to store stress vector at each interface point

double Fx=0;        // x-component of total force on an object
double Fy=0;        // y-component of total force on an object
double Fz=0;        // z-component of total force on an object
double Fmag=0;      // magnitude of total force on an object

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
