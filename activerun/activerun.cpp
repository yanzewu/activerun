#include  <math.h>
#include  <stdio.h>
#include  <ctime>
#include  <fstream>
#include  <stdlib.h>
#include "functions.h"
#define PI 3.14159
  
/*******************VARIABLES*******************/

double kT, dt, zeta, zetaR;

int steps, output, idx1, idx2, idy1, idy2, ax, ay, ax1, ay1, ax2, ay2, cid1, cid2, type1, type2;

double xmax, ymax, delx, dely, delx2, dely2, rsq, r, dr, fbond, fx, fy, KE, PE,
    temp, pressure, t, cellsize_x, cellsize_y, da, db, stepstR, mincellL;  

int tot_colloids, tot_bonds, nactive, swimstart, spstart, ybin, nsizes, cellcount_x, cellcount_y, limit,
    brownrot, repulsion, attraction, morse, swimattract;

double *x, *y, **Fr, *Lr, *theta, **Fa, **Fs, **Fov, **v, *stressPp, *stressPa, *stress_hardsphere, *stress_swim,
  *stressAp, *stressAa, *stress_pair, *tot_stress, **d, *swim_orig_x, *swim_orig_y;

double tauR, LR, U, aswim, apass, overlap, tolerance, PeR, PeS, ksTs, na, phia, Pi_s, maxcellsize,
    Um, Rg, ddep, D, kappa, cutoff;

int *atom_id, *mol_id, *atom_type, **box_coord, **swim_orig_box, **neightest, ***neighbour;

int main()
{
    int a, count, count2, j, i, k, ii, kk, l, ll, 
        check, check2, check3, check4, check5, 
        intdummy, intdummy2, intdummy3; 
    double dummy, dummy2, dummy3, dummy4, dummy5, dummy6;

    std::clock_t start_time;
    double duration;
    start_time = std::clock();
  
    srand(time(0));

    FILE *fparams;
    
    if ((fparams = fopen( "in.gel", "r" )) == NULL){
        printf ("in.gel file could not be opened\n");
        return 1;
    }

    // READ INPUT FILE

    fscanf(fparams, "%*s%lf", &kT); //thermal energy
    fscanf(fparams, "%*s%lf", &dt); //timestep normalized by shortest relaxation time
    fscanf(fparams, "%*s%lf", &stepstR); //simulation duration in terms of tauR
    fscanf(fparams, "%*s%lf", &zeta); //translational drag factor
    fscanf(fparams, "%*s%lf", &overlap); //maximum distance particles can overlap (units of a)
    fscanf(fparams, "%*s%lf", &PeR); //rotational Peclet number = a/(UtauR)
    fscanf(fparams, "%*s%i", &brownrot); //set to integer != 0 to include Brownian Rotation
    fscanf(fparams, "%*s%lf", &PeS); //specify if brownrot == 0
    fscanf(fparams, "%*s%i", &repulsion); //set to integer != 0 to include HS or Morse repulsion
    fscanf(fparams, "%*s%i", &attraction); //set to integer != 0 to include AO or Morse attraction
    fscanf(fparams, "%*s%i", &swimattract); //set to integer != 0 for swimmers to also feel attraction
    fscanf(fparams, "%*s%i", &morse);//set to integer != 0 to use the Morse potential for attraction+repulsion/repulsion (if attraction == 0)
    fscanf(fparams, "%*s%lf", &Rg); //sets the cutoff of AO or Morse attraction
    fscanf(fparams, "%*s%lf", &kappa); //Morse screening length (units of 1/a)
    fscanf(fparams, "%*s%lf", &Um); //sets the well depth of the AO or Morse attraction
    fscanf(fparams, "%*s%i", &output); //output frequency (units of timestep)
    fscanf(fparams, "%*s%i", &swimstart); //timestep that active colloids begin to swim
    fscanf(fparams, "%*s%i", &spstart); //timestep that swim pressure measurement begins
    fscanf(fparams, "%*s%lf", &mincellL); //length of unit cell in cell list normalized by largest interaction length
    fclose(fparams);
  
    // PARSE INPUT FILE

    if (attraction != 0 && morse == 0)
        repulsion = 1;
    if (morse != 0 && attraction == 0 && repulsion == 0) {
        repulsion = 1; 
        attraction = 0; 
        swimattract = 0;
    }

    FILE *ftesting = fopen("test_neighbors", "w");
    setbuf(ftesting, NULL);

    // READ DATA FILE

    FILE *data;
    if ( ( data = fopen( "data.gel", "r" ) ) == NULL){
        printf ("data.gel file could not be opened\n");
        return 1;
    }
    fscanf(data, "%*s%*s%*s");
    fscanf(data, "%i%*s", &tot_colloids);
    fscanf(data, "%i%*s%*s", &nsizes);

    d = create_2d_array<double>(nsizes, 4);

    stressPa   = create_1d_array<double>(3);
    stressPp   = create_1d_array<double>(3);
    stress_hardsphere  = create_1d_array<double>(3);
    stress_swim    = create_1d_array<double>(3);
    stressAa   = create_1d_array<double>(3);
    stressAp   = create_1d_array<double>(3);
    stress_pair  = create_1d_array<double>(3);
    tot_stress = create_1d_array<double>(3);
    for (int i = 0; i < 3; ++i) {
        stress_swim[i] = 0; 
        stressPa[i] = 0; 
        stressPp[i] = 0; 
        stress_hardsphere[i] = 0; 
        stressAa[i] = 0; 
        stressAp[i] = 0; 
        stress_pair[i] = 0; 
    }

    atom_id   = create_1d_array<int>(tot_colloids);
    mol_id    = create_1d_array<int>(tot_colloids);
    atom_type = create_1d_array<int>(tot_colloids);

    x         = create_1d_array<double>(tot_colloids);
    y         = create_1d_array<double>(tot_colloids);
    swim_orig_x        = create_1d_array<double>(tot_colloids);
    swim_orig_y        = create_1d_array<double>(tot_colloids);
    box_coord     = create_2d_array<int>(tot_colloids, 2);
    for(int i = 0; i < tot_colloids; i++){
        box_coord[j][0] = 0;
        box_coord[j][1] = 0;
    }

    swim_orig_box    = create_2d_array<int>(tot_colloids, 2);
    Lr        = create_1d_array<double>(tot_colloids);
    theta     = create_1d_array<double>(tot_colloids);
    Fr        = create_2d_array<double>(tot_colloids, 2);
    Fa        = create_2d_array<double>(tot_colloids, 2);
    Fs        = create_2d_array<double>(tot_colloids, 2);
    Fov       = create_2d_array<double>(tot_colloids, 2);
    v         = create_2d_array<double>(tot_colloids, 2);
    for(int i = 0; i < tot_colloids; i++){
        v[j][0] = 0; 
        v[j][1] = 0;
    }

    fscanf(data, "%*lf%lf%*s%*s", &xmax);
    fscanf(data, "%*lf%lf%*s%*s", &ymax);
    fscanf(data, "%*lf%*lf%*s%*s");
    fscanf(data, "%*lf%*lf%*lf%*s%*s%*s");
    fscanf(data, "%*s");
    for (j = 0; j < nsizes; ++j){
        fscanf(data, "%*i%*lf");
    }

    // READING ATOM TYPE

    fscanf(data, "%*s%*s");
    maxcellsize = 0;
    for (j = 0; j < nsizes; ++j) {
        int atom_idx;
        double atom_size;
        fscanf(data, "%i%lf", &atom_idx, &atom_size);
        d[atom_idx - 1][0] = atom_size;
        maxcellsize = atom_size > maxcellsize ? atom_size : maxcellsize;
    }
    aswim = d[0][0]/2; //swimmer radius
    apass = d[1][0]/2; //average passive colloid radius

    // READING ATOM DATA PART

    fscanf(data, "%*s");
    nactive = 0;
    for (j = 0; j < tot_colloids; ++j) {
        fscanf(data, "%i%i%i%lf%lf%*i", &atom_id[j], &mol_id[j], &atom_type[j], &x[j], &y[j]);

        // counting active particles
        if (atom_type[j] == 1){
          nactive++;
        }
        // fix improper position

        while (x[j] < 0) { box_coord[j][0]--; x[j] += xmax; }
        while (x[j] > xmax) { box_coord[j][0]++; x[j] -= xmax; }
        
        while (y[j] < 0) { box_coord[j][1]--; y[j] += ymax; }
        while (y[j] > ymax) { box_coord[j][1]++; y[j] -= ymax; } 
    }
    fclose(data);

    // SETTING UP ROTATIONAL DATA

    zetaR = zeta*aswim*aswim*4/3;
    tauR = brownrot ? zetaR/kT : aswim*aswim*zeta/kT/PeR/PeS;

    for (j = 0; j < nsizes; ++j) {
        d[j][1] = zeta/apass*d[j][0];
        d[j][2] = zetaR/(apass*apass*apass)*d[j][0]*d[j][0]*d[j][0]/8;
        d[j][3] = brownrot ? d[j][2]/kT : tauR;
    }

    if (nactive != 0){ 
        LR = aswim/PeR;
        U = LR/tauR;
        ksTs = U*U*tauR/2;
        if (brownrot != 0)
            PeS = U*aswim*zeta/kT;
        na = nactive/xmax/ymax;
        phia = na*PI*aswim*aswim;
        Pi_s = na*ksTs;
    
        double dt_scale;
        if (aswim*aswim*zeta/kT < aswim/U)
            dt_scale = aswim*aswim*zeta/kT;
        else if (aswim/U < tauR) {
            dt_scale = aswim/U;
        }
        else {
            dt_scale = tauR;
        }

        dt *= dt_scale;
        steps = stepstR * tauR/dt;
    
        // output for checking
        FILE *swimproperties;
        swimproperties = fopen("swimproperties", "w");
        fprintf(swimproperties, "tauR:\t%f\nLR:\t%f\nU:\t%f\nPe_s:\t%f\nPe_R:\t%f\n", tauR, LR, U, PeS, PeR);
        fprintf(swimproperties, "ksTs:\t%f\nna:\t%f\nphia:\t%f\nPi_s:\t%f\n", ksTs, na, phia, Pi_s);
        fclose(swimproperties);
    }
    else {
        dt *= (apass*apass*zeta/kT);
        steps = stepstR * (apass*apass*zeta/kT)/dt;
    }

    if (attraction != 0) {
        ddep = Rg*apass*2;
        maxcellsize += ddep;
    }
    maxcellsize *= mincellL;

    cellcount_x = xmax/maxcellsize;
    cellcount_y = ymax/maxcellsize;
    cellsize_x = xmax/cellcount_x;
    cellsize_y = ymax/cellcount_y;


    neighbour      = create_3d_array<int>(cellcount_x, cellcount_y, tot_colloids+1);
    //! MEMSET HERE!
    for (j = 0; j < cellcount_x; ++j) {
        for (k = 0; k < cellcount_y; ++k) {
            neighbour[j][k][0] = 0;
        }
    }   
    //neightest  = create_2d_array<int>(cellcount_x, cellcount_y);


    for (j = 0; j < tot_colloids; ++j) {
        int cellid_x  = x[j]/cellsize_x;
        int cellid_y = y[j]/cellsize_y;
        int neighbour_count = neighbour[cellid_x][cellid_y][0];
        neighbour[cellid_x][cellid_y][neighbour_count +1] = j;
        neighbour[cellid_x][cellid_y][0]++;
    }

    // random angle
    for (j = 0; j < tot_colloids; ++j) {
        if (atom_type[i] == 1){
            theta[i] = fRand(0, 2*PI);
        } 
    }

/*********************SIMULATION***********************/
 
    // OUTPUT INITIAL STATE

    FILE *trajec;
    trajec = fopen("outputconfig.lammpstrj", "w");
    setbuf(trajec, NULL);
    fprintf(trajec, "ITEM: TIMESTEP\n");
    fprintf(trajec, "0\n"); 
    fprintf(trajec, "ITEM: NUMBER OF ATOMS\n");
    fprintf(trajec, "%i\n", tot_colloids);
    fprintf(trajec, "ITEM: BOX BOUNDS\n");
    fprintf(trajec, "0 %f\n", xmax);
    fprintf(trajec, "0 %f\n", ymax);
    fprintf(trajec, "-0.25 0.25\n");
    fprintf(trajec, "ITEM: ATOMS id mol type x y z\n");
    for (int i = 0; i < tot_colloids; ++i){
        fprintf(trajec, "%i %i %i %f %f 0\n", i+1, mol_id[i], atom_type[i], x[i], y[i]);
    }
  
    FILE *swimstress, *totalstress;
    swimstress  = fopen("swimstress.output", "w");
    totalstress = fopen("totalstress.output", "w");
    fprintf(swimstress, "Step\tPixx/Pi_s\tPiyy/Pi_s\tPi/Pi_s\n");
    fprintf(totalstress, "Step\tPkt\tPswim\tPHS\tPpot\tPtot\n");


    // the main loop

    while( j < steps ) {


    //Compute Potential Forces
        if (attraction != 0 || morse != 0) {
            for (i = 0; i < cellcount_x; ++i) {
                for (k = 0; k < cellcount_y; ++k) {
                    for (l = 1; l<= neighbour[i][k][0]; ++l) {
                        for (ii = i; ii <= (i+1); ++ii) {
                            for (kk = (k-1); kk <= (k+1); ++kk) {
                                if (ii != i || kk != (k-1)) {
                                    if (ii == cellcount_x) {
                                      idx2 = 0;
                                      ax2 = 1;
                                    }	
                                    else {
                                      idx2 = ii;
                                      ax2 = 0;
                                    }
                                    if (kk == cellcount_y) {
                                      idy2 = 0;
                                      ay2 = 1;
                                    }	
                                    else if (kk == -1)  {
                                      idy2 = cellcount_x - 1;
                                      ay2 = -1;
                                    }
                                    else {
                                      idy2 = kk;
                                      ay2 = 0;
                                    }

                                    if (i == ii && k == kk) limit = (l+1);
                                    else limit = 1;

                                    for (ll = limit; ll<= neighbour[idx2][idy2][0]; ++ll) {
                                      cid1 = neighbour[i][k][l];
                                      cid2 = neighbour[idx2][idy2][ll];
                                      delx = (x[cid2] + ax2*xmax) - (x[cid1]);
                                      dely = (y[cid2] + ay2*ymax) - (y[cid1]);
                                      r = sqrt(delx*delx + dely*dely);
                                  
                                      type1 = atom_type[cid1] - 1;
                                      da = d[type1][0];
                                      type2 = atom_type[cid2] - 1;
                                      db = d[type2][0];
                                      D = (da + db)/2;

                                        if ((type1 == 0 || type2 == 0) && swimattract == 0)
                                            cutoff = 0;
                                        else
                                            cutoff = D + ddep;
                                    
                                        if (morse == 0) {
                                            // ?? potential?
                                            if (r < cutoff && r > D) {		      
                                                dummy  = -delx/r*2*Um/ddep/ddep*(D + ddep - r);
                                                dummy2 = -dely/r*2*Um/ddep/ddep*(D + ddep - r);

                                                Fa[cid2][0] += dummy;
                                                Fa[cid2][1] += dummy2;
                                                Fa[cid1][0] -= dummy;
                                                Fa[cid1][1] -= dummy2;

                                                dummy = dummy/xmax/ymax*delx;
                                                dummy2 = dummy2/xmax/ymax*dely;
                                    
                                                stress_pair[0] += dummy;
                                                stress_pair[1] += dummy2;
                                                stress_pair[2] += (dummy/2 + dummy2/2);
                                            }
                                         }

                                        // morse potential
                                        else {
                                            if (attraction == 0 || (type1 == 0 || type2 == 0) && swimattract == 0)
                                                cutoff = D;
                                    
                                            if (r < cutoff) {
                                                double fx  = -delx/r*2*kappa*Um*(exp(-kappa*(r-D)) - exp(-2*kappa*(r-D)));
                                                double fy = -dely/r*2*kappa*Um*(exp(-kappa*(r-D)) - exp(-2*kappa*(r-D)));
                                        //			dummy = 0; dummy2 = 0;
                                                Fa[cid2][0] += fx;
                                                Fa[cid2][1] += fy;
                                                Fa[cid1][0] -= fx;
                                                Fa[cid1][1] -= fy;
                                
                                                stress_pair[0] += fx / xmax / ymax*delx;
                                                stress_pair[1] += fy / xmax / ymax*dely;
                                                stress_pair[2] += (stress_pair[0] + stress_pair[1]) / 2;
                                            }
                                        }
                                    }
                              }
	                        }
                      }
                    }
	              }
            }
        }
    
        // record swim position?
        if (j - swimstart == spstart) {
            for (int i = 0; i < tot_colloids; ++i) {
                swim_orig_x[i] = x[i];
                swim_orig_y[i] = y[i];
                swim_orig_box[i][0] = box_coord[i][0];
                swim_orig_box[i][1] = box_coord[i][1];
            }
        }

        // random 
        for (i = 0; i < tot_colloids; ++i) {
            Fr[i][0] = sqrt(24 * kT*d[atom_type[i] - 1][1] / dt)*(fRand(0, 1) - 0.5);
            Fr[i][1] = sqrt(24 * kT*d[atom_type[i] - 1][1] / dt)*(fRand(0, 1) - 0.5);
        }

        //Compute Swim Force
        if (j >= swimstart) {
            for(int i = 0; i < tot_colloids; ++i){
                if (atom_type[i] == 1) {
	                Lr[i] = sqrt(24*d[type1][2]*d[type1][2]/d[type1][3]/dt)*(fRand(0,1)-0.5);
                    theta[i] = theta[i] + Lr[i]/d[type1][2]*dt;
                    Fs[i][0] = d[type1][1]*U*cos(theta[i]);
                    Fs[i][1] = d[type1][1]*U*sin(theta[i]);
                }
            }
        }

        // Integrator
        for(int i = 0; i < tot_colloids; i++)
            v[i][0] = 1/d[type1][1]*(Fr[i][0] + Fs[i][0] + Fa[i][0]);
            v[i][1] = 1/d[type1][1]*(Fr[i][1] + Fs[i][1] + Fa[i][1]);

            y[i] = y[i] + dt*v[i][1];
            x[i] = x[i] + dt*v[i][0];
        }

        /*O(N^2) Overlap Resolution ftesting
            // check = 1;
            // while (check != 0) {
            //   check = 0;
            //   for (l = 0; l < tot_colloids; ++l) {
            // 	for (ll=(l+1); ll < tot_colloids; ++ll) {

            // 	  delx = x[l] - x[ll];
            // 	  if (delx <  -xmax/2)
            // 	    delx += xmax;
            // 	  if (delx > xmax/2)
            // 	    delx -= xmax;
            
            // 	  dely = y[l] - y[ll];
            // 	  if (dely < -xmax/2) 
            // 	    dely += ymax;
            // 	  if (dely > ymax/2)
            // 	    dely -= ymax;
            
            // 	  type1 = atom_type[l] - 1;
            // 	  da = d[type1][0];
            // 	  type2 = atom_type[ll] - 1;
            // 	  db = d[type2][0];
            
            // 	  rsq = delx*delx + dely*dely;
            // 	  r = sqrt(rsq);
            
            // 	  tolerance = overlap*(da + db)/2;
            // 	  if (r < ((da + db)/2-tolerance)) {
            // 	    delx2 = ( (da + db)/2 - r)/2*delx/r;
            // 	    dely2 = ( (da + db)/2 - r)/2*dely/r;
            // 	    x[l]  = x[l]  + delx2;
            // 	    x[ll] = x[ll] - delx2;
            // 	    y[l]  = y[l]  + dely2;
            // 	    y[ll] = y[ll] - dely2;
            // 	    check++;
            // 	  }
            // 	}
            //   }
        // }*/
    
        //Resolve Disk Overlaps with Cell List O(N) using Potential Free Algorithm
        check = 1;
        while (check != 0 && repulsion != 0 && morse == 0) {
            check = 0;
            for (i = 0; i < cellcount_x; ++i) {
    	          for (k = 0; k < cellcount_y; ++k) {
                    for (l = 1; l<= neighbour[i][k][0]; ++l) {
                      for (ii = i; ii<=(i+1); ++ii) {
                          for (kk = (k-1); kk<=(k+1); ++kk) {
                              if (ii != i || kk != (k-1)) {
                                  if (ii== cellcount_x) {
                                      idx2 = 0;
                                      ax2 = 1;
                                  }	
                                  else {
                                      idx2 = ii;
                                      ax2 = 0;
                                  }
                                  if (kk== cellcount_y) {
                                    idy2 = 0;
                                    ay2 = 1;
                                  }	
                                  else if (kk == -1)  {
                                    idy2 = cellcount_x - 1;
                                    ay2 = -1;
                                  }
                                  else {
                                    idy2 = kk;
                                    ay2 = 0;
                                  }
                                  if (i == ii && k == kk)
                                    limit = (l+1);
                                  else
                                    limit = 1;
                                  for (ll = limit; ll<= neighbour[idx2][idy2][0]; ++ll) {
                                    cid1 = neighbour[i][k][l];
                                    cid2 = neighbour[idx2][idy2][ll];
                                    delx = (x[cid2] + ax2*xmax) - (x[cid1]);
                                    dely = (y[cid2] + ay2*ymax) - (y[cid1]);
                                    
                                    type1 = atom_type[cid2] - 1;
                                    da = d[type1][0];
                                    type2 = atom_type[cid2] - 1;
                                    db = d[type2][0];
                                    
                                    rsq = delx*delx + dely*dely;
                                    r = sqrt(rsq);
                                    tolerance = overlap*(da + db)/2;
                                    if (r < (da/2 + db/2 - tolerance)) {
                                      delx2 = ( (da + db)/2 - r)/2*delx/r;
                                      dely2 = ( (da + db)/2 - r)/2*dely/r;
                                      x[cid2] += delx2;
                                      x[cid1] -= delx2;
                                      y[cid2] += dely2;
                                      y[cid1] -= dely2;
                                      delx = (x[cid2] + ax2*xmax) - (x[cid1]);
                                      dely = (y[cid2] + ay2*ymax) - (y[cid1]);
                                      stress_hardsphere[0] += (delx*delx2/dt*d[type1][1] + delx*delx2/dt*d[type2][1])/2/xmax/ymax;
                                      stress_hardsphere[1] += (dely*dely2/dt*d[type1][1] + dely*dely2/dt*d[type1][1])/2/xmax/ymax;

                                      //     fprintf(ftesting, "%f\t%f\n", r, sqrt(delx*delx + dely*dely));

                                      stress_hardsphere[2] += (stress_hardsphere[0] /2 + stress_hardsphere[1] /2);
                                      check++;
                                  }
                              }
                          }
                      }
                  }
                }
	          }
        }
    

        //Check interparticle distance 
        check = 0;
        for (l = 0; l < tot_colloids; ++l) {
            for (ll = (l+1); ll < tot_colloids; ++ll) {
                delx = x[l] - x[ll];
                if (delx < -xmax/2) delx += xmax;
                else if (delx > xmax/2) delx -= xmax;

                dely = y[l] - y[ll];
                if (dely < -xmax/2) dely += ymax;
                else if (dely > ymax/2) dely -= ymax;

                rsq = delx*delx + dely*dely;
                r = sqrt(rsq);
                
                type1 = atom_type[l] - 1;
                da = d[type1][0];
                type2 = atom_type[ll] - 1;
                db = d[type2][0];
                tolerance = overlap*(da + db)/2;
                if (r < (da/2 + db/2 - tolerance)) { 		
                    check++;
                    dummy += r;
                }
            }
        }
    
        fprintf(ftesting, "%i\t%i\t%f\n", j+1, check, check ? dummy/check : 0);
 
        //! CHANGE TO MEMSET HERE
        for (i = 0; i < cellcount_x; ++i) {
            for (k = 0; k < cellcount_y; ++k) {
	              for (l = 0; l<= neighbour[i][k][0]; ++l){
                    neighbour[i][k][l] = 0;
                }
            }
        }
    
        // PBC
        for (int i = 0; i < tot_colloids; ++i) {
            if (y[i] < 0) { box_coord[i][1]--; y[i] += ymax;}
            else if (y[i] > ymax) { box_coord[i][1]++; y[i] -= ymax;}
             
            if (x[i] < 0) { box_coord[i][0]--; x[i] += xmax;}
            else if (x[i] > xmax) { box_coord[i][0]++; x[i] -= xmax;}

            intdummy  = x[i]/cellsize_x;
            intdummy2 = y[i]/cellsize_y;
            count = neighbour[intdummy][intdummy2][0];
            neighbour[intdummy][intdummy2][count+1] = i;
            neighbour[intdummy][intdummy2][0]++;
            
        }

        // stress
        if ((j - swimstart) >=  spstart) {
            for (int i = 0; i < tot_colloids; i++) {
                stress_swim[0] += Fs[i][0]*(x[i] - swim_orig_x[i] + xmax*box_coord[i][0] - xmax*swim_orig_box[i][0])/xmax/ymax;
                stress_swim[1] += Fs[i][1]*(y[i] - swim_orig_y[i] + ymax*box_coord[i][1] - ymax*swim_orig_box[i][1])/xmax/ymax;
                stress_swim[2] += (stress_swim[0] + stress_swim[1]) / 2;
            }
        }

    //***********************OUTPUT*************************//


        if ((j+1) % output == 0) {

            // DUMP HEAD
            fprintf(trajec, "ITEM: TIMESTEP\n");
            fprintf(trajec, "%i\n", j+1); 
            fprintf(trajec, "ITEM: NUMBER OF ATOMS\n");
            fprintf(trajec, "%i\n", tot_colloids);
            fprintf(trajec, "ITEM: BOX BOUNDS\n");
            fprintf(trajec, "0 %f\n", xmax);
            fprintf(trajec, "0 %f\n", ymax);
            fprintf(trajec, "-0.25 0.25\n");
            //      fprintf(trajec, "ITEM: ATOMS id mol type x y vx vy Fovx Fovy ix iy iz\n");
            fprintf(trajec, "ITEM: ATOMS id mol type x y z\n");

            for (i = 0; i < tot_colloids; ++i) {
	              //fprintf(trajec, "%i %i %i %f %f 0 %i %i 0\n", i+1, mol_id[i], atom_type[i], x[i], y[i], box_coord[i][0], box_coord[i][1]);
	              fprintf(trajec, "%i %i %i %f %f 0\n", i+1, mol_id[i], atom_type[i], x[i], y[i]);
            }
            
            fprintf(swimstress, "%i\t", (j+1));
            fprintf(totalstress, "%i\t", (j+1));      
            fprintf(totalstress, "%f\t%f\t%f\t%f\t%f\n", 
                (tot_colloids/xmax/ymax*kT), 
                (stress_swim[2]/output), 
                (stress_hardsphere[2]/output), 
                (stress_pair[2]/output), 
                (tot_colloids/xmax/ymax*kT + (stress_swim[2]+stress_hardsphere[2]+stress_pair[2])/output)
            );
            for (i = 0; i < 3; ++i) {
	              fprintf(swimstress, "%f\t", (stress_swim[i]/output/Pi_s));
            }
            fprintf(swimstress, "\n");
        }

        // CLEAR UP

        for(int i = 0; i < 3; i++){
            stress_swim[i] = 0; 
            stress_hardsphere[i] = 0; 
            stress_pair[i] = 0;
        }

        for (int i = 0; i < tot_colloids; i++){
            Fa[i][0] = 0;
            Fa[i][1] = 0;
        }

        j++;
    }

    duration = ( std::clock() - start_time )/(double) CLOCKS_PER_SEC;
  
    FILE *t;
    t = fopen("time.output", "w");
    fprintf(t, "Steps:\t%i\n", steps);
    fprintf(t, "Timestep:\t%f\n", dt);
    fprintf(t, "Output:\t%i\n", output);
    fprintf(t, "Duration:\t%f\n", duration);
    return 0;
}
