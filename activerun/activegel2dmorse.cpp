#include <math.h>
#include <stdio.h>
#include <ctime>
#include <fstream>
#include <stdlib.h>
#include "functions.h"
#define PI 3.14159
  
/*******************VARIABLES*******************/

double kT, dt, zeta, zetaR;

int steps, output, idx1, idx2, idy1, idy2, ax, ay, ax1, ay1, ax2, ay2, cid1, cid2, type1, type2;

double xmax, ymax, delx, dely, delx2, dely2, rsq, r, dr, fbond, fx, fy, KE, PE,
  temp, pressure, t, dnx, dny, da, db, stepstR, mincellL;  

int tot_colloids, tot_bonds, nactive, swimstart, spstart, ybin, nsizes, nx, ny, limit,
  brownrot, repulsion, attraction, morse, swimattract;

double *x, *y, **Fr, *Lr, *theta, **Fa, **Fs, **Fov, **v, *stressPp, *stressPa, *stressPpa, *stressS,
  *stressAp, *stressAa, *stressApa, *tot_stress, **d, *xi, *yi;

double tauR, LR, U, aswim, apass, overlap, tolerance, PeR, PeS, ksTs, na, phia, Pi_s, maxL,
  Um, Rg, ddep, D, kappa, cutoff;

int *atom_id, *mol_id, *atom_type, **image, **imagei, **neightest, ***neigh;

int main()
{
  int a, count, count2, j, i, k, ii, kk, l, ll, 
    check, check2, check3, check4, check5, 
    intdummy, intdummy2, intdummy3; 
  double dummy, dummy2, dummy3, dummy4, dummy5, dummy6;

  std::clock_t start;
  double duration;
  start = std::clock();
  
  srand(time(0));

  FILE *params;
    
  if ( ( params = fopen( "in.gel", "r" ) ) == NULL)
    printf ("in.gel file could not be opened\n");
  else {
    fscanf(params, "%*s%lf", &kT); //thermal energy
    fscanf(params, "%*s%lf", &dt); //timestep normalized by shortest relaxation time
    fscanf(params, "%*s%lf", &stepstR); //simulation duration in terms of tauR
    fscanf(params, "%*s%lf", &zeta); //translational drag factor
    fscanf(params, "%*s%lf", &overlap); //maximum distance particles can overlap (units of a)
    fscanf(params, "%*s%lf", &PeR); //rotational Peclet number = a/(UtauR)
    fscanf(params, "%*s%i", &brownrot); //set to integer != 0 to include Brownian Rotation
    fscanf(params, "%*s%lf", &PeS); //specify if brownrot == 0
    fscanf(params, "%*s%i", &repulsion); //set to integer != 0 to include HS or Morse repulsion
    fscanf(params, "%*s%i", &attraction); //set to integer != 0 to include AO or Morse attraction
    fscanf(params, "%*s%i", &swimattract); //set to integer != 0 for swimmers to also feel attraction
    fscanf(params, "%*s%i", &morse);//set to integer != 0 to use the Morse potential for attraction+repulsion/repulsion (if attraction == 0)
    fscanf(params, "%*s%lf", &Rg); //sets the cutoff of AO or Morse attraction
    fscanf(params, "%*s%lf", &kappa); //Morse screening length (units of 1/a)
    fscanf(params, "%*s%lf", &Um); //sets the well depth of the AO or Morse attraction
    fscanf(params, "%*s%i", &output); //output frequency (units of timestep)
    fscanf(params, "%*s%i", &swimstart); //timestep that active colloids begin to swim
    fscanf(params, "%*s%i", &spstart); //timestep that swim pressure measurement begins
    fscanf(params, "%*s%lf", &mincellL); //length of unit cell in cell list normalized by largest interaction length
    fclose(params);
  }
  if (attraction != 0 && morse == 0)
    repulsion = 1;
  if (morse != 0 && attraction == 0 && repulsion == 0) {
    repulsion = 1; attraction = 0; swimattract = 0;
  }

  FILE *testing = fopen("test_neighbors", "w");
  setbuf(testing, NULL);
  FILE *data;
  if ( ( data = fopen( "data.gel", "r" ) ) == NULL)
    printf ("data.gel file could not be opened\n");
  else {
    fscanf(data, "%*s%*s%*s");
    fscanf(data, "%i%*s", &tot_colloids);
    fscanf(data, "%i%*s%*s", &nsizes);

    d    = create_2d_array<double>(nsizes, 4);
    
    stressPa   = create_1d_array<double>(3);
    stressPp   = create_1d_array<double>(3);
    stressPpa  = create_1d_array<double>(3);
    stressS    = create_1d_array<double>(3);
    stressAa   = create_1d_array<double>(3);
    stressAp   = create_1d_array<double>(3);
    stressApa  = create_1d_array<double>(3);
    tot_stress = create_1d_array<double>(3);
    
    atom_id   = create_1d_array<int>(tot_colloids);
    mol_id    = create_1d_array<int>(tot_colloids);
    atom_type = create_1d_array<int>(tot_colloids);
    
    x         = create_1d_array<double>(tot_colloids);
    y         = create_1d_array<double>(tot_colloids);
    xi        = create_1d_array<double>(tot_colloids);
    yi        = create_1d_array<double>(tot_colloids);
    image     = create_2d_array<int>(tot_colloids, 2);
    imagei    = create_2d_array<int>(tot_colloids, 2);
    Lr        = create_1d_array<double>(tot_colloids);
    theta     = create_1d_array<double>(tot_colloids);
    Fr        = create_2d_array<double>(tot_colloids, 2);
    Fa        = create_2d_array<double>(tot_colloids, 2);
    Fs        = create_2d_array<double>(tot_colloids, 2);
    Fov       = create_2d_array<double>(tot_colloids, 2);
    v         = create_2d_array<double>(tot_colloids, 2);
   
    fscanf(data, "%*lf%lf%*s%*s", &xmax);
    fscanf(data, "%*lf%lf%*s%*s", &ymax);
    fscanf(data, "%*lf%*lf%*s%*s");
    fscanf(data, "%*lf%*lf%*lf%*s%*s%*s");
    fscanf(data, "%*s");
    for (j=0; j<nsizes; ++j) 
      fscanf(data, "%*i%*lf");
    fscanf(data, "%*s%*s");
    maxL = 0;
    for (j=0; j<nsizes; ++j) {
      fscanf(data, "%i%lf", &intdummy, &dummy);
      intdummy -= 1;
      d[intdummy][0] = dummy;
      if (dummy > maxL)
	maxL = dummy;
    }
    aswim = d[0][0]/2; //swimmer radius
    apass = d[1][0]/2; //average passive colloid radius
   
    fscanf(data, "%*s");
    nactive = 0;
    for (j=0; j<tot_colloids; ++j) {
      fscanf(data, "%i%i%i%lf%lf%*i", &intdummy, &intdummy2, &intdummy3, &dummy, &dummy2);
      if (intdummy3 == 1)
	nactive += 1;
      image[j][0] = 0; image[j][1] = 0;
      v[j][0] = 0; v[j][1] = 0;
      atom_id[j] = intdummy;
      mol_id[j] = intdummy2;
      atom_type[j] = intdummy3;
      x[j] = dummy;
      y[j] = dummy2;
      
      while (x[j] < 0) {
	image[j][0] -= 1;
	x[j] += xmax;
      }
      
      while (x[j] > xmax) {
	image[j][0] += 1;
	x[j] -= xmax;
      }
      
      while (y[j] < 0) {
	image[j][1] -= 1;
	y[j] += ymax;
      }
      
      while (y[j] > ymax) {
	image[j][1] += 1;
	y[j] -= ymax;
      } 
    }
  }
  fclose(data);

  zetaR = zeta*aswim*aswim*4/3;
  if (brownrot != 0)
    tauR = zetaR/kT;
  else 
    tauR = aswim*aswim*zeta/kT/PeR/PeS;

  for (j=0; j<nsizes; ++j) {
    d[j][1] = zeta/apass*d[j][0];
    d[j][2] = zetaR/(apass*apass*apass)*d[j][0]*d[j][0]*d[j][0]/8;
    if (brownrot != 0)
      d[j][3] = d[j][2]/kT;
    else
      d[j][3] = tauR;
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
    
    dummy = tauR;
    
    if (aswim/U < tauR) {
      dummy = aswim/U;
      if (aswim*aswim*zeta/kT < aswim/U)
	dummy = aswim*aswim*zeta/kT;
    }
    else {
      if (aswim*aswim*zeta/kT < tauR)
	dummy = aswim*aswim*zeta/kT;
    }
    dt = dt*dummy;
    steps = stepstR*tauR/dt;
    
    FILE *swimproperties;
    swimproperties = fopen("swimproperties", "w");
    fprintf(swimproperties, "tauR:\t%f\nLR:\t%f\nU:\t%f\nPe_s:\t%f\nPe_R:\t%f\n", tauR, LR, U, PeS, PeR);
    fprintf(swimproperties, "ksTs:\t%f\nna:\t%f\nphia:\t%f\nPi_s:\t%f\n", ksTs, na, phia, Pi_s);
  }
  else {
    dt = dt*(apass*apass*zeta/kT);
    steps = stepstR*(apass*apass*zeta/kT)/dt;
  }

  if (attraction != 0) {
    ddep = Rg*apass*2;
    maxL = maxL + ddep;
  }
  maxL = mincellL*maxL;

  nx = xmax/maxL;
  ny = ymax/maxL;
  dnx = xmax/nx;
  dny = ymax/ny;
  neigh      = create_3d_array<int>(nx, ny, tot_colloids+1);
  neightest  = create_2d_array<int>(nx, ny);

  for (j=0; j<nx; ++j) {
    for (k=0; k<ny; ++k) {
      neigh[j][k][0] = 0;
    }
  }

  for (j=0; j<tot_colloids; ++j) {
    intdummy  = x[j]/dnx;
    intdummy2 = y[j]/dny;
    count = neigh[intdummy][intdummy2][0];
    neigh[intdummy][intdummy2][count+1] = j;
    neigh[intdummy][intdummy2][0] += 1;
    if (atom_type[i] == 1) 
      theta[i] = fRand(0, 2*PI);
  }

/*********************SIMULATION***********************/
 
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
  //fprintf(trajec, "ITEM: ATOMS id mol type x y z ix iy iz\n");
  fprintf(trajec, "ITEM: ATOMS id mol type x y z\n");
  for (ii=0; ii<tot_colloids; ++ii){
    //fprintf(trajec, "%i %i %i %f %f 0 %i %i 0\n", ii+1, mol_id[ii], atom_type[ii], x[ii], y[ii], image[ii][0], image[ii][1]);
    fprintf(trajec, "%i %i %i %f %f 0\n", ii+1, mol_id[ii], atom_type[ii], x[ii], y[ii]);
  }
  
  FILE *swimstress, *totalstress;
  swimstress  = fopen("swimstress.output", "w");
  totalstress = fopen("totalstress.output", "w");
  setbuf(swimstress, NULL);
  setbuf(totalstress, NULL);
  fprintf(swimstress, "Step\tPixx/Pi_s\tPiyy/Pi_s\tPi/Pi_s\n");
  fprintf(totalstress, "Step\tPkt\tPswim\tPHS\tPpot\tPtot\n");
  for (i=0; i<3; ++i) {
    stressS[i] = 0; stressPa[i] = 0; stressPp[i] = 0; stressPpa[i] = 0; stressAa[i] = 0; stressAp[i] = 0; stressApa[i] = 0; 
  }

  for (j=0; j<steps; ++j) {
    if ((j+1) % output == 0) {
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
    }

    //Compute Potential Forces
    if (attraction != 0 || morse != 0) {
      for (i=0; i<nx; ++i) {
    	for (k=0; k<ny; ++k) {
    	  for (l=1; l<=neigh[i][k][0]; ++l) {
    	    for (ii=i; ii<=(i+1); ++ii) {
	      for (kk=(k-1); kk<=(k+1); ++kk) {
		if (ii != i || kk != (k-1)) {
		  if (ii==nx) {
		    idx2 = 0;
		    ax2 = 1;
		  }	
		  else {
		    idx2 = ii;
		    ax2 = 0;
		  }
		  if (kk==ny) {
		    idy2 = 0;
		    ay2 = 1;
		  }	
		  else if (kk == -1)  {
		    idy2 = nx - 1;
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
		  for (ll=limit; ll<=neigh[idx2][idy2][0]; ++ll) {
		    cid1 = neigh[i][k][l];
		    cid2 = neigh[idx2][idy2][ll];
		    delx = (x[cid2] + ax2*xmax) - (x[cid1]);
		    dely = (y[cid2] + ay2*ymax) - (y[cid1]);
		    rsq = delx*delx + dely*dely;
		    r = sqrt(rsq);
		    
		    type1 = atom_type[cid2] - 1;
		    da = d[type1][0];
		    type2 = atom_type[cid2] - 1;
		    db = d[type2][0];
		    D = (da + db)/2;
		    rsq = delx*delx + dely*dely;
		    r = sqrt(rsq);

		    if ((type1 == 0 || type2 == 0) && swimattract == 0)
		      cutoff = 0;
		    else
		      cutoff = D + ddep;
		      
		    if (r < cutoff && r > D && morse == 0) {		      
		      dummy  = -1*delx/r*2*Um/ddep/ddep*(D + ddep - r);
		      dummy2 = -1*dely/r*2*Um/ddep/ddep*(D + ddep - r);

		      Fa[cid2][0] += dummy;
		      Fa[cid2][1] += dummy2;
		      Fa[cid1][0] -= dummy;
		      Fa[cid1][1] -= dummy2;

		      dummy = dummy/xmax/ymax*delx;
		      dummy2 = dummy2/xmax/ymax*dely;
		      
		      stressApa[0] += dummy;
		      stressApa[1] += dummy2;
		      stressApa[2] += (dummy/2 + dummy2/2);
		    }
		    if (morse != 0) {
		      if (attraction == 0) 
			cutoff = D;
		      if ((type1 == 0 || type2 == 0) && swimattract == 0)
			cutoff = D;
		      
		      if (r < cutoff) {
			dummy  = -1*delx/r*2*kappa*Um*(exp(-1*kappa*(r-D)) - exp(-2*kappa*(r-D)));
			dummy2 = -1*dely/r*2*kappa*Um*(exp(-1*kappa*(r-D)) - exp(-2*kappa*(r-D)));
			//			dummy = 0; dummy2 = 0;
			Fa[cid2][0] += dummy;
			Fa[cid2][1] += dummy2;
			Fa[cid1][0] -= dummy;
			Fa[cid1][1] -= dummy2;
			
			dummy = dummy/xmax/ymax*delx;
			dummy2 = dummy2/xmax/ymax*dely;
			
			stressApa[0] += dummy;
			stressApa[1] += dummy2;
			stressApa[2] += (dummy/2 + dummy2/2);
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
    
    for (i=0; i<tot_colloids; ++i) {
      type1 = atom_type[i] - 1;
      Fr[i][0] = sqrt(24*kT*d[type1][1]/dt)*(fRand(0,1)-0.5);
      Fr[i][1] = sqrt(24*kT*d[type1][1]/dt)*(fRand(0,1)-0.5);

      //Compute Swim Force
      if (atom_type[i] == 1 && j >= swimstart) {
	if (j - swimstart == spstart) {
	  xi[i] = x[i];
	  yi[i] = y[i];
	  imagei[i][0] = image[i][0];
	  imagei[i][1] = image[i][1];
	}
	Lr[i] = sqrt(24*d[type1][2]*d[type1][2]/d[type1][3]/dt)*(fRand(0,1)-0.5);
      	theta[i] = theta[i] + Lr[i]/d[type1][2]*dt;
      	Fs[i][0] = d[type1][1]*U*cos(theta[i]);
      	Fs[i][1] = d[type1][1]*U*sin(theta[i]);
      }
      else {
       	Fs[i][0] = 0;
       	Fs[i][1] = 0;
      }
      
      v[i][0] = 1/d[type1][1]*(Fr[i][0] + Fs[i][0] + Fa[i][0]);
      v[i][1] = 1/d[type1][1]*(Fr[i][1] + Fs[i][1] + Fa[i][1]);

      y[i] = y[i] + dt*v[i][1];
      x[i] = x[i] + dt*v[i][0];
    }
    ////O(N^2) Overlap Resolution Testing
    // check = 1;
    // while (check != 0) {
    //   check = 0;
    //   for (l=0; l<tot_colloids; ++l) {
    // 	for (ll=(l+1); ll<tot_colloids; ++ll) {

    // 	  delx = x[l] - x[ll];
    // 	  if (delx <  -1*xmax/2)
    // 	    delx += xmax;
    // 	  if (delx > xmax/2)
    // 	    delx -= xmax;
	  
    // 	  dely = y[l] - y[ll];
    // 	  if (dely < -1*xmax/2) 
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
    // 	    check += 1;
    // 	  }
    // 	}
    //   }
    // }
    
    //Resolve Disk Overlaps with Cell List O(N) using Potential Free Algorithm
    check = 1;
    while (check != 0 && repulsion != 0 && morse == 0) {
      check = 0;
      for (i=0; i<nx; ++i) {
    	for (k=0; k<ny; ++k) {
    	  for (l=1; l<=neigh[i][k][0]; ++l) {
    	    for (ii=i; ii<=(i+1); ++ii) {
	      for (kk=(k-1); kk<=(k+1); ++kk) {
		if (ii != i || kk != (k-1)) {
		  if (ii==nx) {
		    idx2 = 0;
		    ax2 = 1;
		  }	
		  else {
		    idx2 = ii;
		    ax2 = 0;
		  }
		  if (kk==ny) {
		    idy2 = 0;
		    ay2 = 1;
		  }	
		  else if (kk == -1)  {
		    idy2 = nx - 1;
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
		  for (ll=limit; ll<=neigh[idx2][idy2][0]; ++ll) {
		    cid1 = neigh[i][k][l];
		    cid2 = neigh[idx2][idy2][ll];
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
		      dummy  = (delx*delx2/dt*d[type1][1] + delx*delx2/dt*d[type2][1])/2/xmax/ymax;
		      dummy2 = (dely*dely2/dt*d[type1][1] + dely*dely2/dt*d[type1][1])/2/xmax/ymax;

		      //     fprintf(testing, "%f\t%f\n", r, sqrt(delx*delx + dely*dely));

		      stressPpa[0] += dummy;
		      stressPpa[1] += dummy2;
		      stressPpa[2] += (dummy/2 + dummy2/2);
		      check += 1;
		    }
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
    for (l=0; l<tot_colloids; ++l) {
      for (ll=(l+1); ll<tot_colloids; ++ll) {
    	delx = x[l] - x[ll];
    	if (delx <  -1*xmax/2)
    	  delx += xmax;
    	if (delx > xmax/2)
    	  delx -= xmax;
	
    	dely = y[l] - y[ll];
    	if (dely < -1*xmax/2) 
    	  dely += ymax;
    	if (dely > ymax/2)
    	  dely -= ymax;
    	rsq = delx*delx + dely*dely;
    	r = sqrt(rsq);
		  
    	type1 = atom_type[l] - 1;
    	da = d[type1][0];
    	type2 = atom_type[ll] - 1;
    	db = d[type2][0];
    	tolerance = overlap*(da + db)/2;
    	if (r < (da/2 + db/2 - tolerance)) { 		
    	  check += 1;
    	  dummy += r;
    	}
      }
    }
    if (check != 0) 
      dummy = dummy/check;
    else
      dummy = 0;
    
    fprintf(testing, "%i\t%i\t%f\n", j+1, check, dummy);
 
    for (i=0; i<nx; ++i) {
      for (k=0; k<ny; ++k) {
	intdummy = neigh[i][k][0];
	for (l=1; l<=intdummy; ++l)
	  neigh[i][k][l] = 0;
	neigh[i][k][0] = 0;
      }
    }
    
    for (i=0; i<tot_colloids; ++i) {
      if (y[i] < 0) {
	image[i][1] -= 1;
	y[i] += ymax;
      }

      if (y[i] > ymax) {
	image[i][1] += 1;
	y[i] -= ymax;
      }
       
      if (x[i] < 0) {
	image[i][0] -= 1;
	x[i] += xmax;
      }
      
      if (x[i] > xmax) {
	image[i][0] += 1;
	x[i] -= xmax;
      }

      Fa[i][0] = 0; Fa[i][1] = 0;

      if ((j - swimstart) >=  spstart) {
	dummy  = Fs[i][0]*(x[i] - xi[i] + xmax*image[i][0] - xmax*imagei[i][0])/xmax/ymax;
	dummy2 = Fs[i][1]*(y[i] - yi[i] + ymax*image[i][1] - ymax*imagei[i][1])/xmax/ymax;
	stressS[0] += dummy;
	stressS[1] += dummy2;
	stressS[2] += (dummy + dummy2)/2;
      }
      intdummy  = x[i]/dnx;
      intdummy2 = y[i]/dny;
      count = neigh[intdummy][intdummy2][0];
      neigh[intdummy][intdummy2][count+1] = i;
      neigh[intdummy][intdummy2][0] += 1;
      
    }

    //***********************OUTPUT*************************//
    if ((j+1) % output == 0) {
      for (i=0; i < tot_colloids; ++i) {
	//fprintf(trajec, "%i %i %i %f %f 0 %i %i 0\n", i+1, mol_id[i], atom_type[i], x[i], y[i], image[i][0], image[i][1]);
	fprintf(trajec, "%i %i %i %f %f 0\n", i+1, mol_id[i], atom_type[i], x[i], y[i]);
      }
      
      fprintf(swimstress, "%i\t", (j+1));
      fprintf(totalstress, "%i\t", (j+1));      
      fprintf(totalstress, "%f\t%f\t%f\t%f\t%f\n", (tot_colloids/xmax/ymax*kT), (stressS[2]/output), (stressPpa[2]/output), (stressApa[2]/output), (tot_colloids/xmax/ymax*kT + (stressS[2]+stressPpa[2]+stressApa[2])/output));
      for (i=0; i<3; ++i) {
	fprintf(swimstress, "%f\t", (stressS[i]/output/Pi_s));
	stressS[i] = 0; stressPpa[i] = 0; stressApa[i] = 0;
      }
      fprintf(swimstress, "\n");

    }
  }

  duration = ( std::clock() - start )/(double) CLOCKS_PER_SEC;
  
  FILE *t;
  t = fopen("time.output", "w");
  fprintf(t, "Steps:\t%i\n", steps);
  fprintf(t, "Timestep:\t%f\n", dt);
  fprintf(t, "Output:\t%i\n", output);
  fprintf(t, "Duration:\t%f\n", duration);
  return 0;
}
