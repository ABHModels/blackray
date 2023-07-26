void raytrace(long double xobs, long double yobs, long double iobs, long double rin, long double disk_length_combined, long double traces[5], int& stop_integration)
{
  long double dobs;
  long double xobs2, yobs2;
  long double atol, rtol;
  long double hstart;
  long double t0, r0, th0, phi0;
  long double kt0, kr0, kth0, kphi0;
  long double r02, s0, s02;
  long double fact1, fact2, fact3;
  long double t, r, th, phi;
  long double kt, kr, kth, kphi;
  long double rau, thau, phiau, kthau;
  long double kyem;
  long double const0, const1;
  long double v1, v2;
  long double h, hnext;
  long double Delta;
  long double spin2 = spin*spin;
	
  long double carter, cosem, c02;
  long double b;
	
  long double met[4][4];
  long double diffs[5], vars[5], vars_temp[5], vars_4th[5], vars_5th[5], k1[5], k2[5], k3[5], k4[5], k5[5], k6[5];
  long double xem[4];
  long double gfactor;
  long double err, errmin, errmax;
    
  int check, check2=0;
  int i;
	
  int div;
	
  /* ----- Set computational parameters ----- */
  dobs = 1.0e+8;    /* distance of the observer */
  errmin = 1.0e-8;
  errmax = 1.0e-6;
  atol = 1.0e-10;
  rtol = 1.0e-10;
  long double thtol = 1.0e-8;
  int count, iter;
        
  hstart = -1.0;
    
  /* ----- compute photon initial conditions ----- */
  xobs2 = xobs*xobs;
  yobs2 = yobs*yobs;
  
  fact1 = yobs*sin(iobs) + dobs*cos(iobs);
  fact2 = dobs*sin(iobs) - yobs*cos(iobs);
  
  r02 = xobs2 + yobs2 + dobs*dobs;
  
  r0 = sqrt(r02);
  th0 = acos(fact1/r0);
  phi0 = atan2(xobs,fact2);
			
  s0  = sin(th0);
  s02 = s0*s0;
			
  kr0 = dobs/r0;
  kth0 = -(cos(iobs) - dobs*fact1/r02)/sqrt(r02-fact1*fact1);
  kphi0 = -xobs*sin(iobs)/(xobs2+fact2*fact2);
  
  metric(r0, th0, met);
  
  fact3 = sqrt(met[0][3]*met[0][3]*kphi0*kphi0-met[0][0]*(met[1][1]*kr0*kr0+met[2][2]*kth0*kth0+met[3][3]*kphi0*kphi0));
  
  kt0 = -(met[0][3]*kphi0+fact3)/met[0][0];
  
  b = -(met[3][3]*kphi0+met[0][3]*kt0)/(met[0][0]*kt0+met[0][3]*kphi0);
      
  kr0 /= fact3;
  kth0 /= fact3;
			
  /* ----- carter constant ----- */
			
  c02 = 1. - s02;
			
  carter = yobs2 - spin2*c02 + xobs2*c02;
  carter = sqrt(carter);
			
  /* ----- solve geodesic equations ----- */
			
  r = r0;
  th = th0;
  phi = phi0;
			
  kr = kr0;
  kth = kth0;
			
			
  const0 = kt0;
  const1 = r02*s02*kphi0/kt0;
			
  stop_integration = 0;
			
  h = hstart; count=0; iter=0;
  
  long double a1 = 1.0/4.0;
  long double b1 = 3.0/32.0;
  long double b2 = 9.0/32.0;
  long double c1 = 1932.0/2197.0;
  long double c2 = -7200.0/2197.0;
  long double c3 = 7296.0/2197.0;
  long double d1 = 439.0/216.0;
  long double d2 = -8.0;
  long double d3 = 3680.0/513.0;
  long double d4 = -845.0/4104.0;
  long double e1 = -8.0/27.0;
  long double e2 = 2.0;
  long double e3 = -3544.0/2565.0;
  long double e4 = 1859.0/4104.0;
  long double e5 = -11.0/40.0;
  long double f1 = 25.0/216.0;
  long double f2 = 0.0;
  long double f3 = 1408.0/2565.0;
  long double f4 = 2197.0/4104.0;
  long double f5 = -1.0/5.0;
  long double g1 = 16.0/135.0;
  long double g2 = 0.0;
  long double g3 = 6656.0/12825.0;
  long double g4 = 28561.0/56430.0;
  long double g5 = -9.0/50.0;
  long double g6 = 2.0/55.0;
  
  do {
    iter++;	
    vars[0] = r;
    vars[1] = th;
    vars[2] = phi;
    vars[3] = kr;
    vars[4] = kth;
				
    do {
      
      check = 0;
					
      /* ----- compute RK1 ----- */
					
      diffeqs(b, vars, diffs);
      for(i = 0; i <= 4; i++)
      {
        k1[i] = h*diffs[i];
        vars_temp[i] = vars[i] + a1*k1[i];
      }
					
      /* ----- compute RK2 ----- */
					
      diffeqs(b, vars_temp, diffs);
      for(i = 0; i <= 4; i++)
      {
        k2[i] = h*diffs[i];
        vars_temp[i] = vars[i] + b1*k1[i] + b2*k2[i];
      }
					
      /* ----- compute RK3 ----- */
      
      diffeqs(b, vars_temp, diffs);
      for(i = 0; i <= 4; i++)
      {
        k3[i] = h*diffs[i];
        vars_temp[i] = vars[i] + c1*k1[i] + c2*k2[i] + c3*k3[i];
      }
					
      /* ----- compute RK4 ----- */
					
      diffeqs(b, vars_temp, diffs);
      for(i = 0; i <= 4; i++)
      {
        k4[i] = h*diffs[i];
        vars_temp[i] = vars[i] + d1*k1[i] + d2*k2[i] + d3*k3[i] + d4*k4[i];
      }
      
      /* ----- compute RK5 ----- */
					
      diffeqs(b, vars_temp, diffs);
      for(i = 0; i <= 4; i++)
      {
        k5[i] = h*diffs[i];
        vars_temp[i] = vars[i] + e1*k1[i] + e2*k2[i] + e3*k3[i] + e4*k4[i] + e5*k5[i];
      }
      
      /* ----- compute RK6 ----- */
					
      diffeqs(b, vars_temp, diffs);
      for(i = 0; i <= 4; i++)
        k6[i] = h*diffs[i];
					
      /* ----- local error ----- */ 
		
      for(i=0; i<= 4; i++)
      {
        vars_4th[i] = vars[i] + f1*k1[i] + f2*k2[i] + f3*k3[i] + f4*k4[i] + f5*k5[i];
        vars_5th[i] = vars[i] + g1*k1[i] + g2*k2[i] + g3*k3[i] + g4*k4[i] + g5*k5[i] + g6*k6[i];
        
        err = fabs((vars_4th[i]-vars_5th[i])/max(vars_4th[i], vars[i]));
        
        if(err > errmax && check2 == 0)
            check = 1;
        else if(err < errmin && check != 1 && check2 == 0)
            check = -1;
      }
      
      if(check == 1)
				h/=2.0;
      else if(check == -1)
				h*=2.0;
          
    } while (check == 1);
				
    /* ----- solutions to the fourth-order RKN method ----- */

    thau = th;
    th = vars_4th[1];
				
    if (cos(th) < 0.0)
    {
      check2=1;
      if(fabs(th-thau)<=thtol)
        count++;
	
      if(count>0)
      {
        rau = r;
        phiau = phi;
                      
        kthau = kth;
                      
        r = vars_4th[0];
        phi = vars_4th[2];
                      
        kr = vars_4th[3];
        kth = vars_4th[4];

        intersection(rau, thau, phiau, r, th, phi, xem);
        
        if (xem[1] > rin && xem[1] < disk_length_combined)
        {  
          long double x1, y1, z1, x2, y2, z2, xyd, zd;
          
          x1 = r*sin(th)*cos(phi); y1 = r*sin(th)*sin(phi); z1 = r*cos(th);
          x2 = rau*sin(thau)*cos(phiau); y2 = rau*sin(thau)*sin(phiau); z2 = rau*cos(thau);
          xyd = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
          zd = fabs(z2 - z1);
          // printf("success\n");
          stop_integration = 1;  break; /* the photon hits the disk */
        }
        else
          stop_integration = 2;   /* the photon misses the disk */
      }
      else
      {
        th = thau;
        h/=2.;
      }
    }
    else
    {
      rau = r;
      phiau = phi;
                      
      kthau = kth;
                      
      r = vars_4th[0];
      phi = vars_4th[2];
                      
      kr = vars_4th[3];
      kth = vars_4th[4];
    }
				
    Delta = r*r - 2.*r + spin2;
				
    if (Delta < 1.e-3) stop_integration = 4; //printf("photon crosses the horizon\n"); /* the photon crosses the horizon */
				
    if (r < 1.) stop_integration = 5; //printf("photon crosses the horizon\n");          /* the photon crosses the horizon */
				
    if (r != r) stop_integration = 6; //printf("numerical problem\n");          /* numerical problems! */
				
    if (r > 1.05*dobs) stop_integration = 7; //printf("photon escaped to infinity\n");   /* the photon escapes to infinity */
    
  } while (stop_integration == 0);
			
  if (stop_integration == 1) { 
    redshift(xem[1], const1, gfactor);
    /*Non Kerr PRD 90, 064002 (2014) Eq. 34*/
    cosem = carter*gfactor/sqrt(xem[1]*xem[1]+epsi3/xem[1]);
  }
  else {
    xem[1]  = 0.0;
    gfactor = 0.0;
  }
    
  traces[0] = xem[1];
  traces[1] = cosem;
  //traces[2] = xem[3];
  traces[3] = gfactor;
}
