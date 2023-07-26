void intersection(long double x_1, long double y_1, long double z_1, long double x_2, long double y_2, long double z_2, long double x_d[])
{
	long double r1, h1, r2, h2;
	long double pp;
	long double r3, r4, c1;
	long double xmax, xmin;
	
	
	r1 = x_1*sin(y_1);
	h1 = x_1*cos(y_1);
	
	r2 = x_2*sin(y_2);
	h2 = - x_2*cos(y_2);
	
	pp = cos(z_2 - z_1);
	
	r3 = sqrt(r1*r1 + r2*r2 - 2*r1*r2*pp);
	
	r4 = h1*r3/(h1 + h2);
	
	c1 = (r1 - r2*pp)/r3;
	
	x_d[0] = 0.;
	x_d[1] = sqrt(r1*r1 + r4*r4);
	
	if (x_d[1] > r1 && x_d[1] > r2) x_d[1] = 0.5*(r1 + r2);
	else if (x_d[1] < r1 && x_d[1] < r2) x_d[1] = 0.5*(r1 + r2);
	else //printf("x_d[1] = %Le, r1 = %Le, r2 = %Le\n",x_d[1],r1,r2);
	
	x_d[2] = 0.5*Pi;
	
	x_d[3] = asin(r4*sqrt(1. - c1*c1)/x_d[1]);
	
	if (z_2 > z_1) {
		x_d[3] = z_1 + x_d[3];
	} else {
		x_d[3] = z_1 - x_d[3];
	}
}
