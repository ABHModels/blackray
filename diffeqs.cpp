void diffeqs(long double b, long double vars[], long double diffs[])
{
	long double r, th;
	long double kt, kphi;
	long double kt2, kr2, kth2, kp2, ktp, krth;
	long double ch_rtt, ch_rtp, ch_rrr, ch_rrth, ch_rthth, ch_rpp;
	long double ch_thtt, ch_thtp, ch_thrr, ch_thrth, ch_ththth, ch_thpp;
	long double g_tt, g_tp, g_pp, gurr, guthth;
	long double dgttdr, dgttdth, dgtpdr, dgtpdth, dgrrdr, dgrrdth, dgththdr, dgththdth, dgppdr, dgppdth;
	long double hgurr, hguthth;
	long double denom;
	long double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47, t48;
	
	r = vars[0];
	th = vars[1];
	
	t1 = cos(th);
	t2 = r * r;
	t3 = pow(t2, 0.2e1);
	t4 = t2 * t3;
	t5 = r * t3;
	t6 = r * t2;
	t7 = spin * spin;
	t8 = t7 * pow(t1, 0.2e1);
	t9 = (t8 + t2) * r + epsi3;
	t10 = a22 + t2;
	t11 = sin(th);
	t12 = t7 + t2;
	t13 = pow(t10, 0.2e1);
	t14 = pow(t11, 0.2e1);
	t15 = t7 * t13 * t14;
	t17 = a13 + t6;
	t18 = t7 * r;
	t19 = t18 * t10 * t14;
	t20 = t12 * t17;
	t21 = t20 - t19;
	t22 = 2.0 * r;
	t23 = -t22 + t12;
	t17 = pow(t12, 0.2e1) * pow(t17, 0.2e1) - t7 * t4 * t23 * t14;
	t24 = a22 + t7;
	t25 = r * a22;
	t26 = t7 * a13;
	t4 = 2.0 * t4 + t2 * (a13 * t24 + (a22 * t7 + (a13 + t25) * r) * r) + t26 * a22;
	t27 = a52 + t2;
	t28 = t23 * t3 - t15;
	t29 = t8 / 0.3e1 + t2;
	t30 = 0.2e1 / 0.3e1 * t7;
	t31 = 0.3e1 / 0.5e1 * t7;
	t32 = (t31 + t2) * r + 0.2e1 / 0.5e1 * a13;
	t31 = r * t32 - t31 * (a22 / 0.3e1 + t2) * t14;
	t33 = 0.1e1 / t21;
	t35 = pow(t33, 0.2e1);
	t36 = t33 * t35;
	t37 = 0.3e1 * t29;
	t38 = t35 * t14;
	t39 = t7 * a52;
	t40 = a52 + t7;
	t42 = 0.1e1 / r;
	t43 = t9 * t42;
	t44 = -2.0 * t11 * t35 * t1 * (-t43 * t17 + (t5 * t9 * t23 + t17) * t14 * t7) + 4.0 * t1 * t7 * t10 * t11 * t14 * t9 * t17 * t36;
	t21 = 0.1e1 / t21;
	t45 = 0.1e1 / t9;
	t46 = 0.1e1 / t27;
	t47 = 0.1e1 / t23;
	t48 = pow(t21, 0.2e1);
	t19 = 2.0 * t4 * t1 * spin * t11 * (t18 * t14 * (t19 + a22 * epsi3 - ((-(a22 - t7) * r + a13 - epsi3) * r - t8 * t10) * r - t26) + t20 * t9) * t21 * t48;
	t21 = 2.0 * t7 * t1;
	t5 = r * t9 * (-t12 * t3 + 2.0 * t5 + t15) * t48;
	t6 = (-2.0 * (t3 * (-t40 + r) + t8 * ((t2 * (r - 0.1e1) + a52) * r - t39)) * r + t3 * (a52 * (-6.0) - 0.3e1 * epsi3) + (-t2 * t40 + t39) * epsi3 + 4.0 * (epsi3 + t39) * t6) * pow(t47, 0.2e1) * pow(t46, 0.2e1);
	
	g_tt = t5;
	g_pp = t43 * t17 * t14 * t48;
	g_tp = -t4 * t9 * spin * t14 * t48;
	
	gurr = t23 * t27 * t45 * t42;
	guthth = r * t45;
	
	dgttdr = t35 * (t28 * (t9 * (0.10e2 * r * t31 * t33 - 0.1e1) - t37 * r) + (-6.0) * (t2 * ((-0.5e1 / 0.3e1 + r) * r + t30) - t30 * t10 * t14) * t2 * t9);
	dgtpdr = t38 * spin * ((0.20e2 * t9 * t31 * t33 + t29 * (-6.0)) * t4 / 0.2e1 - 0.12e2 * r * ((t24 / 0.6e1 + t2 / 0.3e1) * a13 + t3 + t25 * (0.5e1 / 0.12e2 * t2 + t7 / 0.4e1)) * t9);
	dgrrdr = t6;
	dgththdr = -epsi3 * pow(t42, 0.2e1) + t22;
	dgppdr = 0.10e2 * t38 * t9 * (-t31 * t17 * t33 * t42 + t20 * t32 - 0.4e1 / 0.5e1 * t3 * ((-0.7e1 / 0.4e1 + r) * r + 0.3e1 / 0.4e1 * t7) * t7 * t14) + t38 * t17 * t42 * (-t43 + t37);
	
	dgttdth = 2.0 * t18 * t1 * t11 * t35 * (r * t28 + t13 * t9) - 4.0 * t2 * t1 * t7 * t10 * t28 * t11 * t9 * t36;
	dgtpdth = -t19;
	dgrrdth = -t21 * t2 * t11 * t47 * t46;
	dgththdth = -t21 * t11;
	dgppdth = t44;
	
	hgurr = 0.5*gurr;
	hguthth = 0.5*guthth;
	
	ch_rtt = -hgurr*dgttdr;
	ch_rtp = -hgurr*dgtpdr;
	ch_rrr = hgurr*dgrrdr;
	ch_rrth = hgurr*dgrrdth;
	ch_rthth = -hgurr*dgththdr;
	ch_rpp = -hgurr*dgppdr;
	ch_thtt = -hguthth*dgttdth;
	ch_thtp = -hguthth*dgtpdth;
	ch_thrr = -hguthth*dgrrdth;
	ch_thrth = hguthth*dgththdr;
	ch_ththth = hguthth*dgththdth;
	ch_thpp = -hguthth*dgppdth;
	
	denom = (g_tt*g_pp-g_tp*g_tp);
	
	kt = -(g_pp+b*g_tp)/denom;
	kphi = (g_tp+b*g_tt)/denom;

	diffs[0] = vars[3];
	diffs[1] = vars[4];
	diffs[2] = kphi;
	
	kt2 = kt*kt;
	kr2 = vars[3]*vars[3];
	kth2 = vars[4]*vars[4];
	kp2 = kphi*kphi;
	ktp = kt*kphi;
	krth = vars[3]*vars[4];
	
	diffs[3] = -(ch_rtt*kt2+ch_rrr*kr2+ch_rthth*kth2+ch_rpp*kp2+2.0*(ch_rtp*ktp+ch_rrth*krth));
	diffs[4] = -(ch_thtt*kt2+ch_thrr*kr2+ch_ththth*kth2+ch_thpp*kp2+2.0*(ch_thtp*ktp+ch_thrth*krth));
}
