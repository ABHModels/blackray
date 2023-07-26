void metric_rderivatives(long double r, long double th, long double dmn[][4])
{
    long double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23;
    long double dgttdr, dgtpdr, dgppdr;
    
    t1 = cos(th);
    t2 = spin * spin;
    t3 = r * r;
    t4 = pow(t3, 0.2e1);
    t5 = t3 * t4;
    t1 = t2 * pow(t1, 0.2e1);
    t6 = (t1 + t3) * r + epsi3;
    t7 = a22 + t3;
    t8 = sin(th);
    t9 = t3 + t2;
    t10 = -0.2e1 * r + t9;
    t8 = pow(t8, 0.2e1);
    t11 = r * t3 + a13;
    t12 = t9 * t11;
    t13 = -t2 * r * t7 * t8 + t12;
    t1 = t1 / 0.3e1 + t3;
    t14 = 0.2e1 / 0.3e1 * t2;
    t15 = 0.3e1 / 0.5e1 * t2;
    t16 = (t15 + t3) * r + 0.2e1 / 0.5e1 * a13;
    t15 = r * t16 - t15 * (a22 / 0.3e1 + t3) * t8;
    t13 = 0.1e1 / t13;
    t17 = pow(t13, 0.2e1);
    t18 = 0.3e1 * t1;
    t19 = a22 + t2;
    t20 = r * a22;
    t21 = t2 * a22;
    t22 = 0.1e1 / 0.2e1;
    t23 = t17 * t8;
    t9 = -t2 * t5 * t10 * t8 + pow(t9, 0.2e1) * pow(t11, 0.2e1);
    t11 = 0.1e1 / r;
    
    dgttdr = t17 * ((t6 * (0.10e2 * r * t15 * t13 - 0.1e1) - t18 * r) * (-t2 * pow(t7, 0.2e1) * t8 + t10 * t4) - 0.6e1 * t6 * (t3 * ((-0.5e1 / 0.3e1 + r) * r + t14) - t14 * t7 * t8) * t3);
    dgtpdr = t23 * spin * ((0.20e2 * t6 * t15 * t13 - 0.6e1 * t1) * (t22 * (t3 * (a13 * t19 + ((a13 + t20) * r + t21) * r) + t21 * a13) + t5) - 0.12e2 * t6 * ((t19 / 0.6e1 + t3 / 0.3e1) * a13 + t4 + t20 * (0.5e1 / 0.12e2 * t3 + t2 / 0.4e1)) * r);
    dgppdr = 0.10e2 * t23 * t6 * (-t15 * t9 * t11 * t13 + t12 * t16 - 0.4e1 / 0.5e1 * t4 * ((-0.7e1 / 0.4e1 + r) * r + 0.3e1 / 0.4e1 * t2) * t2 * t8) + t23 * t9 * t11 * (-t11 * t6 + t18);
    
    dmn[0][0] = dgttdr;
    dmn[0][3] = dgtpdr;
    dmn[3][0] = dmn[0][3];
    dmn[3][3] = dgppdr;
}
