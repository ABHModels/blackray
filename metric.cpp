
void metric(long double r, long double th, long double mn[][4])
{
    long double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16;
    long double g_tt, g_rr, g_thth, g_pp, g_tp; /*metric components with lowered indicies*/
    
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
    t8 = pow(t8, 0.2e1);
    t10 = r * t3 + a13;
    t11 = -t2 * r * t7 * t8 + t10 * t9;
    t12 = -0.2e1 * r + t9;
    t13 = a52 + t3;
    t14 = 0.1e1 / r;
    t15 = t2 * a22;
    t16 = 0.1e1 / t12;
    t13 = 0.1e1 / t13;
    t11 = pow(t11, -0.2e1);
    
    g_tt = t6 * (t2 * pow(t7, 0.2e1) * t8 + 0.2e1 * r * t4 - t4 * t9) * r * t11;
    g_rr = t6 * r * t16 * t13;
    g_thth = t14 * epsi3 + t1 + t3;
    g_pp = (pow(t9, 0.2e1) * pow(t10, 0.2e1) - t2 * t5 * t12 * t8) * t6 * t8 * t11 * t14;
    g_tp = -(0.2e1 * t5 + t3 * (a13 * (a22 + t2) + ((r * a22 + a13) * r + t15) * r) + t15 * a13) * spin * t6 * t8 * t11;
    
    mn[0][0] = g_tt;
    mn[0][3] = g_tp;
    mn[1][1] = g_rr;
    mn[2][2] = g_thth;
    mn[3][0] = mn[0][3];
    mn[3][3] = g_pp;
}
