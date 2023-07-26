void find_isco(long double z1, long double& isco)
{
    
    int i, j, casenum=1, stop=0, count=0;
    long double detol = 1.0e-8;
    long double rll,rul,rinit=z1,rnew,rold,rstep=1.e-5;
    long double mn[4][4],dmn[4][4];
    long double mnp[4][4],mnm[4][4],dmnp[4][4],dmnm[4][4];
    long double ep,em,eold,enew,omp,omm,omold,omnew;
    long double epsq,emsq;
    long double deold,denew;
    long double dr=1.0e-5;
    long double sqspin=spin*spin;
    
    if(spin>0.)
        rll = 1.+sqrt(1.-sqspin);
    else if(spin<0.)
        rll = 1.+sqrt(1.-sqspin);
    else
        rll = 1.;
    
    if(a13>-5.7){
        rul = rinit; rold = rul;
        metric(rold+dr,Pi/2.,mnp);
        metric(rold-dr,Pi/2.,mnm);
        metric_rderivatives(rold+dr,Pi/2.,dmnp);
        metric_rderivatives(rold-dr,Pi/2.,dmnm);
        omp = (-dmnp[0][3]+sqrt(dmnp[0][3]*dmnp[0][3]-dmnp[0][0]*dmnp[3][3]))/dmnp[3][3];
        omm = (-dmnm[0][3]+sqrt(dmnm[0][3]*dmnm[0][3]-dmnm[0][0]*dmnm[3][3]))/dmnm[3][3];
        ep = - (mnp[0][0]+mnp[0][3]*omp)/sqrt(-mnp[0][0]-2.*mnp[0][3]*omp-mnp[3][3]*omp*omp);
        em = - (mnm[0][0]+mnm[0][3]*omm)/sqrt(-mnm[0][0]-2.*mnm[0][3]*omm-mnm[3][3]*omm*omm);
        deold = 0.5*(ep-em)/dr;
        
        do{
            count++;
            if(count>40){
                printf("No convergence after %i iterations.\n",count);
                break;
            }
            
            rnew = (rll+rul)/2.;
            metric(rnew+dr,Pi/2.,mnp);
            metric(rnew-dr,Pi/2.,mnm);
            metric_rderivatives(rnew+dr,Pi/2.,dmnp);
            metric_rderivatives(rnew-dr,Pi/2.,dmnm);
            omp = (-dmnp[0][3]+sqrt(dmnp[0][3]*dmnp[0][3]-dmnp[0][0]*dmnp[3][3]))/dmnp[3][3];
            omm = (-dmnm[0][3]+sqrt(dmnm[0][3]*dmnm[0][3]-dmnm[0][0]*dmnm[3][3]))/dmnm[3][3];
            ep = - (mnp[0][0]+mnp[0][3]*omp)/sqrt(-mnp[0][0]-2.*mnp[0][3]*omp-mnp[3][3]*omp*omp);
            em = - (mnm[0][0]+mnm[0][3]*omm)/sqrt(-mnm[0][0]-2.*mnm[0][3]*omm-mnm[3][3]*omm*omm);
            denew = 0.5*(ep-em)/dr;
            
            if(fabs(denew)<fabs(detol)) {
                //printf("denew = %Le, deold = %Le\n",denew,deold);
                stop = 1;
            }
            else if((denew*deold)>0.0) {
                if(rnew<rold)
                    rul = rnew;
                else if(rnew>rold)
                    rll = rnew;
                else
                    printf("rold=rnew? rold = %Le, rnew = %Le",rold,rnew);
            }
            else if((denew*deold)<0.0) {
                if(rnew<rold)
                    rll = rnew;
                else if(rnew>rold)
                    rul = rnew;
                else
                    printf("rold=rnew? rold = %Le, rnew = %Le",rold,rnew);
            }
            else
                printf("Compare enew and eold. eold = %Le, enew = %Le\n",deold,denew);
            
            
            rold = rnew;
            omold = omnew;
            deold = denew;
            
        }while(stop==0);
    }
    else {
        rold = 4.0; // ISCO radius for these special cases is never above this, so we start here. -SN
        
        do{
            rnew = rold-rstep;
            
            if(rnew<rll){
                printf("Couldn't find ISCO? rnew = %Le, rll = %Le\n",rnew,rll);
                break;
            }
            
            metric(rnew+dr,Pi/2.,mnp);
            metric(rnew-dr,Pi/2.,mnm);
            metric_rderivatives(rnew+dr,Pi/2.,dmnp);
            metric_rderivatives(rnew-dr,Pi/2.,dmnm);
            omp = (-dmnp[0][3]+sqrt(dmnp[0][3]*dmnp[0][3]-dmnp[0][0]*dmnp[3][3]))/dmnp[3][3];
            omm = (-dmnm[0][3]+sqrt(dmnm[0][3]*dmnm[0][3]-dmnm[0][0]*dmnm[3][3]))/dmnm[3][3];
            epsq = -mnp[0][0]-2.*mnp[0][3]*omp-mnp[3][3]*omp*omp;
            emsq = -mnm[0][0]-2.*mnm[0][3]*omm-mnm[3][3]*omm*omm;
            
            if(epsq>0. && emsq>0.){
                ep = - (mnp[0][0]+mnp[0][3]*omp)/sqrt(epsq);
                em = - (mnm[0][0]+mnm[0][3]*omm)/sqrt(emsq);
                denew = 0.5*(ep-em)/dr;
                
                if(fabs(denew)<fabs(detol)) {
                    //printf("denew = %Le, deold = %Le\n",denew,deold);
                    stop = 1;
                    //break;
                }
                else {
                    //printf("%Le\n",rnew);
                    rold = rnew;
                }
            }
            else{
                printf("epsq = %Le, emsq = %Le, rnew = %Le\n",epsq,emsq,rnew);
                stop = 1;
            }
            
        }while(stop==0);
    }
    
    
    if(stop==1)
        isco=rnew;
}
