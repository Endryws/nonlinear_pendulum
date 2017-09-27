# -*- coding: utf-8 -*-
"""
Created on Fri Apr 01 13:55:59 2016

@author: Dimitri
"""

import scipy as sci
import matplotlib.pyplot as plt
from biblioteca_caos_final import rk4
from biblioteca_caos_final import lpnov
from pendulo import pendulo

D=9.5e-2; m=1.47e-2; I=1.738e-4; g=9.81;
w0=sci.sqrt(m*g*D/I);


def control_tdf(func,x0,nstep,nsap,omsap,mode=1,lyapu=0,con1=0,con2=0,r=0,tob=1):
    '''
    "parametro k11 con1" 
    "parametro k22 con2"
    "parâmetro ETDF controle r"
    "xt paramtro traçar deslocamento ou nao"
    '''
    x=[]; xp=[]; pmapx=[]; pmapxp=[]; t=[]; contr=[]; ee=[]; eep=[]; contrde=[]; tpoin=[];
    pst=''; xst=''; Jaco=sci.matrix(sci.identity(len(x0))); ene=0; tst=-1; prb=1;
        
    ti=0.;
    t.append(ti); x.append(x0[0]); xp.append(x0[1]);
   
    """parametros"""
    dt=sci.double(omsap)/sci.double(nstep);
    alfa=0; "reservado para x"
    beta=0;  "reservado para v"

    control=sci.matrix([[0],[0]]); adm=1; delta=1.5e-5; cut=2; cut1=2; xt=1;
    """"""
    contr.append(control[0,0]); contrde.append(control[1,0]);
     
    lin1=sci.matrix([[1/sci.sqrt(2)],[-1/sci.sqrt(2)]]);
    lin2=sci.matrix([[1/sci.sqrt(2)],[1/sci.sqrt(2)]]);
    l1=0; lv1=[]; tl=[];
    lv1.append(l1); tl.append(cut*omsap);    
#    lv2.append(l2); lv2=[]; l2=0; 
    """Print answers"""
    f = open('caos.txt', 'w'); pmap=open('poincaremap.txt','w');

    lasor=open('lastorbit.txt','w'); lynov=open('lyapunov_it.txt','w');
    f.write('t(s*w)\tx\ty\tenec\tdesloc\r\n'); pmap.write('x\ty\r\n');
    st='{0}\t{1}\t{2}\r\n'.format(ti,x0[0],x0[1]); f.write(st);
    """"""
    
    for i in range(0,nsap):
        for l in range(0,nstep):
            if i<cut:
                   x0=rk4(func,dt,x0,ti,alfa,adm,mode,xt);
            if i>=cut:
#                     alfa=x[l+(i-tob)*(nstep)];
                     beta=xp[l+(i-tob)*(nstep)];
                     if i>(tob*prb+cut):
                         if prb<6:
                             prb=prb+1;
                     for como in range (2,prb):
#                         alfa=r**(como-1)*x[l+(i-como*tob)*(nstep)]+alfa;
                         beta=r**sci.double(como-1)*(xp[l+(i-como*tob)*(nstep)])+beta;
#                     alfa=(1-r)*alfa;
                     beta=(1-r)*beta;
                     [x0,lin1,control,J]=lpnov(func,dt,x0,ti,alfa,6,1,xt,beta,lin1,lin2,delta,con1,con2,control[1,0],Jaco);
                     if i>cut1:
                        
                        l1=l1+sci.log2(sci.linalg.norm(lin1/delta));
#                        l2=l2+sci.log2(sci.linalg.norm(lin2/delta));
                        lamb1=l1*8.8783055563101847/(ti+dt-cut*omsap);
#                        lamb2=l2*8.8783055563101847/(ti+dt-cut*omsap); lv2.append(control[0]);
                        lv1.append(lamb1); tl.append(ti+dt);
                        st='{0}\t{1}\r\n'.format(ti,lamb1);
                        lynov.write(st);
                        if i==nsap-tob:
                           Jaco=sci.matrix(J);
    
            ene=sci.real(control[0,0])+ene;
            ti=ti+dt;
                     
            x.append(x0[0]); xp.append(x0[1]); t.append(ti); contr.append(control[0,0]); contrde.append(sci.real(control[1,0]));
            """escrever dados"""
            xst=xst+'{0}\t{1}\t{2}\t{3}\t{4}\r\n'.format(ti,x0[0],x0[1],sci.real(control[0,0]),sci.real(control[1,0]));
            
            """outros controles e objetivos"""
            if 300<i<600:
               con2=0.42; r=0.1;tob=2; '''0.39/0.05'''
            if 600<i<900:
               con2=0.1;r=0.1;tob=3;

                
        if i>cut1:    
            pmapx.append(x0[0]);
            pmapxp.append(x0[1]);
            tpoin.append(ti);
            if i>(cut1+2*tob+1):
              if sci.sqrt((pmapxp[-1-tob]-x0[1])**2+(pmapx[-1-tob]-x0[0])**2)+sci.sqrt((pmapxp[-1-2*tob]-x0[1])**2+(pmapx[-1-2*tob]-x0[0])**2)<0.01 and tst<0:
                tst=ti;
#            print(x0)
            pst=pst+'{0}\t{1}\t{2}\r\n'.format(x0[0],x0[1],ti);
        f.write(xst);
        xst='';
    pmap.write(pst);
    f.close(); pmap.close(); lynov.close();
    """Configure plots"""
    plt.plot(t,x); plt.xlabel('Time (s)'); plt.ylabel('Position (rad)'); 
    plt.show();  
    plt.plot(x,xp); plt.xlabel('Position (rad)'); plt.ylabel('Velocity(rad/s)'); 
    plt.show();

    """Configure plots"""
    plt.plot(t,contr); plt.xlabel('Time s'); plt.ylabel('Control N'); plt.show();   
    plt.plot(t,contrde); plt.xlabel('Time s'); plt.ylabel('Control m'); plt.show();  
    
    plt.plot(tl,lv1);
    plt.xlabel('tempo (s)'); plt.ylabel('coeficiente de lyapunov'); plt.show();
    print(lamb1);
    plt.show();
  
    plt.plot(pmapx,pmapxp,'ko',markersize=2); plt.xlabel('Position (rad)'); plt.ylabel('Velocity(rad/s)'); 
    plt.show();
    
    plt.plot(tpoin,pmapx,'ko',markersize=2); plt.xlabel('Position (rad)'); plt.ylabel('Velocity(rad/s)'); 
    plt.show();
    """"""
    
    """Last orbits verification"""
    ex=[]; ep=[];
    st='x_utlima\ty_ultim\tx_primeira\ty_primeira\r\n';
    lasor.write(st);
    for h in range (0,2*nstep*tob):
        ex.append(x.pop()); ep.append(xp.pop());
        ee.append(x.pop(0)); eep.append(xp.pop(0));
        st='{0}\t{1}\t{2}\t{3}\r\n'.format(ex[h],ep[h],ee[h],eep[h]);
        lasor.write(st);
    lasor.close();
    plt.plot(ex,ep,'r-',ee,eep,'b-'); plt.xlabel('Position (rad)'); plt.ylabel('Velocity(rad/s)'); plt.show();
    """"""
    print(pmapx.pop())
    z=sci.amax(sci.real((sci.linalg.eig(Jaco)[0]))/(omsap*tob));
    return z,Jaco,tst,ene
    
def pend(x,t,alfa=0,adm=1,xt=0,beta=0,lin1=0,lin2=0,con1=0,con2=0,precon=0,jaco=0):
    """parametros fixos"""
    am=sci.double(1.4585); bm=sci.double(500); c=0.3; d=4.8e-2; D=9.5e-2;
    a=1.6e-1; m=1.47e-2; I=1.738e-4; csi=2.368e-5; mu=1.272e-4; g=9.81;
    comp=1; w0=sci.sqrt(m*g*D/I);
    """"""
    if adm==1:
        
        b=1.5e-2; '''6e-2'''
        T=273.15+12.;
        Tm=273.15+9.3; """287.0;"""
        TA=273.15+16.2; """337.0;"""
        omega=8.5;
        
        sma=am*(T-Tm)*(x[0]*d/(2.*comp))-bm*(x[0]*d/(2.*comp))**3+((bm**2)/(4.*am*(TA-Tm)))*(x[0]*d/(2.*comp))**5.;
        y=sci.sqrt(a**2.+b**2.-2.*a*b*sci.cos((omega/w0)*t))-(a-b)-d*x[0]/2.;
        u=am*(T-Tm)*y/comp-bm*((y/comp)**3.)+((bm**2.)/(4.*am*(TA-Tm)))*(y/comp)**5.;
        x=[x[1],(-csi*x[1]/(w0*I)-mu*sci.sign(x[1])/(m*g*D)-sci.sin(x[0])/2.+(c*d/(2.*D*m*g))*(-sma+u))];
        return x    
    elif adm==6:
        b=1.5e-2; '''6e-2'''
        T=273.15+12.;
        Tm=273.15+9.3; """287.0;"""
        TA=273.15+16.2; """337.0;"""
        omega=8.5;
        maxc=1;
#        err=1;
        if abs(con2*(beta-x[1]))>maxc:
             con2=maxc/(beta-x[1]); control=0; control1=0;
        else:
            if xt>=1.5:
             control=sci.absolute(con1*(alfa-x[0])+con2*(beta-x[1])); control1=0;
#             control=sci.absolute((1+err)+(1-err)*sci.cos(omega*t/w0))*(con1*(alfa-x[0])+con2*(beta-x[1]));control1=0;
            if xt<1.5:
                control=sci.absolute(con1*(alfa-x[0])+con2*(beta-x[1]));
                ap=c*am*(T-Tm);bp=c*bm ;cp=c*bm**2/(4*am*(TA-Tm)); N=(con1*(alfa-x[0])+con2*(beta-x[1]))*(2.*D*m*g)/d;
                cf1=ap+3*(x[0]*d/2)**2*bp+5*(x[0]*d/2)**4*cp;cf2=3*(x[0]*d/2)*bp+10*(x[0]*d/2)**3*cp;cf3=bp+10*(x[0]*d/2)**2*cp;cf4=5*cp*(x[0]*d/2);cf5=cp;
                rot=sci.roots([cf5,cf4,cf3,cf2,cf1,-N]); control1=rot[sci.argmin(abs(rot-precon))]; 
                
#                if control>0.05 :
#                    control=0.05;
#                control1=abs(cf5/6*(control-precon)**6+cf4/5*(control-precon)**5+cf3/4*(control-precon)**4+cf2/3*(control-precon)**3+cf1/2*(control-precon)**2+(ap*(x[0]*d/2)+bp*(x[0]*d/2)**3+cp*(x[0]*d/2)**5)*(control-precon))
#                if control1<0 or type(control)=='numpy.complex128':
#                   control1=0;
               
        sma=am*(T-Tm)*(x[0]*d/(2.*comp))-bm*(x[0]*d/(2.*comp))**3+((bm**2)/(4.*am*(TA-Tm)))*(x[0]*d/(2.*comp))**5.;
        y=sci.sqrt(a**2.+b**2.-2.*a*b*sci.cos((omega/w0)*t))-(a-b)-d*x[0]/2.;
        u=am*(T-Tm)*y/comp-bm*((y/comp)**3.)+((bm**2.)/(4.*am*(TA-Tm)))*(y/comp)**5.;
              
        smal=(d/2)*(am*(T-Tm)-3*bm*(x[0]**2)*(d**2)/4+((5*bm**2)/(4*am*(TA-Tm)))*(x[0]**4)*(d/2)**4);
        ul=(-d/2)*(am*(T-Tm)-3*bm*(y**2)+((5*bm**2)/(4*am*(TA-Tm)))*(y**4));
        
        jac=sci.matrix([[0, 1],[-sci.cos(x[0])/2+(c*d/(2.*D*m*g))*(ul-smal)-con1,-csi/(I*w0)-con2]]);        
        
#        x=[x[1],-csi*x[1]/(w0*I)-mu*sci.sign(x[1])/(m*g*D)-sci.sin(x[0])/2.+(c*d/(2.*D*m*g))*(-sma+u)+(1-(err)*sci.cos(omega*t/w0))*(con1*(alfa-x[0])+con2*(beta-x[1]))];
        x=[x[1],-csi*x[1]/(w0*I)-mu*sci.sign(x[1])/(m*g*D)-sci.sin(x[0])/2.+(c*d/(2.*D*m*g))*(-sma+u)+(con1*(alfa-x[0])+con2*(beta-x[1]))];
        
        lin1=sci.dot(jac,lin1);
        jaco=jac;
        control=sci.matrix([[control1],[control]])
        return x,lin1,control,jaco
