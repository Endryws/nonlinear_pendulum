import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"

mpl.rc('xtick', labelsize=22) 
mpl.rc('ytick', labelsize=22)

#==========================================================================

w = 5.61

#==========================================================================
def aprend(fil,tl):
    poin=open(fil,'r');
    lol=open('poin_orb.txt','w');
    g=np.loadtxt(poin,'float','#',None,None,8,(0,1),False,0);
    g=list(g);
    r1=0.15;
    r2=0.3;

    ag={}; at={}; p=0;
	
    for ta in range(1,tl):
            p=0;
            n=0;

            while n < (len(g)-ta):
                if ( (g[n + ta][0] - g[n][0] )**2 + (g[n + ta][1] - g[n][1])**2)**(0.5)  <= r1:
                    ag[(p,ta)]=g[n];
                    p=p+1;
                    g.pop(n);
                    n=n-1;

                n=n+1;
    #print ag
    for n in range(1,tl):
      for m in range(0,len(g)):
          po=[0,0];
          z=1;
          for k in range(m+1,len(g)):
                if ag.get((m,n),[0])[0]!=0 and ag.get((k,n),[0])[0]!=0:
                    if (ag[(m,n)][0] - ag[(k,n)][0])**2 + ((ag[(m,n)][1] - ag[(k,n)][1])**2)**0.5 <= r2:
                        po=[ag[(k,n)][0]+po[0],ag[(k,n)][1]+po[1]];
                        del ag[(k,n)]
                        z=z+1
          if po!=[0,0]:
             at[(n,m)]=[po[0]/(z-1),po[1]/(z-1)];
    
    lol.write('x(rad)\ty(rad/s)\tnorbita\r\n');
    
    for ta in range(1,tl):
        for m in range(0,len(at)):
            if at.get((ta,m),[0,0])[0]!=0:
             print '{0}\t{1}\t{2}\r\n'.format(at.get((ta,m),[0,0])[0],at.get((ta,m),[0,0])[1],ta)
             lol.write('{0}\t{1}\t{2}\r\n'.format(at.get((ta,m),[0,0])[0],at.get((ta,m),[0,0])[1],ta))
         
    lol.close()
    poin.close()
    return at

#================================================================

input_name = "poincare_orb.lis"

aprend(input_name,6)
