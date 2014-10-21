import matplotlib
matplotlib.rc('text', usetex = True)
from pylab import *
import os

#Use the dt distributions from Crescent City
#tide gauge to make example plots
d = loadtxt('cumulative_probs.yearly.100.txt')

figure(1,(12,9))
clf()
axes((.1,.1,.8,.38))

#Add 1.13 to get s referenced to MSL
s = d[:,0] - 2. +1.13 

plot([-2.6,0.6,0.6,1.0],[1,1,0,0],'k')
for i in [7,6,5,3,1]:
    plot(s,d[:,i])
legend(['dt = infinity','dt = 24 hours', 'dt = 12 hours', 'dt = 4 hours', \
    'dt = 2 hour', 'dt = 0'],'lower left')
title('Cumulative Distributions for dt-Method',fontsize=18)
plot([-1.08,-1.08,-3],[1.05,0.776,0.776],'r--',linewidth=3)
plot([-1.08],[0.776],'ro')
text(-2.8,0.8,'desired probability',fontsize=17)

ylim(-0.1,1.1)
xticks([-2,-1,0,0.6,1.],['MLW','MSL','MHW','HHW',r'${\bf \hat{\xi}}$'],fontsize=18)
yticks(fontsize=18)
ylabel('probability ',fontsize=18)
xlabel(r'tide stage ${\bf \hat{\xi}}$',fontsize=18)

axes((.1,.55,.8,.38))
s = linspace(-3,0.6,101)
def h(s):
    h = (s+1.4) + (s+1.6)**2
    h = where(s<-1.4, 0, h)
    return h
srun = linspace(-2.,0,6)
plot([-1.4],[0.0],'bo')  
plot(srun,h(srun),'bo')
 
#Now plot the black line
srun2=zeros(9)
srun2[0]=-3.; srun2[1:3]=srun[0:2]; srun2[3]=-1.4; srun2[4:8]=srun[2:6];
srun2[8]=.6;
plot(srun2,h(srun2),'k')
 
plot([-3,-1.08,-1.08],[0.59,0.59,-0.5],'r--',linewidth=3)
ylim(-0.5,5)
xticks([-2,-1,0,0.6,1.],['MLW','MSL','MHW','HHW',r'${\bf \hat{\xi}}$'],fontsize=18)
yticks([0,2,4],['0','1','2'],fontsize=18)

plot([-1.08],[0.59],'ro')
text(-2.8,0.8,r'exceedance level ${\bf \zeta_i}$',fontsize=17)
text(-1.2,-.45,r'${\bf \hat{\xi_i}}$',fontsize=15)
xlim(-3,1)
ylabel(r'Quantity of Interest ${\bf \zeta}$',fontsize=18)
title(r"${\bf E_{jk}}$'s GeoClaw Simulation Curve at one Location",fontsize=18)
savefig('tidepofzeta_dt.png')
