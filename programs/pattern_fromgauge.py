###################################################################
#
#   Purpose: Find the pattern for a particular tsunami and output
#            lines to be used for pattern_yearly_run.py.
#
#   Example: CSZBe01r07 tsunami called L1.
#
#   Run:     python pattern_fromgauge.py CSZBe01r07
#
#   input:  User makes selections in the code clearly marked below
#           gauge101_t_eta.txt (tsunami gauge time series)
#
#   output: pattern.txt  job run printing, last two lines are the
#                        pattern (waves_indices and waves_adjust)
#                        to be input to pattern_yearly_run.py
#                        There may be two more lines preceeding
#                        these that give the amplitude A and
#                        the Gplot command to be used with
#                        pattern_yearly_run.py as well.
#
#           pattern.png  plot of the pattern on top of the tsunami
#                        time series in gauge101_t_eta.txt
####################################################################

import pylab
from pylab import *
import numpy
from numpy import *
from matplotlib.mlab import find
import matplotlib.pyplot as plt
import sys

thisrun_directory=sys.argv[1] #name of the tsunami (e.g., CSZBe01r07)

#############   User sets these parameters ##########################
#
patternfile='pattern.txt'     #last two lines give the pattern
patternpng='pattern.png'      #plot of the pattern
#
#               To create a new pattern you will most likely want to
#               leave the 6 parameters below as they are presently 
#               set.  That is iexist=False, idummy=False, 
#               igaugeread=True, icreate=True, ipopfind=True, and
#               ipop_do=True.
#
# 
iexist=False
#iexist=True    #working with exisiting pattern
                #True if existing, False if creating here
                #If True, set icreate=False, igaugeread=False,
                #         idummy=False. Then set ipopfind and
                #         ipop_do as you like, probably both True.
                

idummy=False
#idummy=True    #dummy sin function tsunami test
                #If True, takes the time series as know sine
                #         function, so set iexist=False, 
                #         igaugeread=False, icreate=True. You 
                #         can set ipopfind, and ipop_do as
                #         you like, probably both to True.

#igaugeread=False
igaugeread=True #pattern doesn't exist yet, read tsunami gauge data
                #If True, you must set icreate=True, iexist=False,
                #         idummy=False.  You can set ipopfind and
                #         ipop_do as you like, probably both to True.

#icreate=False
icreate=True    #create the pattern for the gauge data

ipopfind=True   #find the unnecessary waves and print them
ipop_do=True    #get rid of the unnecessary waves

### Crescent City specific parameters
HAT=1.5;        #Highest Tide  (at Crescent City)
LAT=-1.83;      #Lowest Tide   (at Crescent City)
MHW=.77;        #MHW           (at Crescent City)
MHHW=.97;       #MHHW          (at Crescent City)

if(igaugeread==True): #User should check the following
    bg=MHHW    
    gaugefile= 'gauge101_t_eta.txt'

############End of User Input ######################################

outfile=open(patternfile,'w')
outfile.write('thisrun_directory is: %s \n\n'  %thisrun_directory)

if (iexist==True): 
    if 0:
        #USING EXISTING PATTERN BELOW 
        #This was the old info for CSZBe01r07 ################################
        #Note that the first wave arrived at 20 minutes after earthquake######
        #and the waves_indices have not been shifted by 20 yet to start the ##
        #first wave at time 0. This resetting will be done in the code.     ##
        #We print again after reprocessing to see if any waves pop.   ########
        #    
        maxht_aboveMHW=6.5180
        uplift_subside=1.46
        waves_indices=[[20,37],[42,65],[86,110],[118,128],[150,159],[166,180]]
        waves_max=[6.51801801802,2.26126126126,4.07432432432,\
                     2.02477477477,0.763513513514,\
                     1.19707207207]
   
        waves_adjust=[0.0,6.5180-2.26126126126,6.5180-4.07432432432,\
                      6.5180-2.02477477477,6.5180-0.763513513514,\
                      6.5180-1.19707207207]
        bg=MHHW
        bgminus=bg-MHW 
        A=maxht_aboveMHW - bgminus + uplift_subside
        outfile.write('MHW was: %s and bg for job run was %s \n\n' %(MHW,bg))
        outfile.write('The max gauge value above MHW was: ')
        outfile.write('%s \n\n' %maxht_aboveMHW)
        outfile.write('The Uplift-Subsidence was: %s \n\n' %uplift_subside)
        outfile.write('The Tsunami Amplitude was: %s \n\n' %A)
        xlim1=5; xlim2=185;   #for plotting limits on the time axis
        nwaves=len(waves_indices)
    if 1:
        #USING EXISTING PATTERN BELOW 
        #This was the info for AASZe02r01 used in the paper     ##############
        #Note that the first wave arrived at 263 minutes after earthquake#####
        #and the waves_indices have not been shifted by 263 yet to start the #
        #first wave at time 0. This resetting will be done in the code below.#
        #We print again after reprocessing to see if any waves pop.          #
        #We see that no waves popped, since the variation of tide is 3.33 m  #    
        #at Crescent City. -1.83m is lowest, and 1.5m the highest tide and   #
        #the popping is conservative based on these limits.                  #

        maxht_aboveMHW=1.7
        uplift_subside=0.0
        waves_indices=[[263,305],[347,387],[423,465],[506,538],[572,588],\
                       [607,613],[635,659]]
        waves_adjust=[1.7-1.14,1.7-1.20,1.7-1.18,1.7-.92,1.7-.82,\
                      1.7-.25,0.0]
        waves_max=[1.14,1.20,1.18,.92,.82,.25,1.7]
        bg=MHHW
        bgminus=bg-MHW 
        A=maxht_aboveMHW - bgminus + uplift_subside
        outfile.write('MHW was: %s and bg for job run was %s \n\n' %(MHW,bg))
        outfile.write('The max gauge value above MHW was: ')
        outfile.write('%s \n\n' %maxht_aboveMHW)
        outfile.write('The Uplift-Subsidence was: %s \n\n' %uplift_subside)
        outfile.write('The Tsunami Amplitude was: %s \n\n' %A)
        xlim1=260; xlim2=665;           #for plotting limits on the time axis
        nwaves=len(waves_indices)

if (idummy==True):
    #SELECT FOR A DUMMY sin function Tsunami #######
    bg=.77    
    bgminus=bg-MHW;
    value=zeros((361,2))
    value[:,0]=linspace(0,180,361)
    tmin=value[:,0]
    value[:,1]=6.0*exp(-tmin/45.) * sin(4.0*pi*tmin/60.0)
    gauge=value[:,1];
    uplift_subside=0.0;
    pylab.plot(tmin,gauge,'r-')
    pylab.hold(True)
    outfile.write('MHW was: %s and bg for job run was %s \n\n' %(MHW,bg))
    outfile.write('The Uplift-Subsidence was: %s \n\n' %uplift_subside)
    xlim1=tmin[0]-5; xlim2=tmin[-1]+5;
 
if (igaugeread==True):
    bgminus=bg-MHW;
    value=loadtxt(gaugefile,skiprows=0)

    #determine if there was uplift or subsidence
    #by seeing if the value changed from the bg
    #sealevel after one second of time elaspsed.
    #The gauge data is relative to MHW, so check
    #if the value is different than bgminus.
    #Only CSZB earthquake zone will have this so
    #only have to check for these tsunamis

    if (thisrun_directory[3]=='B'):
        #indices where time <= 1sec
        before_indices=find(value[:,0] <= 1.0)

        #first index where time > 1sec
        the_index=before_indices[-1] + 1

        if (value[the_index,1] > bgminus):
            #subtracts from amplitude
            uplift_subside=bgminus-value[the_index,1]

        else:
            #adds to amplitude
            uplift_subside=bgminus-value[the_index,1]
        outfile.write('Uplift(<0) or Subsidence(>0) was: ')
        outfile.write('%s \n\n' %uplift_subside)


    #gauge file will have time in sec, so convert to minutes 
    tmin=value[:,0]/60.
    gauge=value[:,1];
    pylab.plot(tmin,gauge,'r-')
    pylab.hold(True)
    xlim1=tmin[0]-5; xlim2=tmin[-1]+5;
 
if (icreate==True):
    #CREATE THE PATTERN

    #Find the maximum height above MHW, gauges report above MHW
    maxht_aboveMHW=amax(gauge)
    outfile.write('The max gauge value above MHW was: ')
    outfile.write('%s \n\n' %maxht_aboveMHW)

    #Amplitude of the tsunami
    A=maxht_aboveMHW - bgminus + uplift_subside
    outfile.write('Tsunami amplitude was: %s \n\n' %A)

    pos=find(gauge -bgminus > 0.0)

    #if a differ location is 1, contiguous positive
    differ=zeros(len(pos),dtype=int)
    differ[0]=0
    differ[1:]=pos[1:]-pos[0:-1]

    istart=find(differ != 1)
    starter=pos[istart]
    nwaves=len(starter)

    iend=zeros(nwaves,dtype=int)
    iend[0:-1]=istart[1:]-1
    iend[-1]=len(pos)-1
    ender=pos[iend]

    waves_indices=[]; waves_max=[]; waves_adjust=[]; popping=[];
    for k in range(nwaves):
        maxval=amax(gauge[starter[k]:ender[k]+1])
        dval=maxht_aboveMHW-maxval
        waves_max.append(maxval)
        waves_adjust.append(dval)
        wavetimes=[int(round(tmin[starter[k]])),int(round(tmin[ender[k]]))]
        waves_indices.append(wavetimes)

    #check for overlapping wave times (if equal, that is ok)
    #if found, set left end equal to last wave's right endpoint
    for k in range(1,nwaves):
         if(waves_indices[k][0] < waves_indices[k-1][1]):
             waves_indices[k][0]=waves_indices[k-1][1]
             outfile.write(' wave %s left endpoint reset \n\n' %k)

# Plot the Tsunami and the Pattern before starting index set to 0
# for the first wave and before popping unnecessary waves.
t=zeros(4*nwaves); pat=zeros(4*nwaves);
for iwave in range(nwaves):
    starting=iwave*4
    t[starting]=waves_indices[iwave][0]
    t[starting+1]=waves_indices[iwave][0]
    t[starting+2]=waves_indices[iwave][1]
    t[starting+3]=waves_indices[iwave][1]
    pat[starting]=bgminus
    pat[starting+1]=waves_max[iwave]
    pat[starting+2]=waves_max[iwave]
    pat[starting+3]=bgminus
pylab.plot(t,pat,'b-')

# Write out the results so far before popping and before setting starting
# index to 0 for the first wave.
outfile.write('The number of waves before popping is : %s \n\n' %nwaves)
outfile.write('The Pattern graph is: %s \n\n'  %patternpng)
outfile.write('The sealevel used for the Gauge data ')
outfile.write('to make Pattern was: %s \n\n' %bg)
outfile.write('The number of waves found was: %s \n\n' %nwaves)
outfile.write('The highest Gauge data above MHW in each wave: \n')
outfile.write(' %s \n\n' %waves_max)
outfile.write('The waves_indices (before popping and starting at 0) are:')
outfile.write('\n')
outfile.write('waves_indices=')
outfile.write('%s \n\n' %waves_indices)
outfile.write('The waves_adjust before popping is: \n')
outfile.write('waves_adjust=')
outfile.write('%s \n\n' %waves_adjust)

#Find and print the unnecessary waves
if (ipopfind==True):
    poplist=[];
    for k in range(nwaves):
        dval=maxht_aboveMHW-waves_max[k]
        if (dval > (HAT-LAT)):      #add k to poplist
            poplist.append(k)
    outfile.write('poplist is: \n')
    outfile.write('%s \n\n' %poplist)

#Get rid of the unnecessary waves (by popping them off)
#and then plot the new pattern on same graph as the data and
#the pattern before popping.
if (ipop_do == True): 
    npop=len(poplist)
    for k in range(npop-1,-1,-1):
        waves_max.pop(poplist[k])
        waves_indices.pop(poplist[k])
        waves_adjust.pop(poplist[k])
    nwaves=len(waves_indices)

    #Plot the graph using the necessary waves before
    #resetting the index of the waves to start at 0
    #
    t=zeros(4*nwaves); pat=zeros(4*nwaves);
    for iwave in range(nwaves):
        starting=iwave*4
        t[starting]=waves_indices[iwave][0]
        t[starting+1]=waves_indices[iwave][0]
        t[starting+2]=waves_indices[iwave][1]
        t[starting+3]=waves_indices[iwave][1]
        pat[starting]=bgminus
        pat[starting+1]=waves_max[iwave]
        pat[starting+2]=waves_max[iwave]
        pat[starting+3]=bgminus
    pylab.plot(t,pat,'k-',linewidth=3)

pylab.xlabel('Time in minutes after tsunami')
pylab.ylabel('Tsunami Pattern (Black), Tsunami in Red')
pylab.xlim(xlim1,xlim2)
pylab.ylim(-maxht_aboveMHW-.5,maxht_aboveMHW+.5)
titletsunami=' %s ' %thisrun_directory
pylab.title(titletsunami)
plt.savefig(patternpng)

# For the surviving waves after popping (not zapped by HAT-LAT),
# start waves at 0 as required by the program pattern_yearly_run.py
# that creates the cumulative probability distribution.
#
wavestart_minute=waves_indices[0][0]
for j in range(nwaves):
    waves_indices[j][0]=waves_indices[j][0]- wavestart_minute
    waves_indices[j][1]=waves_indices[j][1]- wavestart_minute

outfile.write('FINAL RESULTS AFTER POPPING UNNECESSARY WAVES: ')
outfile.write('\n\n')
if (ipop_do == True):
    outfile.write('Popping unnecessary waves was done: \n\n')
else:
    outfile.write('Popping unnecessary waves not done: \n\n')

outfile.write('The number of waves used was: %s \n\n' %nwaves)
outfile.write('The first wave started at ')
outfile.write('%s minutes after earthquake \n\n' %wavestart_minute)
outfile.write('The highest Gauge data above MHW in each wave: \n')
outfile.write('%s \n\n' %waves_max)
outfile.write('The data for %s ' %thisrun_directory)
outfile.write('to use with pattern_yearly_run.py: A and Gplot optional \n')
outfile.write('A=')
outfile.write('%s \n' %A)
outfile.write('Gplot=True \n')
outfile.write('waves_indices=')
outfile.write('%s \n' %waves_indices)
outfile.write('waves_adjust=')
outfile.write('%s \n ' %waves_adjust)
