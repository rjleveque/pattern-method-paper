
""" 
Read yearly predicted.npy gauge data and create the pdf and cummulative 
distributions and plot these in png files. These will be for
a specific pattern, and the output files will have only one
column of probabilities.

Usage: 

    >>> from pattern_yearly import pattern_yearly

    where for example

      INPUT:
          fnames='predicted.npy'  (ex, yearly 1 min tide-data)

          waves are modelled as square

          waves_adjust = [del0,del1,del2,del3] say for 4 waves
                    would be amount to add to each wave to get up
                    to the level of the biggest.  If the third
                    wave has the biggest pattern height, del2=0, e.g.

          waves_indices = [(0,60),(80,140),(160,220),(240,300)] would be 4
                    1 hr waves with a time gap of 20 minutes between each.

          nbins = number of bins for the pdf and cummulative distributions

          hmin,hmax = [-.75,2.75] for the tide heights at CC gauge in MLLW

          Note: Will be converting the CC data in terms of MSL
                hmin-1.13,hmax-1.13 = [-1.88,1.62] tide heights at
                CC in MSL.  The highest tide seen in the year's worth of
                data was 2.63 MLLW which is 1.50 MSL
                toMSL=-1.13
      
          commentsout: file where job printout is written
  
          thisrun_directory: name of the tsunami being processed
 
          A is the amplitude of the tsunami being considered at Gauge 101. It
          is used to compute the PDF and cumulative dist. for the G-method for
          plotting comparison purposes.  If no A is set, it's taken as None.

          Gplot set to True requests the G-method plots are included on the
          plots as comparisons. If not set, it's taken as None.

      def pattern_yearly(fnames,nbins,hmin,hmax,\
                       waves_indices,waves_adjust,\
                       gaugeno,pdffile,cfile,pdfpng,cpng,\
                       commentsout,thisrun_directory,A=None,Gplot=None)

      where for example

      OUTPUT: 
          pdffile='pdf_pattern_testing.txt'        (where heights are ref to MSL)

          cfile='cumulative_pattern_testing.txt'   (where heights are ref to MSL)

          pdfpng='pdf_pattern_testing.png'         (where heights are ref to MSL)

          cpng='cumulative_pattern_testing.png'    (where heights are ref to MSL)

          commentsout  (passed from calling program, where to do the printing)
See:   
"""

import matplotlib
matplotlib.rc('text', usetex = True)
import pylab
from pylab import savefig
from numpy import exp,pi,nanmax,nanmin,load,linspace, nan, isnan, zeros
from numpy import array, shape, sqrt
from scipy.special import erf

def makebins(predicted,pdfhist,chist,hmin,hmax,nbins,waves_indices,\
             waves_adjust,commentsout):
    hvals=linspace(hmin,hmax,nbins+1)
    binsize=(hmax-hmin)/float(nbins)

    #length of the tide data file
    length_predicted=len(predicted) 

    no_waves=len(waves_adjust)
    commentsout.write('\n')
    commentsout.write('The number of waves found was: %s \n\n' %no_waves)
    adjusted_maxes=zeros(no_waves)

    interval_length=waves_indices[-1][1] - waves_indices[0][0] + 1
    commentsout.write('The interval_length used in makebins is: ')
    commentsout.write('%s \n\n' %interval_length)
    commentsout.write('The binsize is: %s \n\n' %binsize)

    # The last 120 values will need to be zapped when the pattern
    # interval_length is 121, for example, since there there will 
    # not be enough lines in the file to make a window of length 121.
    # The interval_length=121 is when the pattern lasts for 2 hours.

    valid_length_predicted=length_predicted - (interval_length -1)  
    valid_length_predicted_start=valid_length_predicted
    for i in range(valid_length_predicted_start):                        
        testpre=predicted[i]
        if (isnan(testpre)):
            valid_length_predicted=valid_length_predicted-1
            #print 'blank: i,predicted[i]: ',i,testpre
        else:  
            # i is a valid pattern starting point
            patternstart=i
            for iwave in range(no_waves):
                wavestart=patternstart+waves_indices[iwave][0]                
                waveend=patternstart+waves_indices[iwave][1]  #to include
                adjusted_maxes[iwave]=nanmax(predicted[wavestart:waveend+1]-\
                                             waves_adjust[iwave])
            testpre_max=nanmax(adjusted_maxes)
            
            # find the correct bin to increment for testpre_max
            itest=int((testpre_max - hmin)/binsize)
            if ((testpre_max > hvals[itest+1])\
               and (testpre_max <= hvals[itest+2])):
               binput=itest+1
               pdfhist[binput]   +=1
               chist[0:binput+1] +=1
            elif ((testpre_max > hvals[itest])\
               and (testpre_max <= hvals[itest+1])):
               binput=itest
               pdfhist[binput]   +=1
               chist[0:binput+1] +=1
            elif ((testpre_max > hvals[itest-1])\
               and (testpre_max <= hvals[itest])):
               binput=itest-1
               pdfhist[binput]   +=1
               chist[0:binput+1] +=1
            else:
               commentsout.write('\n')
               commentsout.write('bin determination failed \n\n')
               commentsout.write('predicted value is %s \n\n ' %testpre_max)
               commentsout.write('itest is %s: \n\n ' %itest)

    return pdfhist,chist,valid_length_predicted
            
def pattern_yearly(fnames,nbins,hmin,hmax,toMSL,\
                     waves_indices,waves_adjust,\
                     gaugeno,pdffile,cfile,pdfpng,cpng,\
                     commentsout,thisrun_directory,A=None,Gplot=None):

    hvals=linspace(hmin,hmax,nbins+1)         # endpoint values for the bins
    pdfhist=zeros(nbins)                      # space for the pdf histogram
    chist=zeros(nbins)                        # space for cumulative histogram 

    #predicted is an array of the predicted tide every 60 sec. in [hmin, hmax)
    #We have a particular pattern of a tsunami described by the waves_indices
    #and waves_adjust lists as described earlier. We pass this pattern over the
    #tide data to determine the water level that can be exceeded and build
    #the following histograms and ultimately pdf and cumulative distributions.
    #
    #pdfhist is a histogram of the number of tidal value heights above MLLW in 
    #each bin  and chist is a cumulative histogram of tidal value heights above
    #above MLLW in each bin 
 
    predicted=load(fnames)            #relative to MLLW so far
    max_predicted=nanmax(predicted)
    min_predicted=nanmin(predicted)
    commentsout.write('\n')
    commentsout.write('max and min of predicted values ref to MLLW are ')
    commentsout.write('%s and %s \n\n' %(max_predicted,min_predicted))
    commentsout.write('max and min of predicted values ref to MSL are ')
    commentsout.write('%s and %s \n\n ' %(max_predicted+toMSL,min_predicted+toMSL))

    pdfhist,chist,valid_length_predicted=makebins(predicted,pdfhist,\
                    chist,hmin,hmax,nbins,waves_indices,\
                    waves_adjust,commentsout)

    #valid_length_predicted is length of colm of probabilities for this pattern

    totalnumber = valid_length_predicted  
    commentsout.write('totalnumber of non-nan data is: %s \n\n' %totalnumber)

    #Turn the histograms above into a pdf and a cumulative distribution
    #by dividing by the total number of valid data entries.
    #Plot the pdf and cumulative distribution
    #Write the pdf and cumulative distribution to a file
 
    #pdfhist and chist are arrays of one colm and totalnumber is a scalar

    pdf=pdfhist/totalnumber 
    cpdf=chist/totalnumber
    commentsout.write('sum of pdfhist is: %s \n\n' %sum(pdfhist))
    commentsout.write('sum of chist is: %s \n\n' %sum(chist))

    #double check reasonableness of pdf and cpdf
    bin_midpoints=(hvals[0:-1]+hvals[1:])/2.0
    meanp=sum(bin_midpoints*pdf[:])
    Eofxsqrd=sum((bin_midpoints)*(bin_midpoints)*pdf[:])    
    sigmap=sqrt(Eofxsqrd-meanp*meanp) 
    commentsout.write('sum of pdf[:] is: %s \n\n' %sum(pdf[:]))
    commentsout.write('totalnumber is: %s \n\n' %totalnumber)
    commentsout.write('chist[0] is (matches totalnumber): %s \n\n' %chist[0])
    commentsout.write('mean(MLLW) and standard deviation of pdf are ')
    commentsout.write('%s and %s \n\n' %(meanp,sigmap))
    commentsout.write('mean(MSL) and standard deviation of pdf are ')
    commentsout.write('%s and %s \n\n' %(meanp+toMSL,sigmap))

    ############ Convert the probability mass function above to real pdf ###
    ############ height of mass function must be divided by the binsize  ###
    #####        so integration of the pdf is one.                       ###
    #
    binsize=(hmax-hmin)/float(nbins)
    pdf[:]=pdf[:]/binsize
    ####
  
    ####  Calculate the pdf as the negative of the derivative of the cpdf ##
    #     for plotting purposes
    #
    #   Find derivatives at the center of the bins
    #   Recall the cpdf values are at the bin left endpoints, and the
    #   right endpoint of the last bin is always guaranteed to have no
    #   value exceeding it, so the probability is 0 at this absent endpoint
    #   Use 2.5 binsizes on each side of the bin-midpoint for smooth curve
    #
    deriv=zeros(nbins);                       #deriv[0],...deriv[nbins-1]
    deriv[0]=(cpdf[1]-cpdf[0])/binsize        #deriv at first bin center
    deriv[1]=(cpdf[2]-cpdf[1])/binsize        #deriv at second bin center
    deriv[2]=(cpdf[3]-cpdf[2])/binsize        #deriv at third bin center
    deriv[-1]=(0.0-cpdf[-1])/binsize          #deriv at nbin-1 bin center
    deriv[-2]=(cpdf[-1]-cpdf[-2])/binsize     #deriv at nbin-2 bin center
    deriv[-3]=(cpdf[-2]-cpdf[-1])/binsize     #deriv at nbin-3 bin center
    for i in range(3,nbins-3):                #deriv for bins 4,..,nbin-4 centers
        deriv[i]=(cpdf[i+3]-cpdf[i-2])/(5*binsize);
    deriv=-deriv
    ####################

    ##### Calculate a normal pdf and associated cumulative distribution
    ##### using the same mean and standard deviation as our pattern pdf 
    ##### Call these pdf_normal and cum_erf. They could be plotted. Not
    ##### being plotted or used at the moment.
    #####
    w0=meanp+toMSL
    w=hvals[:-1]+toMSL
    wmid=bin_midpoints+toMSL
    const=1.0/(sqrt(2*pi)*sigmap)
    pdf_normal=const*exp(-(wmid-w0)*(wmid-w0)/(2*sigmap*sigmap))
    erfarg=(w-w0)/(sqrt(2)*sigmap)
    cum_erf=.5*(1.0-erf(erfarg))
    ####
    ####

    if (A != None):
        #The A value is the amplitude for tsuanmi of interest. Use it to
        #calculate the pdf and cumulative distribution for the G-method
        #A is passed to this routine as a parameter. Below are typical
        #values that A may be:

        #The following 8 values are specific to Crescent City
        C=1.044; Cprime=.707; alpha=.17; beta=.858; alphaprime=.056;
        betaprime=1.119; sigma0=.638; MHHW=.97;

        w0_MOF=C*MHHW*exp(-alpha*(A/sigma0)**beta)
        sigmap_MOF=sigma0*(1.0-Cprime*exp(-alphaprime*(A/sigma0)**betaprime))
        commentsout.write('A used was: %s \n\n' %A)
        commentsout.write('w0_MOF mean of G-method tide cumulative ')
        commentsout.write('was: %s \n\n' %w0_MOF)
        commentsout.write('sigmap_MOF std of G-method tide cumulative ')
        commentsout.write('was: %s \n\n' %sigmap_MOF)
        const_MOF=1.0/(sqrt(2*pi)*sigmap_MOF)   
        pdf_MOF=const_MOF*exp(-(wmid-w0_MOF)*(wmid-w0_MOF)\
                                /(2*sigmap_MOF*sigmap_MOF))
        erfarg_MOF=(w-w0_MOF)/(sqrt(2)*sigmap_MOF)
        cum_erf_MOF=.5*(1.0-erf(erfarg_MOF))

    #Plot the pdf at Gauge 101 for the thisrun_directory tsunami
    figno=0
    fig = pylab.figure(figno+1)
    patterngauge=101

    commentsout.write('Gplot was: %s \n\n' %Gplot)
    if ((Gplot==None) | (A==None)):
        pylab.plot(wmid,deriv[:],'r-')          
        pylab.legend(['Pattern'],'upper left')
        pylab.title('Gauge %s:  Probability Density Function for %s ' \
                           %(patterngauge,thisrun_directory))
        pylab.xlabel('tidal stage above MSL')

    if ((Gplot==True) & (A!=None)):
        pylab.plot(wmid,deriv[:],'r-')          
        pylab.plot(wmid,pdf_MOF[:],'g--',linewidth=2)
        pylab.legend(['Pattern','G-Method'], \
                      'upper left')
        pylab.title('Gauge %s:  Probability Density Functions for %s ' \
                           %(patterngauge,thisrun_directory))
        pylab.xlabel('tidal stage above MSL')
    savefig(pdfpng)

    #Plot the cumulative dist. at Gauge 101 for thisrun_directory tsunami
    fig = pylab.figure(figno+2)
    pylab.axis([-2.0,2.0,0.0,1.1])

    if ((Gplot==None) | (A==None)):
        pylab.plot(w,cpdf[:],'k-')      
        pylab.legend(['Pattern'],'lower left')
        pylab.title('Gauge %s: Cumulative Probability of Exceedance \
                     for %s ' %(patterngauge,thisrun_directory))
        pylab.xlabel('tidal stage above MSL')
        pylab.ylabel('probability')

    if ((Gplot==True) & (A!=None)):
        pylab.plot(w,cpdf[:],'k-')      
        pylab.plot(w,cum_erf_MOF[:],'g--',linewidth=2)
        pylab.legend(['Pattern','G-Method'], 'lower left')
        pylab.title('Gauge %s: Cumulative Probability of Exceedance \
                     for %s ' %(patterngauge,thisrun_directory))
        pylab.xlabel('tidal stage above MSL')
        pylab.ylabel('probability')
    savefig(cpng)
 
    #Write the pdf and cumulative distribution to a file
    outfile=open(pdffile,'w')
    sh=shape(pdf);           
    for i in range(sh[0]):
        tupleprint=(hvals[i]+toMSL,hvals[i+1]+toMSL,pdf[i])
        outfile.write(len(tupleprint)*'%7.3f ' %tupleprint)
        outfile.write('\n')
    outfile.close()
    outfile=open(cfile,'w')
    for i in range(sh[0]):
        tupleprint=(hvals[i]+toMSL,cpdf[i])
        outfile.write(len(tupleprint)*'%7.3f ' %tupleprint)
        outfile.write('\n')
    outfile.close()
#
if __name__ == '__main__':

################ DATA TO CHECK #####################################

    A=1.5;                            # Amplitude of the tsunami
                                      # Optional, A=None is in code
                                      # Set to override.

    Gplot=True                        # Put G-method graphs on plots
                                      # for comparison. This is optional.
                                      # Gplot=False is in code. Set to
                                      # override.

    nbins=100; hmin=-.75;             # nbins, hmin ref to MLLW
    hmax=2.75;                        # hmax ref to MLLW

    toMSL=-1.13                       # MSL=MLLW-1.13

    gaugeno=9419750                   # CC gauge number
    fnames='predicted.npy'            # tidal yearly record at gaugeno

    waves_indices=[(0,60),(70,130)]   # two 1 hr waves, 10 min gap
    waves_adjust=[0.0, 0.3]           # second wave is .3m less 1st

    #waves_indices=[(0,60),(120,180)] # two equal 1hr waves, 1hr. gap
    #waves_adjust=[0.0,0.0]           # no adjustment

    #waves_indices=[(0,120)]          # one two hour wave
    #waves_adjust=[0.0]               # no adjustment

    pdffile='pdf_pattern_testing.txt'
    cfile='cumulative_pattern_testing.txt'
    pdfpng='pdf_pattern_testing.png'
    cpng='cumulative_pattern_testing.png'
    commentsfile='pattern_yearly.output'
    ################################## END OF DATA SETTING ###########
 
 
    ################################## DATA ECHO #####################
    commentsout=open(commentsfile,'w')
    commentsout.write('\n')
    commentsout.write('minimum tidal stage above MLLW (meters) used ')
    commentsout.write('for bins: %s \n\n' %hmin)
    commentsout.write('maximum tidal stage above MLLW (meters) used ')
    commentsout.write('for bins: %s \n\n' %hmax)
    hminadj=hmin-1.13
    commentsout.write('minimum tidal stage above MSL  (meters) used ')
    commentsout.write('for bins: %s \n\n' %hminadj)
    hmaxadj=hmax-1.13
    commentsout.write('maximum tidal stage above MSL  (meters) used ')
    commentsout.write('for bins: %s \n\n' %hmaxadj)
    commentsout.write('MSL = MLLW + %s \n\n ' %toMSL)
    commentsout.write('number of bins used: %s \n\n ' %nbins)
    commentsout.write('file where predicted gauge data is: \n ')
    commentsout.write('%s \n\n' %fnames)
    commentsout.write('The file for Probability Pattern Density: \n ')
    commentsout.write('%s \n\n' %pdffile)
    commentsout.write('The file for Cumulative Pattern Distribution of  ')
    commentsout.write('Exceedance: \n ')
    commentsout.write('%s \n\n' %cfile)
    commentsout.write('The plot for Probability Pattern Density: \n ')
    commentsout.write('%s \n\n' %pdfpng)
    commentsout.write('The plot for Cumulative Pattern Distribution of ')
    commentsout.write('Exceedance: \n ')
    commentsout.write('%s \n\n' %cpng)
    commentsout.write('waves_indices are: \n ')
    commentsout.write('%s \n\n' %waves_indices)
    commentsout.write('waves_adjust: \n ')
    commentsout.write('%s \n\n' %waves_adjust)
    ############################ END OF DATA ECHO ########################
 
    pattern_yearly(fnames,nbins,hmin,hmax,toMSL,\
                       waves_indices,waves_adjust,\
                       gaugeno,pdffile,cfile,pdfpng,cpng,\
                       commentsout,thisrun_directory,A=A,Gplot=Gplot)
