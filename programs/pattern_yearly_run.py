#########################################################################
##       Purpose:  Create cumulative distribution for the pattern       
##
##       Run:    python pattern_yearly_run.py CSZBe01r07
##
##       Calls:  function pattern_yearly from pattern_yearly.py
##
##               pattern_yearly(fnames,nbins,hmin,hmax,toMSL,\
##                     waves_indices,waves_adjust,\
##                     gaugeno,pdffile,cfile,pdfpng,cpng,\
##                     commentsout,thisrun_directory,A=A,Gplot=Gplot)
##                     #commentsout,thisrun_directory)
##
##               (The call above to pattern_yearly is the last statement
##                in this program.  Edit it to select whether or not A 
##                and G are being set. If not included in the call, they
##                are used as None in the function pattern_yearly.)
##
##       Input:  fnames='predicted.npy'   #file for tidal yearly record
##
##               User specifies which tsunami, and other data including
##               the pattern in waves_indices and waves_adjust by editing
##               this code as indicated below.
##
##               Optionally, one can specify the amplitude and whether
##               the G-Method plots should be included on the plots for
##               comparison. See the present setup under the CSZBe01r07
##               tsunami in the code pattern_yearly_run.py. The amplitude
##               A is not needed to make the pattern, but it is required
##               to find the pdf and cumulative distributions for the
##               G-Method. These optional parameters are called A and
##               Gplot in the code.  If used, A should be in meters and
##               Gplot=True.
##
##       Output:
##            pdffile='pdf_pattern.txt'                 #the pdf
##            cfile='cumulative_pattern.txt'            #the cumulative dist.
##            pdfpng='pdf_pattern.png'                  #pdf plot
##            cpng='cumulative_pattern.png'             #cumulative plot  
##            commentsfile='pattern_yearly_run.output'  #job run comments
##
##       The sys.argv[1] in the python call above is for
##       the Cascadia tsunami CSZBe01r07 of the paper and was
##       the L1 Bandon source. The name CSZBe01r07 will be used
##       for the plot titles. 
##########################################################################

from  pattern_yearly import pattern_yearly
from numpy import array,exp
import sys

#tsunami name, e.g. CSZBe01r07 read from command line
thisrun_directory=sys.argv[1]            

#########################USER SETS DATA, STARTING HERE ####################
pdffile='pdf_pattern.txt'                 #file for the pdf
cfile='cumulative_pattern.txt'            #file for the cumulative dist.
pdfpng='pdf_pattern.png'                  #file for pdf plot
cpng='cumulative_pattern.png'             #file for cumulative plot  
commentsfile='pattern_yearly_run.output'  #file for job run comments
fnames='predicted.npy'                    #file for tidal yearly record

nbins=70; hmin=-.75; hmax=2.75;          # number of bins 
                                          # hmin,hmax ref to MLLW

toMSL=-1.13                               # MSL=MLLW-1.13
gaugeno=9419750                           # CC gauge number

### Tsunami examples, select with (if 1:) the one you want to try.
### Or create your own here.

#waves info for Toy example 1
if 0:
    waves_indices=[(0,60),(120,180)]      # two waves, 1hr. gap, 3hr total
    waves_adjust=[0.0,0.0]                # waves equal heights
    A=1.5                                 # Amp. of highest wave is 1.5m.
    Gplot=True                            # Plot G-method results as well

#waves info for Toy example 2
if 0:
    waves_indices=[(0,60),(70,130)]       # two 1 hr waves, 10 min. gap
    waves_adjust=[0.0,0.3]                # 2nd wave .3 meters less than 1st
                                          # A and Gplot not specified.

#waves info for Toy example 3 
if 0:
    waves_indices=[(0,60)]                # one 1hr wave
    waves_adjust=[0.0]                    # no adjusting
                                          # A and Gplot not specified

#waves info for Toy example 4
if 0:
    waves_indices=[(0,120)]              # one two hour wave
    waves_adjust=[0.0]                   # no adjustment
                                         # A and Gplot not specified.

# The next three are 3 tsunamis similar to those used in the paper
#
#waves info for CSZBe01r07 ############# a large Bandon source tsunami
if 1:  
    #These values are output from pattern_fromgauge.py and are the last four
    #lines in the file pattern.txt
    A=9.547774 
    Gplot=True 
    waves_indices=[[0, 17]]
    waves_adjust=[0.0]
##############################################################################

#waves info for AASZe02r01 ############ In the paper, illustrating the pattern
if 0:
    waves_indices=[[263,305],[347,387],[423,465],[506,538],[572,588],[607,613],\
                   [635,659]]
    for j in range(len(waves_indices)):
        waves_indices[j][0]=waves_indices[j][0]-263
        waves_indices[j][1]=waves_indices[j][1]-263
    waves_adjust=[
    1.69992795389-1.13886887608,
    1.69992795389-1.2019092219,
    1.69992795389-1.18299711816,
    1.69992795389-0.918227665706,
    1.69992795389-0.823667146974,
    1.69992795389-0.25,0.]
    A=1.5                 #Optional, sets amplitude, water above MHHW
    Gplot=True            #Optional, requests G-method pdf and cumulative 
                          #to be put on graphs as comparison.

#waves info for AASZe03r01, similar to 1964 Alaska event ######################
if 0:
    waves_indices=[[243,266],[277,287],[316,339],[355,364],[382,414],\
                   [422,434],[449,459],[476,510]]
    for j in range(len(waves_indices)):
        waves_indices[j][0]=waves_indices[j][0]-243
        waves_indices[j][1]=waves_indices[j][1]-243
    waves_adjust=[0.0,4.11740890688-1.14170040486,4.11740890688-3.995951417,\
                  4.11740890688-0.95951417004,4.11740890688-2.09311740891,\
                  4.11740890688-1.2024291498,4.11740890688-0.655870445344,\
                  4.11740890688-1.56680161943]
    A=3.91740890688   #water above MHHW, no subsidence
    Gplot=True

##
#Waves info for all the sources used in the CC study available
#upon request.                    
#
################################## END OF USER DATA #######################
 
################################## DATA ECHO ##############################
commentsout=open(commentsfile,'w')
commentsout.write('thisrun_directory is: %s \n\n' %thisrun_directory)
commentsout.write('minimum tidal stage above MLLW (meters) used ')
commentsout.write('for bins: %s \n\n' %hmin)
commentsout.write('maximum tidal stage above MLLW (meters) used ')
commentsout.write('for bins: %s \n\n' %hmax)
hminadj=hmin-1.13
commentsout.write('minimum tidal stage above MSL (meters) used ')
commentsout.write('for bins: %s \n\n' %hminadj)
hmaxadj=hmax-1.13
commentsout.write('maximum tidal stage above MSL (meters) used ')
commentsout.write('for bins: %s \n\n' %hmaxadj)
commentsout.write('MSL = MLLW + %s \n\n' %toMSL)
commentsout.write('number of bins used: %s \n\n' %nbins)
commentsout.write('The file where predicted gauge data is: \n ')
commentsout.write('%s \n\n' %fnames)
commentsout.write('The file for Probability Pattern Density: \n ')
commentsout.write('%s \n\n' %pdffile)
commentsout.write('The file for Cumulative Pattern Distribution of ')
commentsout.write('Exceedance: \n ')
commentsout.write('%s \n\n' %cfile)
commentsout.write('The plot file for Probability Pattern Density: \n ')
commentsout.write('%s \n\n' %pdfpng)
commentsout.write('The plot file for Cumulative Pattern Distribution of ')
commentsout.write('Exceedance: \n ')
commentsout.write('%s \n\n' %cpng)
commentsout.write('waves_indices are: \n ')
commentsout.write('%s \n\n' %waves_indices)
commentsout.write('waves_adjust: \n ')
commentsout.write('%s \n\n' %waves_adjust)
############################ END OF DATA ECHO #############################
 
pattern_yearly(fnames,nbins,hmin,hmax,toMSL,\
                   waves_indices,waves_adjust,\
                   gaugeno,pdffile,cfile,pdfpng,cpng,\
                   commentsout,thisrun_directory,A=A,Gplot=Gplot)
                   #commentsout,thisrun_directory)
