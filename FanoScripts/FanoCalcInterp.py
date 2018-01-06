###############################################################################
#
# Allocate energy though interpolation of a probability distribution and then
# confirm that we can reproduce the total Fano factor.
#
# Start from BJPB's FanoCalc_bjpj then add P (missing quanta not properly
# accounted for) with different approximations.
#                                                  psihas@fnal.gov
###############################################################################

# Running with atom's Hydrogen requires cell 'markers'.
# Please do not delete cell comments # <codecell>

import os
import numpy
import matplotlib
import pylab
import scipy.interpolate as interp
#%matplotlib inline
matplotlib.rcParams.update({'font.size': 14})

print ("imported modules..")

#.........................
# Load data
#.........................

# Check file path  # <codecell>
#
FileName = "SF6Modes.csv"
FilePath = "FanoData/"

# Load file # <codecell>
#
Data=numpy.loadtxt(FilePath+FileName,delimiter=',')
# Energy per mode
Es=Data[:,0]
# Relative prob per mode
Ns=Data[:,1]
# Is it observable energy?
IsObs=Data[:,2]

# List of mode number indices # <codecell>
#
ModeNumsPrime=numpy.arange(0,len(Ns)+1)
ModeNums=numpy.arange(0,len(Ns))

# Normalized, cumulative, and lookup
NormNs=Ns/sum(Ns)
CumtvN=numpy.concatenate([[0],numpy.cumsum(NormNs)])
LookupFunction=interp.interp1d(CumtvN,ModeNumsPrime,kind='linear')


# Plot cumulative distribution # <codecell>
#
pylab.figure(figsize=(5,5))
pylab.grid()
pylab.plot(CumtvN,linewidth=2)
pylab.xlabel("Mode number")
pylab.ylabel("Cumulative sum P")


#.........................
# One Run Through
#.........................
# Generate one event at the Q value and spend energy into random modes.
# <codecell>

# TODO move these vars upstairs
EnergyToSpend=3e6
EThresh=15

ELeft=EnergyToSpend
ExcitationsSpent=numpy.zeros_like(Ns)

count=0
MaxCount=1e9
while((ELeft>EThresh) and (count<MaxCount)):
    count=count+1
    ModeNum=int(LookupFunction(numpy.random.rand()))
    EThisMode=Es[ModeNum]
    ELeft=ELeft-EThisMode
    ExcitationsSpent[ModeNum]+=1
#    print ("Allocated energy, except for " + str(ELeft))
# TODO: Make progress bar. Printouts make my laptop cry

# Make plots. # <codecell>
pylab.figure(figsize=(5,5))
pylab.semilogy(ModeNums,ExcitationsSpent,'o',color='blue')
pylab.semilogy(ModeNums,ExcitationsSpent*IsObs,'o',color='red')
pylab.xlabel("Mode Num")
pylab.ylabel("N Spent")
pylab.grid()
pylab.show()
pylab.figure(figsize=(5,5))
pylab.semilogy(ModeNums,ExcitationsSpent*Es,'o',color='blue')
pylab.semilogy(ModeNums,ExcitationsSpent*Es*IsObs,'o',color='red')
pylab.xlabel("Mode Num")
pylab.ylabel("E Spent / eV")
pylab.grid()
pylab.show()

#.........................
# With Trials
#.........................
# Generate many low energy events to study event-to-event fluctuations
# <codecell>

NEvents=100
EventEnergy=100000
VisibleEnergy=[]
VisibleQuanta=[]


for i in range(0,NEvents):
    ExcitationsSpent=numpy.zeros_like(Ns)
    ELeft=EventEnergy
    while((ELeft>EThresh)):
        ModeNum=int(LookupFunction(numpy.random.rand()))
        ELeft=ELeft-Es[ModeNum]
        ExcitationsSpent[ModeNum]+=1
    VisibleEnergy.append(sum(ExcitationsSpent*IsObs*Es))
    VisibleQuanta.append(sum(ExcitationsSpent*IsObs))
    if(i%10==0):
        print(i,)
        #TODO add progress bar

# Print plots for Visible N and Visible E #TODO legends
# <codecell>
pylab.hist(VisibleQuanta)
# <codecell>
pylab.hist(VisibleEnergy)

# Calculate resolution
# <codecell>
ResE=numpy.std(VisibleEnergy)/numpy.average(VisibleEnergy)
print("Resn on visible energy",ResE)

ResQ=numpy.std(VisibleQuanta)/numpy.average(VisibleQuanta)
print("Resn on quanta",ResQ)

binomial=1./numpy.sqrt(numpy.average(VisibleQuanta))
print("Resn on quanta from binomial",binomial)

# <codecell>
F1=ResQ/binomial
F2=ResE/binomial

print("Quanta Fano", F1**2)
print("Energy Fano", F2**2)
