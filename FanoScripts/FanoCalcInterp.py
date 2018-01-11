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
ApproxMessyStates=0

FilePath = "FanoData/"
FileName = "SF6Modes_100keV_B.csv"
if ApproxMessyStates > 0:
    FileName = "SF6Modes_100keV_A.csv"

# Load file # <codecell>
#
Data=numpy.loadtxt(FilePath+FileName,delimiter=',')
# Energy per mode
Es=Data[:,0]
# Avg Quanta (used for relative prob per mode)
Ns=Data[:,1]
# Is it observable energy?
IsObs=Data[:,2]
# Missing Quanta
MissQ=Data[:,3]
# Fixed # Quanta to add
XTRAs=Data[:,4]
# Mode names
#NAMEs=Data[:,5]

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



# <codecell>
#............................................
# Approx A never worked but I kept it as a
# side calculation.
#............................................
AvgMissingQuanta=130 # From degrad output

NinKshells=Ns*MissQ
EinKshells=Ns*Es*MissQ
EinKshells = EinKshells/(sum(EinKshells))

print (">>\n Adding ",AvgMissingQuanta," for side calculation. Main calculation not affected ")
# print (EinKshells)
# print("\n or ",AvgMissingQuanta/sum(NinKshells), " quanta per kshell N")

# Weight by N in each mode
PsToAdd=EinKshells*(AvgMissingQuanta/sum(NinKshells))
#............................................


# <codecell>
#.........................
# One Run Through
#.........................
# Generate one event at the Q value and spend energy into random modes.

# TODO move these vars upstairs
EnergyToSpend=3e6
EThresh=15
EOneQ=15.67
FractionToSpend=0.85

ELeft=EnergyToSpend*FractionToSpend
ExcitationsSpent=numpy.zeros_like(Ns)
AddedQuantaSpent=numpy.zeros_like(Ns)

count=0
MaxCount=1e9
while((ELeft>EThresh) and (count<MaxCount)):

    count=count+1
    ModeNum=int(LookupFunction(numpy.random.rand()))
    EThisMode=Es[ModeNum]

    if ( ApproxMessyStates > 0 ) && ( MissQ[ModeNum] > 0 ):
         EThisMode=EOneQ

    ELeft=ELeft-EThisMode
    ExcitationsSpent[ModeNum]+=1
    AddedQuantaSpent[ModeNum]+=PsToAdd[ModeNum]

#    print ("Allocated energy, except for " + str(ELeft))
# TODO: Make progress bar. Printouts make my laptop cry

# Make plots. # <codecell>
pylab.figure(figsize=(5,5))
pylab.semilogy(ModeNums,ExcitationsSpent,'o',color='blue')
pylab.semilogy(ModeNums,ExcitationsSpent*IsObs,'o',color='red')
pylab.semilogy(ModeNums,AddedQuantaSpent*IsObs,'o',color='green')
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

NEvents        = 50000
EventEnergy    = 100000
EnergyToAssign = EventEnergy*FractionToSpend
VisibleEnergy  = []
VisibleQuanta  = []
VisibleQuantaCorr = []

for i in range(0,NEvents):

    ExcitationsSpent=numpy.zeros_like(Ns)
    AddedQuantaSpent=numpy.zeros_like(Ns)
    ELeft=EnergyToAssign
    while((ELeft>EThresh)):

        ModeNum=int(LookupFunction(numpy.random.rand()))
        EThisMode=Es[ModeNum]
        if ( ApproxMessyStates > 0 ) && ( MissQ[ModeNum] > 0 ):
            EThisMode=EOneQ

        ELeft=ELeft-EThisMode
        ExcitationsSpent[ModeNum]+=( 1+XTRAs[ModeNum] )
        AddedQuantaSpent[ModeNum]+=(1+XTRAs[ModeNum]+PsToAdd[ModeNum])

    VisibleEnergy.append(sum(ExcitationsSpent*IsObs*Es))
    VisibleQuanta.append(sum(ExcitationsSpent*IsObs))
    VisibleQuantaCorr.append(sum(AddedQuantaSpent*IsObs))

    if( i%10 == 0 ):
        print(i,)
        #TODO add progress bar

print("N_evt = ",VisibleQuanta)
# print("N_evt +p2= ",VisibleQuantaCorr)

# Print plots for Visible N and Visible E #TODO legends

# <codecell>
pylab.hist(VisibleQuanta)
#<codecell>
# pylab.hist(VisibleQuantaCorr)
# <codecell>
pylab.hist(VisibleEnergy)

# Calculate resolution & Fano factor
# <codecell>
print ("\n>>\n>> For ", NEvents, " events, allocating ",FractionToSpend , "\n>>")
print ("Avg E = ",numpy.average(VisibleEnergy),"  +- ",numpy.std(VisibleEnergy))

ResE=numpy.std(VisibleEnergy)/numpy.average(VisibleEnergy)
print("Resn on visible energy",ResE)

print ("Avg Q = ",numpy.average(VisibleQuanta),"  +- ",numpy.std(VisibleQuanta))
ResQ=numpy.std(VisibleQuanta)/numpy.average(VisibleQuanta)
print("Resn on quanta",ResQ)

binomial=1./numpy.sqrt(numpy.average(VisibleQuanta))
print("Resn on quanta from binomial",binomial)

F1=ResQ/binomial
F2=ResE/binomial

print("Quanta Fano", F1**2)
print("Energy Fano", F2**2)
