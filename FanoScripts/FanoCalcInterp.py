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
import matplotlib.pyplot as plt

#%matplotlib inline
matplotlib.rcParams.update({'font.size': 14})

print ("imported modules..")

def make_histo(mylist,mycol,space,mylabel,myx,mytxt,outname):
    cols=['gray','xkcd:bubblegum pink','black','DarkRed','cornflowerblue','hotpink']
    NP=numpy.array(mylist)
    binN = numpy.linspace(NP[0]*(1-space), NP[0]*(1+space), 20)
    Sig = plt
    Sig.hist(NP, binN, alpha=0.5, label=mylabel, color=cols[mycol])
    Sig.legend(loc='upper right')
    Sig.xlabel(myx)
    Sig.text(numpy.average(NP)*(1-space), 250, mytxt)
    Sig.grid(True)
    Sig.savefig(outname+'.png', dpi=100)
    Sig.show()



#.........................
# Load data
#.........................
# Check file path  # <codecell>
#
ApproxMessyStates=1

FilePath = "FanoData/"
FileName = "SF6Modes_100keV_B.csv"
OutPath = "test/AppB/"
if ApproxMessyStates > 0:
    FileName = "SF6Modes_100keV_A.csv"
    OutPath = "test/AppC/"
OutFile =OutPath+"testlog.log"

f=open(OutFile, "a+")
f.write("\n___ NEW ITERATION___\n")
f.write("\nFrom ",FileName)
f.write("\nApprox ",ApproxMessyStates)

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

f.write("\nIter 1:\nESpent = ",FractionToSpend," * ",EnergyToSpend)
if ( ApproxMessyStates > 0 ):
    f.write("\n >> MinE for shell states = ",EOneQ)
else :
    f.write("\n >> Ommitting shell states ")

ELeft=EnergyToSpend*FractionToSpend
ExcitationsSpent=numpy.zeros_like(Ns)
AddedQuantaSpent=numpy.zeros_like(Ns)

count=0
MaxCount=1e9
while((ELeft>EThresh) and (count<MaxCount)):

    count=count+1
    ModeNum=int(LookupFunction(numpy.random.rand()))
    EThisMode=Es[ModeNum]

    if ( ApproxMessyStates > 0 ) and ( MissQ[ModeNum] > 0 ):
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

NEvents        = 10000
EventEnergy    = 100000
EnergyToAssign = EventEnergy*FractionToSpend
VisibleEnergy  = []
VisibleQuanta  = []
VisibleQuantaCorr = []
ExcitationsSpentPerEvt=[]

f.write("\n\n--- ",NEvents," of E = ",EventEnergy," assigned = ",FractionToSpend," ---\n")


for i in range(0,NEvents):

    ExcitationsSpent=numpy.zeros_like(Ns)
    AddedQuantaSpent=numpy.zeros_like(Ns)
    ELeft=EnergyToAssign
    while((ELeft>EThresh)):

        ModeNum=int(LookupFunction(numpy.random.rand()))
        EThisMode=Es[ModeNum]
        if ( ApproxMessyStates > 0 ) and ( MissQ[ModeNum] > 0 ):
            EThisMode=EOneQ

        ELeft=ELeft-EThisMode
        ExcitationsSpent[ModeNum]+=( 1+XTRAs[ModeNum] )
        AddedQuantaSpent[ModeNum]+=(1+XTRAs[ModeNum]+PsToAdd[ModeNum])

    ExcitationsSpentPerEvt.append(ExcitationsSpent)
    VisibleEnergy.append(sum(ExcitationsSpent*IsObs*Es))
    VisibleQuanta.append(sum(ExcitationsSpent*IsObs))
    VisibleQuantaCorr.append(sum(AddedQuantaSpent*IsObs))

    if( i%10 == 0 ):
        print(i,)
        #TODO add progress bar

print("N_evt = ",VisibleQuanta)
f.write("\n Visible Quanta in each event: \n",VisibleQuanta)# print("N_evt +p2= ",VisibleQuantaCorr)

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
f.write("Avg E = ",numpy.average(VisibleEnergy),"  +- ",numpy.std(VisibleEnergy))


ResE=numpy.std(VisibleEnergy)/numpy.average(VisibleEnergy)
print("Resn on visible energy",ResE)
f.write("Resn on visible energy",ResE)

print ("Avg Q = ",numpy.average(VisibleQuanta),"  +- ",numpy.std(VisibleQuanta))
ResQ=numpy.std(VisibleQuanta)/numpy.average(VisibleQuanta)
print("Resn on quanta",ResQ)
f.write("Resn on quanta",ResQ)

binomial=1./numpy.sqrt(numpy.average(VisibleQuanta))
print("Resn on quanta from binomial",binomial)
f.write("Resn on quanta from binomial",binomial)

F1=ResQ/binomial
F2=ResE/binomial

print("Quanta Fano", F1**2)
print("Energy Fano", F2**2)
f.write("Quanta Fano", F1**2)
f.write("Energy Fano", F2**2)
# <codecell>


txt=r'Res: %.3f'%(ResQ)+'\n'+r'F: %.3f'%(F1**2)
make_histo(VisibleQuanta,1,'Total','Visible Quanta',txt,OutPath+'VisibleQuanta')
#make_histo(mylist,mycol,mylabel,myx,mytxt,outname)
txt=r'Res: %.3f'%(ResE)+'\n'+r'F: %.3f'%(F2**2)
make_histo(VisibleEnergy,1,'Total','Visible Quanta',txt,OutPath+'VisibleEnergy')

#print(ExcitationsSpentPerEvt)
# <codecell>
OutPath = "AppB/"
if ApproxMessyStates > 0:
    FileName = "SF6Modes_100keV_A.csv"
    OutPath = "AppC/"
print("make N arrays")
AllExcitationsSpent=numpy.zeros_like(Ns)
colors=['gray','blue','black','DarkRed','cornflowerblue','hotpink']

for mIdx,m in enumerate(AllExcitationsSpent):
    AllExcitationsN=[]
    NiHistos=[]
    for evt in ExcitationsSpentPerEvt:
        AllExcitationsN.append(evt[mIdx])
        AllExcitationsSpent[mIdx]+=evt[mIdx]
    if IsObs[mIdx]==1:
        print ("Mode ",mIdx," ... Total Qs=",AllExcitationsSpent[mIdx])
        ResN=numpy.std(AllExcitationsN)/numpy.average(AllExcitationsN)
        binomialN=1./numpy.sqrt(numpy.average(AllExcitationsN))
        FN=ResN/binomialN

        AllExcitationsNP=numpy.array(AllExcitationsN)
        binN = numpy.linspace(AllExcitationsNP[0]*0.85, AllExcitationsNP[0]*1.15, 20)
        Sig = plt
        Sig.hist(AllExcitationsNP, binN, alpha=0.5, label=r'Mode  %d'%(mIdx), color=colors[0])
        Sig.legend(loc='upper right')
        Sig.xlabel(r'Quanta in Mode  %d'%(mIdx))
        Sig.text(AllExcitationsNP[0]*0.85, .025, r'Res: %.3f'%(ResN)+'\n'+r'F: %.3f'%(FN**2))
        Sig.grid(True)
        Sig.savefig(r'ModesQ_%d.png' %(mIdx), dpi=100)
        Sig.show()

        figname=(r'ModesQ_%d' %(mIdx))
        lab=(r'Mode  %d'%(mIdx))
        xlab=(r'Quanta in Mode  %d'%(mIdx))
        txtN=(r'Res: %.3f'%(ResN)+'\n'+r'F: %.3f'%(FN**2))
        make_histo(AllExcitationsNP,0,0.25,lab,xlab,txtN,OutPath+figname)
        make_histo(AllExcitationsNP,0,0.35,lab,xlab,txtN,OutPath+"zoom"+figname)

# <codecell>
print ("...")
