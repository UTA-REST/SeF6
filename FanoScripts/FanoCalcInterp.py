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
cols=['xkcd:medium blue','xkcd:pinkish red','xkcd:bubblegum pink','black','DarkRed','cornflowerblue','hotpink']

def make_histo(mylist,mycol,space,mylabel,myx,mytxt,outname):
    NP=numpy.array(mylist)
    binN = numpy.linspace(NP[0]*(1-space), NP[0]*(1+space), 20)
    Sig = plt
    Sig.hist(NP, binN, alpha=0.5, label=mylabel, color=cols[mycol])
    Sig.legend(loc='upper right')
    Sig.xlabel(myx)
    Sig.text(numpy.mean(NP)*(1-space), 250, mytxt)
    Sig.grid(True)
    Sig.savefig(outname+'.pdf', dpi=100)
    Sig.show()



#.........................
# Load data
#.........................
# Check file path  # <codecell>
#
ApproxMessyStates=0

FilePath = "FanoData/"
FileName = "SF6Modes_100keV_B.csv"
OutPath = "test/AppB/"
if ApproxMessyStates > 0:
    FileName = "SF6Modes_100keV_A.csv"
    OutPath = "test/AppC/"
OutFile =OutPath+"testlog.log"

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

# List of mode number indices
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
pylab.savefig(OutPath+"interp.pdf")


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

    if ( ApproxMessyStates > 0 ) and ( MissQ[ModeNum] > 0 ):
         EThisMode=EOneQ

    ELeft=ELeft-EThisMode
    ExcitationsSpent[ModeNum]+=1
    AddedQuantaSpent[ModeNum]+=PsToAdd[ModeNum]

#    print ("Allocated energy, except for " + str(ELeft))
# TODO: Make progress bar. Printouts make my laptop cry

# Make plots. # <codecell>
pylab.figure(figsize=(5,5))
pylab.semilogy(ModeNums,ExcitationsSpent,'o',color='xkcd:metallic blue',label='All Modes')
pylab.semilogy(ModeNums,ExcitationsSpent*IsObs,'o',color='xkcd:darkish red',label='Visible Modes')
pylab.semilogy(ModeNums,AddedQuantaSpent*IsObs,'o',color='xkcd:sunflower',label='Added Q to Modes')
pylab.legend(loc='lower left')
pylab.xlabel("Mode Num")
pylab.ylabel("N Spent")
pylab.grid()
pylab.savefig(OutPath+"TestEvent_Nspent.pdf")
pylab.show()
pylab.figure(figsize=(5,5))
pylab.semilogy(ModeNums,ExcitationsSpent*Es,'o',color='xkcd:metallic blue',label='All Modes')
pylab.semilogy(ModeNums,ExcitationsSpent*Es*IsObs,'o',color='xkcd:darkish red',label='Visible Modes')
pylab.legend(loc='lower left')
pylab.xlabel("Mode Num")
pylab.ylabel("E Spent [eV]")
pylab.grid()
pylab.savefig(OutPath+"TestEvent_Espent.pdf")
pylab.show()

#.........................
# With Trials
#.........................
# Generate many low energy events to study event-to-event fluctuations
# <codecell>

NEvents        = 1000
EventEnergy    = 100000
EnergyToAssign = EventEnergy*FractionToSpend
VisibleEnergy  = []
VisibleQuanta  = []
VisibleQuantaCorr = []
ExcitationsSpentPerEvt=[]


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
# Print plots for Visible N and Visible E #TODO legends
#<codecell>
#pylab.hist(VisibleQuanta)
#<codecell>
# pylab.hist(VisibleQuantaCorr)
#<codecell>
#pylab.hist(VisibleEnergy)

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

print("Quanta Fano ", F1**2)
print("Energy Fano ", F2**2)

# <codecell>

f=open(OutFile, "a+")
f.write("\n____________________\n____________________\n")
f.write("\n___ NEW ITERATION___\n")
f.write("\n____________________\n____________________\n")
f.write("\nFrom \n"+FileName)
f.write("\n____________________\n")
f.write("\nApprox %d\n\n"%ApproxMessyStates)
if ( ApproxMessyStates > 0 ):
    f.write("\n >> MinE for shell states = %.3f"%EOneQ)
else :
    f.write("\n >> Ommitting shell states \n\n")
f.write("\n____________________\n")

f.write("\nIter 1:\nESpent = %.3f * %.3f"%(FractionToSpend,EnergyToSpend))

f.write("\nEvents:\n--- %d of E = %.3f assigned = %.3f  ---\n"%(NEvents,EventEnergy,FractionToSpend))

f.write("Avg E = %.3f +- %.3f\n\n"
%(numpy.average(VisibleEnergy),numpy.std(VisibleEnergy)))
f.write("Avg Q = %.3f +- %.3f\n\n"
%(numpy.average(VisibleQuanta),numpy.std(VisibleQuanta)))
f.write("Resn on visible energy %.5f\n"%ResE)
f.write("Resn on quanta %.5f\n"%ResQ)
f.write("Resn on quanta from binomial%.5f\n"%binomial)
f.write("\n____________________\n")

f.write("Quanta Fano %.5f\n"%(F1**2))
f.write("Energy Fano %.5f\n"%(F2**2))
f.write("\n____________________\n")


# <codecell>
# Print Plots of Event E and Q
#
txt=r'Res: %.3f'%(ResQ)+'\n'+r'F: %.3f'%(F1**2)
make_histo(VisibleQuanta,0,0.05,'Total','Visible Quanta',txt,OutPath+'VisibleQuanta')
#make_histo(mylist,mycol,mylabel,myx,mytxt,outname)
txt=r'Res: %.3f'%(ResE)+'\n'+r'F: %.3f'%(F2**2)
make_histo(VisibleEnergy,0,0.1,'Total','Visible Energy [eV]',txt,OutPath+'VisibleEnergy')

#print(ExcitationsSpentPerEvt)
# <codecell>
print("make N arrays")
AllExcitationsSpent=numpy.zeros_like(Ns)
colors=['gray','blue','black','DarkRed','cornflowerblue','hotpink']

f=open(OutFile, "a+")
f.write("\n____________________\n____________________\n")
f.write("\n\n --- VISIBLE MODES: ---\n\n")
f.write("\n____________________\n____________________\n")

for mIdx,m in enumerate(AllExcitationsSpent):
    AllExcitationsN=[]
    NiHistos=[]
    for evt in ExcitationsSpentPerEvt:
        AllExcitationsN.append(evt[mIdx])
        AllExcitationsSpent[mIdx]+=evt[mIdx]
    if IsObs[mIdx]==1:
        print ("\nMode ",mIdx," ... Total Qs=",AllExcitationsSpent[mIdx])
        f.write("...\n Mode %d ...\n"%mIdx)

        ResN=numpy.std(AllExcitationsN)/numpy.average(AllExcitationsN)
        binomialN=1./numpy.sqrt(numpy.average(AllExcitationsN))
        FN=ResN/binomialN

        AllExcitationsNP=numpy.array(AllExcitationsN)
        binN = numpy.linspace(AllExcitationsNP[0]*0.85, AllExcitationsNP[0]*1.15, 20)
        # Sig = plt
        # Sig.hist(AllExcitationsNP, binN, alpha=0.5, label=r'Mode  %d'%(mIdx), color=colors[0])
        # Sig.legend(loc='upper right')
        # Sig.xlabel(r'Quanta in Mode  %d'%(mIdx))
        # Sig.text(AllExcitationsNP[0]*0.85, .025, r'Res: %.3f'%(ResN)+'\n'+r'F: %.3f'%(FN**2))
        # Sig.grid(True)
        # Sig.savefig(r'ModesQ_%d.png' %(mIdx), dpi=100)
        # Sig.show()

        figname=(r'ModesQ_%d' %(mIdx))
        lab=(r'Mode  %d'%(mIdx))
        xlab=(r'Quanta in Mode  %d'%(mIdx))
        txtN=(r'Res: %.3f'%(ResN)+'\n'+r'F: %.3f'%(FN**2))
        f.write(" ... ... Res: %.3f F: %.3f"%(ResN,FN**2))
        make_histo(AllExcitationsNP,1,0.25,lab,xlab,txtN,OutPath+figname)
        make_histo(AllExcitationsNP,1,0.35,lab,xlab,txtN,OutPath+"zoom"+figname)

# <codecell>
print ("...")
f=open(OutFile, "a+")
f.write("...\n")
