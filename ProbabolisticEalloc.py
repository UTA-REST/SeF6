#######################################################
#  Probabilistic Energy Allocation
#
#  Generate events at a Q value and distribute the
#  energy probabilistically according to the assigned
#  number of quanta form MAGBOLTZ.
#
#                                    psihas@fnal.gov
#######################################################

import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.stats import norm

# Generate Nevents
Nevents = 5
# Allocate Qvalue energy to each
Qvalue = 2995e3
#Qvalue =26934.61795 #test from orig. output

verb = 0 # TODO: learn true/false syntax in python

## The following info is available for each mode:
#    name, shortname, Eloss, avgQuanta, fraction of N I calculated

class modeinfo:
    def __init__(self,name,shortname,Eloss,avgQuanta,fracN):
        self.name          = name
        self.shortname     = shortname
        self.Eloss         = Eloss
        self.avgQuanta     = avgQuanta
        self.avgEnergy     = avgQuanta * Eloss
        self.fracN         = fracN # from excel
        self.probThisModeE = 0.0
        self.probThisModeN = 0.0
#        self.probThisModeN  = avgQuanta / TOTALavgN
        self.probCutHigh   = -1.0
        self.probCutLow    = -1.0
        self.accumulatedN  = 0
        self.accumulatedE  = 0.0
        self.ALLaccumulatedN = list()
        self.ALLaccumulatedE = list()

        def echo(self):
         print "Mode is " + self.name
         print self.Eloss,self.avgQuanta

         def thisN(self):
             print self.accumulatedN

         def thisN(self):
             print self.accumulatedE


def define_modes(globalvar):
    modes = list()
    modes.append(modeinfo( "IONISATION SF5 +", " ", 15.67, 520.6, 0.054713035070572))
    modes.append(modeinfo( "IONISATION SF4 +", " ", 18.5, 36.982, 0.00388666435455224))
    modes.append(modeinfo( "IONISATION SF3 +", " ", 18.8, 125.86, 0.0132273964540572))
    modes.append(modeinfo( "IONISATION SF2 +", " ", 27, 30.134, 0.00316696619058129))
    modes.append(modeinfo( "IONISATION SF +", " ", 31, 58.77, 0.00617649840779393))
    modes.append(modeinfo( "IONISATION SUM OF: S + AND F +", " ", 37, 101.24, 0.0106399302161827))
    modes.append(modeinfo( "IONISATION SUM OF (SF3 SF2 SF)2+", " ", 46.5, 13.628, 0.00143224979243518))
    modes.append(modeinfo( "IONISATION SULFUR L3 SHELL", " ", 164.16, 4.508, 0.000473773265651438))
    modes.append(modeinfo( "IONISATION SULFUR L2 SHELL", " ", 165.36, 2.176, 0.000228689136215069))
    modes.append(modeinfo( "IONISATION SULFUR L1 SHELL", " ", 230.9, 1.012, 0.00010635726371767))
    modes.append(modeinfo( "IONISATION SULFUR K SHELL", " ", 2472, 0.044, 4.62422885728999E-06))
    modes.append(modeinfo( "IONISATION FLUORINE K SHELL", " ", 685.4, 0.298, 3.13186408971004E-05))
    # modes.append(modeinfo( "ATTACHMENT (VALID FOR T=300KELVIN) 30.4040", " ", , , 0))
    modes.append(modeinfo( "VIBRATION V4 SUPERELASTIC", " ", -0.076253, 15.846, 0.00166535296528675))
    modes.append(modeinfo( "VIBRATION V4 ANISOTROPIC", " ", 0.076253, 392.84, 0.0412859560067682))
    modes.append(modeinfo( "VIBRATION V1 SUPERELASTIC", " ", -0.096032, 1.392, 0.000146293785666992))
    modes.append(modeinfo( "VIBRATION V1 ISOTROPIC", " ", 0.096032, 70.832, 0.00744416769135374))
    modes.append(modeinfo( "VIBRATION V3 SUPERELASTIC", " ", -0.11754, 59.136, 0.00621496358419775))
    modes.append(modeinfo( "VIBRATION V3 ANISOTROPIC", " ", 0.11754, 7544.368, 0.792883732173072))
    modes.append(modeinfo( "VIBRATION 2V1", " ", 0.192064, 29.024, 0.00305030950804511))
    modes.append(modeinfo( "VIBRATION 3V1", " ", 0.288096, 14.63, 0.00153755609504892))
    modes.append(modeinfo( "VIBRATION 4V1", " ", 0.384128, 8.476, 0.000890794631690681))
    modes.append(modeinfo( "VIBRATION 5V1 + HIGHER HARMONICS", " ", 0.48016, 7.36, 0.000773507372492144))
    modes.append(modeinfo( "EXC. TRIPLET DISSOCIATION", " ", 9.6, 28.376, 0.00298220722851047))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.0443", " ", 10, 23.034, 0.00242078380679131))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.0642", " ", 10.5, 26.658, 0.00280165211085538))
    modes.append(modeinfo( "EXC. TRIPLET DISSOCIATION", " ", 10.9, 48.24, 0.00506983636535612))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.1839", " ", 11, 64.976, 0.0068287248688926))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.1073", " ", 11.5, 31.494, 0.0033098969007157))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.0880", " ", 12, 22.558, 0.00237075805824426))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.0304", " ", 12.5, 6.644, 0.000698258557450788))
    modes.append(modeinfo( "EXC. ION PAIR F- + SF5+ F=0.0648", " ", 13, 12.446, 0.00130802618995071))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.1067", " ", 13.5, 17.776, 0.00186818845834516))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.1047", " ", 14, 15.216, 0.00159914241573919))
    modes.append(modeinfo( "EXC. TRIPLET DISSOCIATION", " ", 14.4, 59.75, 0.00627949259597902))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.1211", " ", 14.5, 1.62, 0.000170255698836586))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.2225", " ", 15, 26.256, 0.00275940347447741))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.2731", " ", 15.5, 28.662, 0.00301226471608286))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.1514", " ", 16, 14.674, 0.00154218032390621))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.1831", " ", 16.5, 16.46, 0.00172988197706803))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.1678", " ", 17, 13.58, 0.00142720517913632))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.1098", " ", 17.5, 8.174, 0.000859055606352008))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.0623", " ", 18, 4.322, 0.000454225389118349))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.0361", " ", 18.5, 2.326, 0.000244453552774012))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.0107", " ", 19, 0.7, 0.000073567277275068))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.0129", " ", 19.5, 0.726, 7.62997761452848E-05))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.0053", " ", 20, 0.26, 2.73249887021681E-05))
    modes.append(modeinfo( "EXC. SINGLET DISSOC. F=0.0290", " ", 23, 1.016, 0.000106777648159242))
    # Count TOTALavgN_ALL
    TOTALavgN_ALL = 0.0
    TOTALavgE_ALL = 0.0
    TOTALfracN    = 0.0

    # Add N for each mode
    for mIdx, m in enumerate(modes):
        TOTALavgN_ALL = TOTALavgN_ALL + m.avgQuanta
        TOTALavgE_ALL = TOTALavgE_ALL + m.avgEnergy
        TOTALfracN    = TOTALfracN + m.fracN

    # Calculate probability for each mode
    ProbCutoff = 1.0
    for mIdx, m in enumerate(modes):
        m.probThisModeN = m.avgQuanta/TOTALavgN_ALL
        m.probThisModeE  = m.avgEnergy / TOTALavgE_ALL

        m.probCutHigh     = ProbCutoff
#        m.probCutLow      = ProbCutoff - m.probThisModeN # Same as fracN :)
        m.probCutLow      = ProbCutoff - m.fracN
        ProbCutoff        = m.probCutLow


    if globalvar:
        AVGvisN = 0.0
        AVGvisE = 0.0
        for mIdx, m in enumerate(modes):
            print m.probCutHigh, m.probCutLow

            if mIdx < 12:
                AVGvisN = AVGvisN + m.avgQuanta
                AVGvisE = AVGvisE + m.avgEnergy

        print "\n ... \n ... Testing ...\n ... \n"
        print "Made object with " + str(len(modes)) + " modes...."

        print "Made obect: list of modes. 2 Probabilities are calculated"
        print "Total E = " + str(TOTALavgE_ALL)
        print "Total N = " + str(TOTALavgN_ALL)
        print "Total frac N = " + str(TOTALfracN)

        print "AVG Visible E fraction is = " + str(AVGvisE/TOTALavgE_ALL)
        print "AVG Visible N fraction is = " + str(AVGvisN/TOTALavgN_ALL)


    return modes
#################################################################

#################################################################
#################################################################

# One object to hold info for all events...
AllEventModes = define_modes(1)

probdist = list()

print "Making that the prob is ..."
for j in range (1,Nevents):

    prob = np.random.random()
    probdist.append(prob)


print "Making " + str(Nevents) + "events..."

# lists to account the totals of energy
AllVisibleEs = list()
AllVisibleNs = list()
AllEs = list()
AllNs = list()

for i in range (1,Nevents):

    evt_modes = define_modes(0)

    # assign energy
    evt_energy = 0.0
    evt_remaining = Qvalue

    while evt_energy < Qvalue and evt_remaining > 15:
        # roll die
        prob = np.random.random()
        # print prob
        # assign to one of the modes
        for mIdx, m in enumerate(evt_modes):
            if (prob < m.probCutHigh) and (prob > m.probCutLow):
                # check you have enough energy to assign here
                if ( m.Eloss < evt_remaining ):
                    # add this energy to this mode
                    m.accumulatedN  = m.accumulatedN + 1
                    m.accumulatedE  = m.accumulatedE + m.Eloss
                    # add to the total counts
                    evt_energy    = evt_energy + m.Eloss
                    evt_remaining = evt_remaining - m.Eloss

    VisibleE = 0.0
    VisibleN = 0
    ALLN = 0.0
    ALLE = 0.0

    VisibleEstates = 12 # the 1st 12 are visible E

    for mIdx, m in enumerate(evt_modes):
        if mIdx < VisibleEstates:

            #            print  " now E = " + str(m.accumulatedE) + " now N = " + str(m.accumulatedN) + " - " + m.name + str(m.Eloss)#str(TOTALassignedE)
            VisibleE = VisibleE + m.accumulatedE
            VisibleN = VisibleN + m.accumulatedN
            AllEventModes[mIdx].ALLaccumulatedN.append(m.accumulatedN)
            AllEventModes[mIdx].ALLaccumulatedE.append(m.accumulatedE)
        ALLN = ALLN + m.accumulatedN
        ALLE = ALLE + m.accumulatedE

    AllVisibleEs.append(VisibleE)
    AllVisibleNs.append(VisibleN)

    AllEs.append(ALLE)
    AllNs.append(ALLN)

    if verb:
        for mIdx, m in enumerate(evt_modes):
            print " N_" + str(mIdx) + ": true = " + str(m.fracN)+", sim = " + str(m.accumulatedN/ALLN)

            for mIdx, m in enumerate(evt_modes):
                print " E_" + str(mIdx) + ": true = " + str(m.probThisModeE)+", sim = " + str(m.accumulatedE/ALLE)


            print "............ Visible N = " + str(VisibleN) + \
            " - that is  "+str(VisibleN/ALLN)
            print "............ Visible E = " + str(VisibleE) + \
            " - that is " +str(VisibleE/evt_energy) + " ALLE = " +\
             str(ALLE) + ", evt_energy = " + str(evt_energy)

# #########################################
# # ploting
# #########################################
# Turn interactive plotting off
plt.ioff()

probdist = np.array(probdist)

bin = np.linspace(0,1,Nevents)
plt.hist(probdist, bin)
plt.title("If this is flat you're ok ... now close this with the x")
plt.show()


import datetime
today = datetime.date.today()
date = today.strftime('%Y-%m-%d-%H-%M-%S')
#f"{datetime.datetime.now():%Y-%m-%d-%H-%M-%S}"
# TODO get time to work... now printing 00-00-00

# pngs are empty :(, TODO fix them
# fig = plt.figure()
# fig.savefig("prob_%s.png"%date)
# plt.close(fig)


##############################################
# NOTE: ALLEs and ALLNs can also be plotted
##############################################

AllVisibleEs = np.array(AllVisibleEs)
AllVisibleNs = np.array(AllVisibleNs)

(muEs, sigmaEs) = norm.fit(AllVisibleEs)
(muNs, sigmaNs) = norm.fit(AllVisibleNs)

print "Es[]"
print AllVisibleEs
print "Ns[]"
print AllVisibleNs

binEs = np.linspace(AllVisibleEs[1]*0.995, AllVisibleEs[1]*1.005, 100)
Es = plt
Es.title(r'$\mathrm{Histogram\ of\ E total:}\ \mu=%.3f,\ \sigma=%.3f,\ \delta E/E=%.5f$' %(muEs, sigmaEs, sigmaEs/muEs))
Es.hist(AllVisibleEs, binEs, alpha=0.5, label='Total E', color='hotpink')
Es.legend(loc='upper right')
Es.show()

# fig = Es.figure()
# fig.savefig("E_%s.png"%date)
# Es.close(fig)


binNs = np.linspace(AllVisibleNs[1]*0.95, AllVisibleNs[1]*1.05, 100)
Ns = plt
Ns.title(r'$\mathrm{Histogram\ of\ N total:}\ \mu=%.3f,\ \sigma=%.3f,\ \delta N/N=%.5f, F=%.3f$' %(muNs, sigmaNs, sigmaNs/muNs, (sigmaNs*sigmaNs/muNs)))
Ns.hist(AllVisibleNs, binNs, alpha=0.5, label='Total N', color='hotpink')
Ns.legend(loc='upper right')
Ns.show()

# fig = Ns.figure()
# fig.savefig("N_%s.png"%date)
# Ns.close(fig)

print str(muEs)+", "+str(sigmaEs)+", "+str(muNs)+", "+str(sigmaNs)+", F= "+str(sigmaNs*sigmaNs/muNs)

#########################################
# Plot Ns in fractions of N1
#########################################

(muN0s, sigmaN0s) = norm.fit(AllEventModes[0].ALLaccumulatedN)
(muE0s, sigmaE0s) = norm.fit(AllEventModes[0].ALLaccumulatedE)
for mIdx, ALLm in enumerate(AllEventModes):
    if mIdx < VisibleEstates:
        AllNIsNormalized = ( ALLm.ALLaccumulatedN )/muN0s
        AllEIsNormalized = ( ALLm.ALLaccumulatedE )/muE0s

        (muNis, sigmaNis) = norm.fit(AllNIsNormalized)

        binNi = np.linspace(0, 2, 300)
        Sig = plt
        Sig.title(r'$\mathrm{Histogram\ of\ N \mu=%.3f \ total:}\ \mu=%.3f,\ \sigma=%.3f,\ \delta N/N=%.5f, F=%.3f$' %(mIdx, muNis, sigmaNis, sigmaNis/muNis, (sigmaNis*sigmaNis/muNis)))
        Sig.hist(AllNIsNormalized, binNi, alpha=0.5, label='N/N_0 for this mode', color='gray')
        Sig.legend(loc='upper right')
        Sig.show()
