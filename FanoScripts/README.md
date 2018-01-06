# Scripts for Energy Res and Fano Factor Studies

This is part of a description of what I've learned about the Fano factor and energy resolution, and a description of the code in this directory.
<TODO add Fano reference>

**<span style="color:HotPink" > Fano Factor </span>**

 The Fano factor describes the variation in the ion-electron pairs. Ben describes it as coming from "our ignorance of how the energy is partitioned". It is `F` in the expression of the variation in the energy:

 ```latex
 $$\sigma = \sqrt{F N}$$
```
where `N` is the number of quanta, `E` is the energy spent and `w` the energy \
cost.

For example, for Xe and SF6 (which we expect to be similar to SeF6):

```latex
$F_{Xe} = 0.15 \pm 0.02$

$F_{SF6} = 0.22 $  
```

For SeF6, the energy in the electrons from the decay will go into different ions. Specifically, `$SeF^{+}_{0,..,6}$`, each with different mobility values.
This means that the ions from each of the signal "modes" will reach the readout plane at a different time, which makes it possible to resolve them separately.

Because each mode has a different energy cost (also fraction of the total energy going into it), an "Effective Fano Factor" (`$F_i$`) for each mode can be calculated.
Furthermore, it is possible that accounting the energy going into each state may yield a better energy resolution.

**<span style="color:HotPink" > How do we determine whether this is the case? </span>**


- Understand (i.e. reproduce) the energy allocation into the different states, given by the degrad MC, and reproduce the correct total `F`.

- Extract the `$F_i$` (effective `F`) for each mode from our model.
Uncertainty difference between using  `$\sigma = \sqrt{F N}$` and a sum of the `$\sigma_i = \sqrt{F_i N_i}$`.


# Progress:

**<span style="color:SteelBlue" > ProbabilisticEalloc </span>**
Allocates energy according to probabilities extracted from the number of quanta in each state (according to the degrad output). Produces histograms for `N` and `E` distributions as well as `N_i` and `E_i`. Calculates `F` and `F_i`.

**Issues:**
- The total `F` reproduced with this method does not look correct. It is way too big, sometimes larger than 1.
- Hard-coded tables for modes from Degrad output

---
**<span style="color:HotPink" > Ns and Ps </span>**

We have since learned from the experts that not all the quanta may be accounted for in the processes reported by degrad.

Specifically, Auger cascades produced by K-shell excitations. These will produce multiple electrons. The N reported by those modes might be smaller than the total number of interactions, thus, to get the correct number `N_{expected}` we can try modeling a distribution of unaccounted quanta `P` with no energy cost, to incorporate to these states such that `N_{expected} = N_{reported} + P`

---

**<span style="color:SteelBlue" > FanoCalcInterp </span>**
Allocates energy like **ProbabilisticEalloc**, but using an interp function. Produces distribution of pdf, `E` and `N` per state for one iteration. Produces histograms for `N` and `E`. Calculates `F` and resolution. Other fancy python changes that are way cleaner.

Yields a more sensible Fano factor.

**Issues:**
- Still not reproducing `F` from Degrad.
- Assumes the missing energy (The missing 25% not allocated by Degrad) goes into elastic collisions (which Degrad models with w=0). TODO: check whether this approximation has a sizable impact.
- Needs calculation of `N_i` and `E_i`
