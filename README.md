# ZW_DPS Project — Complete Line-by-Line Documentation

## How This Document Is Organised

Each file is explained in order of execution. For every block of code the explanation covers:
- what the line literally does in the language (C++, Python, bash, Make)
- why that particular syntax was chosen
- the physics phenomenon or engineering decision it encodes
- how it connects to the other files

---

## DATA FLOW OVERVIEW

```
Settings.h
    │
    ├─ GenerateZ.cc  ──→  root/Z_N.root
    ├─ GenerateW.cc  ──→  root/W_N.root
    │
    ├─ GenerateSPS.cc     ──→  root/SPS_N.root
    ├─ GenerateDPS.cc     ──→  root/DPS_N.root   (reads Z.root + W.root)
    ├─ GenerateDPSReal.cc ──→  root/DPSReal_N.root
    │
    ├─ Compare.py         reads SPS + DPS     → plots/
    ├─ CompareDPSModels.py reads DPS + DPSReal → plots/DPSModelComparison/
    ├─ Convergence.py     reads SPS or DPS    → plots/convergence/
    └─ Analyze.py         reads SPS + DPS     → plots/analysis/
```

The Makefile orchestrates compilation and execution. settings_editor.py modifies
Settings.h before a run. run_all.sh is an alternative to make all.

---

# FILE 1 — Settings.h

## Purpose

This is a **C++ header file** shared by all four generators. Instead of
hard-coding beam energy or event counts inside each generator, every
configurable number lives here and is included via `#include "Settings.h"`.
When you change Settings.h and recompile, every generator automatically
picks up the new values.

## Line-by-line

```cpp
#pragma once
```
A **header guard** — ensures this file is included only once per translation
unit even if multiple files include it. Without it, redefinition errors
occur at link time. `#pragma once` is a non-standard but universally
supported shorthand for the classic `#ifndef / #define / #endif` guard.

---

```cpp
constexpr double ECM_GEV = 13600.;
```
**Physics**: The LHC Run 3 centre-of-mass energy is 13.6 TeV = 13,600 GeV.
`constexpr` means the value is evaluated at compile time, not runtime — the
compiler can embed it directly into machine code wherever it is used, giving
zero overhead versus a literal constant. The `.` suffix makes it a `double`,
avoiding integer division issues if it is ever used in arithmetic.

---

```cpp
constexpr double Z_MASS_MIN = 75.;
constexpr double Z_MASS_MAX = 120.;
```
**Physics**: The Z boson has mass 91.1876 GeV. At low invariant mass (below
~75 GeV) the Drell-Yan process is dominated not by the Z resonance but by
the virtual photon γ* — a different physics process. Restricting the hard
process to 75–120 GeV (the "Z pole region") removes γ* contamination and
ensures both the SPS and DPS samples have the same Z mass spectrum. This is
critical for a fair comparison: without the window the DPS-Mixed Z sample
(from GenerateZ alone) would include low-mass γ* events that the SPS sample
(which includes a W) does not. The window is passed to PYTHIA as
`PhaseSpace:mHatMin` and `PhaseSpace:mHatMax`.

---

```cpp
constexpr int N_EVENTS_Z   = 30000;
constexpr int N_EVENTS_W   = 30000;
constexpr int N_EVENTS_SPS = 30000;
constexpr int N_EVENTS_DPS = 30000;
```
Default sample sizes. The Makefile reads these at make-time using grep to
stamp file and binary names, so `Z_30000.root` and `SPS_30000.root` never
collide when you run different N. The int type is appropriate since event
counts are always whole numbers. settings_editor.py can modify these without
manual text editing.

---

```cpp
constexpr bool ENABLE_MPI = true;
```
**Physics**: Multi-Parton Interactions (MPI) are additional soft parton-parton
scatterings that occur in the same pp collision alongside the primary hard
scatter. They create the "underlying event" — extra charged particle activity
that fills the detector. For a physics study they must be ON to produce
realistic event topologies. Setting this to false speeds up generation by
~20% and is only useful for debugging.

---

```cpp
#include "Pythia8/Pythia.h"
#include <cmath>
#include <string>
```
These includes are placed inside Settings.h rather than in each .cc file
because the inline helper functions below use Pythia8 types. Since Settings.h
uses `#pragma once`, including it multiple times is safe.

---

```cpp
inline void setQuiet(Pythia8::Pythia& p) {
    p.readString("Print:quiet      = on");
    p.readString("Next:numberCount = 0");
}
```
**PYTHIA mechanism**: `Pythia::readString()` passes configuration strings to
PYTHIA's internal parameter database before `pythia.init()`. These two
settings suppress: (1) the verbose per-event header lines PYTHIA prints by
default, and (2) the progress messages like "1000 events generated". Without
them, terminals fill with tens of thousands of lines during a 30,000-event
run. `inline` tells the compiler to expand the function body at the call site
rather than making a function call — appropriate for small utility functions
called once per program.

---

```cpp
inline void setBeam(Pythia8::Pythia& p) {
    p.readString("Beams:idA = 2212");  // proton PDG id
    p.readString("Beams:idB = 2212");  // proton PDG id
    p.readString("Beams:eCM = " + std::to_string(ECM_GEV));
}
```
**Physics**: PDG id 2212 is the proton. `eCM` is the total centre-of-mass
energy in GeV. `std::to_string(ECM_GEV)` converts the double 13600.0 to the
string "13600.000000" — PYTHIA's `readString` expects a string, not a number.
Centralising this in Settings.h ensures every generator uses identical beam
conditions, which is essential for comparing SPS vs DPS: if their beams
differed, the comparison would be meaningless.

---

```cpp
inline void setMPI(Pythia8::Pythia& p) {
    p.readString(std::string("PartonLevel:MPI = ") + (ENABLE_MPI ? "on" : "off"));
}
```
The ternary operator selects "on" or "off" based on the compile-time constant
`ENABLE_MPI`. String concatenation builds the full PYTHIA configuration
string. `PartonLevel:MPI` controls whether PYTHIA generates additional soft
scatters ("MPI" = Multi-Parton Interactions) on top of the primary hard
process. This is distinct from `SecondHard`, which generates a second HARD
scatter — MPI generates many soft scatters.

---

```cpp
inline void setZDecayMuMu(Pythia8::Pythia& p) {
    p.readString("23:onMode    = off");
    p.readString("23:onIfMatch = 13 -13");
}
```
**Physics**: PDG id 23 is the Z boson. By default PYTHIA lets it decay
according to Standard Model branching ratios: ~20% to charged leptons, ~70%
to hadrons, ~20% to neutrinos. For this analysis we need only the muon
channel (Z → μ⁺μ⁻) because:
1. Muons are detected with high efficiency at LHC experiments
2. Both daughters are fully reconstructable (unlike the neutrino in W→μν)
3. The invariant mass reconstruction m(μ⁺μ⁻) directly gives the Z mass

`23:onMode = off` turns OFF all Z decay modes. `23:onIfMatch = 13 -13`
re-enables only the channel where the two daughters have PDG ids 13 (μ⁻)
and -13 (μ⁺). This forces 100% of Z's to decay this way.

---

```cpp
inline void setWDecayMuNu(Pythia8::Pythia& p) {
    p.readString("24:onMode      = off");
    p.readString("24:onIfMatch   = 13 -14");   // W+  → μ⁻ ν̄_μ  WAIT: see note
    p.readString("-24:onMode     = off");
    p.readString("-24:onIfMatch  = -13 14");   // W-  → μ⁺ ν_μ
}
```
PDG id 24 = W⁺, -24 = W⁻. The decay channels are:
- W⁺ → μ⁺ ν_μ means daughters have ids -13 (μ⁺) and 14 (ν_μ)
- W⁻ → μ⁻ ν̄_μ means daughters have ids 13 (μ⁻) and -14 (ν̄_μ)
Note PYTHIA uses `onIfMatch` with the daughter PDG IDs. `13 -14` means
one daughter is a μ⁻ (13) and the other is an anti-muon-neutrino (-14) —
that is the W⁺ decay. The muon-channel W decay is preferred because the
muon is identified and its pT measured, giving access to the W transverse
mass MT = √(2 pT(μ) pT(ν) (1 − cos Δφ(μ,ν))).

---

```cpp
inline double dPhi(double phi1, double phi2) {
    double d = std::fabs(phi1 - phi2);
    return (d > M_PI) ? 2.0 * M_PI - d : d;
}
```
**Physics**: Azimuthal angle φ is defined on [0, 2π) or equivalently
on (−π, π]. The raw difference φ₁ − φ₂ can be any value in (−2π, 2π).
The physically meaningful angular separation is the MINIMUM arc between the
two directions, which lies in [0, π]. If the raw |Δφ| > π, we subtract it
from 2π to get the shorter arc. `M_PI` is π from `<cmath>`. This function is
used in every generator to compute |Δφ(W,Z)|.

---

# FILE 2 — GenerateZ.cc

## Purpose

Generates standalone Z bosons via Drell-Yan, decays them to μ⁺μ⁻, and saves
their kinematics to `root/Z_N.root`. This file serves two roles:
1. Provides the Z sample for GenerateDPS.cc (event mixing)
2. Serves as a validation sample to confirm the Z mass distribution peaks at 91 GeV

## Line-by-line

```cpp
#include "Settings.h"
```
Pulls in all shared constants (ECM_GEV, N_EVENTS_Z, PDG ids, helper functions)
without duplicating them. Because Settings.h includes Pythia8/Pythia.h, this
single include gives access to the full PYTHIA API.

---

```cpp
#include "TFile.h"
#include "TTree.h"
```
ROOT classes for writing output. TFile is an abstract ROOT file container.
TTree is ROOT's columnar data structure — analogous to a database table where
each column (branch) stores one variable per event, and each row is one event.

---

```cpp
int nEvents = N_EVENTS_Z;
if (argc > 1) nEvents = atoi(argv[1]);
```
`argc` counts command-line arguments. If the user runs `./GenerateZ 50000`,
then `argc = 2` and `argv[1] = "50000"`. `atoi` converts the string to int.
This allows the Makefile to compile one binary stamped with the default N
but allows `make all N=50000` style overrides. If no argument is given,
the compile-time constant N_EVENTS_Z is used.

---

```cpp
string outPath = "root/Z_" + to_string(nEvents) + ".root";
```
The output filename encodes N: `root/Z_30000.root`. This means running with
different N values never overwrites previous results. The `root/` prefix means
the file is created in the `root/` subdirectory, which the Makefile creates
with `mkdir -p $(ROOT)` before running the binary.

---

```cpp
pythia.readString("WeakSingleBoson:all           = off");
pythia.readString("WeakSingleBoson:ffbar2gmZ     = on");
pythia.readString("WeakBosonAndParton:qqbar2gmZg = on");
pythia.readString("WeakBosonAndParton:qg2gmZq    = on");
```
**Physics**: This configures the hard process. `WeakSingleBoson:all = off`
disables all single-boson production (both W and Z). Then three specific
channels are re-enabled:

- `ffbar2gmZ`: q q̄ → γ*/Z (s-channel) — the primary Drell-Yan process.
  An antiquark and quark annihilate through a virtual photon/Z propagator.
  This is the dominant contribution near the Z pole.

- `qqbar2gmZg`: q q̄ → γ*/Z + g — adds a gluon in the final state.
  This gives the Z non-zero transverse momentum at leading order, which is
  necessary for realistic pT(Z) distributions.

- `qg2gmZq`: q g → γ*/Z + q — a "Compton-like" process where a gluon from
  the proton sea converts to a quark while emitting a Z. Also gives Z pT.

Without the last two, all Z bosons would have pT = 0 at leading order —
unphysical and in contradiction with LHC data.

---

```cpp
pythia.readString("PhaseSpace:mHatMin = " + to_string(Z_MASS_MIN));
pythia.readString("PhaseSpace:mHatMax = " + to_string(Z_MASS_MAX));
```
**Physics**: `mHat` is the invariant mass of the partonic centre-of-mass
system, equal to the Z/γ* mass. Restricting to 75–120 GeV:
- removes the low-mass Drell-Yan continuum from virtual photons
- ensures this Z sample matches the mass window used by GenerateDPSReal

The string `+ to_string(Z_MASS_MIN)` converts 75.0 to "75.000000" at
runtime and concatenates it with the parameter name string.

---

```cpp
TTree* tree = new TTree("ZTree", "Standalone Z boson sample");
```
Creates a ROOT TTree named "ZTree". The first argument is the internal ROOT
object name (used when loading: `f->Get("ZTree")`). The second is a
human-readable description stored in the file metadata. Named "ZTree" rather
than "Events" to distinguish it from the WZ event trees (which use "Events")
— GenerateDPS.cc loads it specifically as `fZ->Get("ZTree")`.

---

```cpp
double Z_pt, Z_eta, Z_phi, Z_mass, Z_rapidity;
double Z_px, Z_py, Z_pz, Z_E, Z_mumu_mass;
double x1, x2, Q2, weight;
```
One C++ variable per quantity to be saved. ROOT TTrees write the value of
these variables at the time of `tree->Fill()`. The variables are declared
here (not inside the loop) because ROOT's `Branch()` stores a POINTER to
the variable — that pointer must remain valid for the lifetime of the tree.

**Physics meaning of each**:
- `Z_pt`: transverse momentum = √(px²+py²) [GeV]
- `Z_eta`: pseudorapidity = −ln(tan(θ/2)) where θ is the polar angle
- `Z_phi`: azimuthal angle in the transverse plane [rad]
- `Z_mass`: true invariant mass from PYTHIA's internal bookkeeping [GeV]
- `Z_rapidity`: true rapidity = ½ ln((E+pz)/(E−pz))
- `Z_px, Z_py, Z_pz, Z_E`: four-momentum components needed to reconstruct
  the Z from its daughters or to combine with W in DPS mixing
- `Z_mumu_mass`: reconstructed mass from (p_μ⁺ + p_μ⁻)² — should equal
  Z_mass as a sanity check on the four-momentum reconstruction
- `x1, x2`: Bjorken-x values of the two incoming partons — fraction of
  proton momentum carried by the parton. x1 × x2 × s = ŝ (partonic CM energy²)
- `Q2`: factorisation scale squared — the scale at which the PDF is evaluated
- `weight`: PYTHIA event weight, needed for weighted histograms

---

```cpp
tree->Branch("Z_pt", &Z_pt);
```
Registers the variable `Z_pt` as a branch of the TTree. `&Z_pt` is the
memory address of the variable. ROOT will read the current value at that
address each time `tree->Fill()` is called. The branch name "Z_pt" is the
column name used when loading with uproot: `tree["Z_pt"].array()`.

---

```cpp
int iZ = -1;
for (int i = 0; i < pythia.event.size(); ++i)
    if (pythia.event[i].id() == ID_Z) iZ = i;
if (iZ < 0) { ++nNoZ; continue; }
```
**PYTHIA event record**: `pythia.event` is a list of all particles in the
event — incoming protons, their partons, intermediate states, shower
products, and final-state hadrons. Each particle has an `id()` (PDG id) and
a `status()`. The Z appears multiple times: once when produced, once for
each copy made by colour reconnection / shower bookkeeping. The LAST copy
with id == 23 is the Z just before it decays. We overwrite `iZ` in each
iteration, so after the loop it holds the index of the last Z — which is
what we want. If no Z is found (`iZ < 0`), we skip the event.

---

```cpp
int d1 = pythia.event[iZ].daughter1();
int d2 = pythia.event[iZ].daughter2();
if (d1 <= 0 || d2 <= 0) { ++nNoDau; continue; }
```
`daughter1()` and `daughter2()` return the indices in the event record of the
Z's decay products. We forced Z → μ⁺μ⁻, so d1 and d2 should point to the
muon and antimuon. A return value ≤ 0 means no daughters were found (rare,
but possible if the event record was truncated or the Z is from an earlier
intermediate state). We skip such events.

---

```cpp
Vec4 pReco = pythia.event[d1].p() + pythia.event[d2].p();
Z_mumu_mass = pReco.mCalc();
```
`Vec4` is PYTHIA's four-vector class. `p()` returns the four-momentum
(px, py, pz, E) of particle d1. Adding two Vec4 objects adds their
components. `mCalc()` computes the invariant mass = √(E²−|p|²). This
reconstructed mass should match `pZ.m()` to within floating-point precision —
if it doesn't, something is wrong with the event record bookkeeping.

---

```cpp
x1     = pythia.info.x1();
x2     = pythia.info.x2();
Q2     = pythia.info.Q2Fac();
weight = pythia.info.weight();
```
`pythia.info` is a PYTHIA object that stores event-level metadata for the
current event. `x1()` and `x2()` give the Bjorken-x of the two incoming
partons. `Q2Fac()` gives the factorisation scale squared — the scale at
which the PDF was evaluated for this event. These quantities are standard
outputs of any parton-level event generator.

---

```cpp
pythia.stat();
```
After the event loop, prints PYTHIA's statistics table to stdout: total
cross section, number of events accepted and rejected, and any errors/warnings
that occurred. The cross section (σ) is printed in mb (millibarns) and can
be used to normalise distributions to an absolute yield at a given luminosity.

---

# FILE 3 — GenerateW.cc

Essentially identical structure to GenerateZ.cc. Key differences:

```cpp
pythia.readString("WeakSingleBoson:ffbar2W     = on");
pythia.readString("WeakBosonAndParton:qqbar2Wg = on");
pythia.readString("WeakBosonAndParton:qg2Wq    = on");
```
**Physics**: `ffbar2W` = q q̄' → W (charged-current Drell-Yan). Note the
prime on q̄' — this is a flavour-changing weak interaction: an up-type quark
can annihilate with an anti-down-type quark (and similarly for other CKM
combinations) to produce a W. No mass window is set because the W mass
(~80.4 GeV) is the only relevant scale and PYTHIA generates it correctly
from the Breit-Wigner lineshape without additional cuts.

```cpp
double MT, x1, x2, Q2, weight;
```
**Physics**: MT is the W transverse mass:
  MT = √(2 pT(μ) pT(ν) (1 − cos Δφ(μ,ν)))
This is the experimentally accessible approximation to the W mass when the
neutrino pT is inferred from missing transverse energy. It provides a Jacobian
peak at MT = MW.

```cpp
int d1 = pW.daughter1(), d2 = pW.daughter2();
for (int di : {d1, d2}) {
    int aid = abs(pythia.event[di].id());
    if (aid == 13) { /* fill mu_pt etc */ }
    if (aid == 14) { /* fill nu_pt etc */ }
}
```
The W daughters are a muon (|id|=13) and a neutrino (|id|=14). Taking the
absolute value handles both W⁺ (daughters: μ⁺=-13, ν=14) and W⁻ (daughters:
μ⁻=13, ν̄=-14). The `{d1, d2}` initialiser list in the range-for loop is a
C++11 feature that iterates over exactly two values without allocating a
vector.

The output tree is named "WTree" to distinguish it from "Events" trees.
GenerateDPS.cc loads it as `fW->Get("WTree")`.

---

# FILE 4 — GenerateSPS.cc

## Purpose

Generates WZ events via Single Parton Scattering: one hard q q̄' → WZ
matrix element per pp collision. This is the Standard Model prediction at
LO. Both W and Z come from the same parton pair — they are kinematically
correlated by momentum conservation.

## Key configuration

```cpp
pythia.readString("WeakDoubleBoson:all      = off");
pythia.readString("WeakDoubleBoson:ffbar2ZW = on");
pythia.readString("SecondHard:generate      = off");
```
**Physics**: `WeakDoubleBoson` is the class of processes that produce two
electroweak bosons from one partonic interaction. `ffbar2ZW` activates the
channel q q̄' → Z W — the single most important WZ production channel at the
LHC. The W and Z are both present from the beginning in this channel, not
from two separate scatterings.

`SecondHard:generate = off` is the critical SPS switch — it ensures no
second hard scatter is generated. This is what makes it "single parton
scattering": exactly one hard partonic interaction per pp collision.

No mass window is set here (unlike GenerateZ) because `ffbar2ZW` already
places both bosons at their physical masses through the WZ propagator
structure, and cutting on mHat would remove events with large boson recoil.

## Event loop: finding W and Z

```cpp
int iZ = -1;
for (int i = 0; i < pythia.event.size(); ++i)
    if (pythia.event[i].id() == ID_Z) iZ = i;

int iW = -1;
W_charge = 0;
for (int i = 0; i < pythia.event.size(); ++i) {
    int id = pythia.event[i].id();
    if (id == ID_WPLUS || id == ID_WMINUS) {
        iW = i;
        W_charge = (id == ID_WPLUS) ? +1 : -1;
    }
}
```
In the SPS event record both Z (id=23) and W (id=24 or -24) appear together.
Two separate loops find the last copy of each. The charge is extracted from
the sign of the W id.

## Derived observables

```cpp
dEta_WZ = fabs(W_eta - Z_eta);
dPhi_WZ = dPhi(W_phi, Z_phi);
dY_WZ   = fabs(W_rapidity - Z_rapidity);

Vec4 pWZ = pW.p() + pZ.p();
M_WZ  = pWZ.mCalc();
pT_WZ = pWZ.pT();
sT    = W_pt + Z_pt;
```
**Physics of each**:

- `dEta_WZ = |η_W − η_Z|`: Pseudorapidity separation. In SPS, both bosons
  come from the same partonic system, so they share the same longitudinal
  boost y_pair = ½ln(x₁/x₂). This means η_W and η_Z are correlated →
  small dEta is preferred in SPS.

- `dPhi_WZ = |φ_W − φ_Z|` (folded to [0,π]): Azimuthal separation. At LO
  the SPS process q q̄' → WZ is a 2→2 process. Momentum conservation in the
  transverse plane requires the WZ system to have zero pT. But ISR (initial
  state radiation) gives the WZ system a net pT, and the W and Z must
  balance it. This pushes them toward back-to-back (dPhi → π). This is the
  most powerful SPS-vs-DPS discriminator.

- `dY_WZ = |y_W − y_Z|`: True rapidity separation. Complementary to dEta.

- `M_WZ = mCalc()`: Invariant mass of the WZ system. `Vec4::mCalc()` computes
  √((E_W+E_Z)² − |p_W+p_Z|²). In SPS this has a threshold at M_W + M_Z ≈
  170 GeV and falls steeply — the full partonic energy must go into making
  both bosons.

- `pT_WZ = pWZ.pT()`: Transverse momentum of the WZ system = √(px_WZ² + py_WZ²).
  At LO this is zero (momentum conservation in 2→2), but ISR gives it
  a non-zero value.

- `sT = W_pt + Z_pt`: Scalar pT sum, a DPS discriminator. In DPS the two
  bosons are independent, so sT is the sum of two independent pT values.
  In SPS at LO, pT(W) ≈ pT(Z) (back-to-back), so sT ≈ 2 × pT_WZ.

```cpp
nJets30 = 0;
for (int i = 0; i < pythia.event.size(); ++i) {
    const Particle& p = pythia.event[i];
    if (!p.isFinal()) continue;
    if (abs(p.id()) > 8 && p.id() != 21) continue;
    if (p.pT() > 30.) ++nJets30;
}
```
Counts final-state quarks (|id| ≤ 8) and gluons (id=21) with pT > 30 GeV
as a crude jet proxy. `isFinal()` returns true for particles that are not
decayed or rescattered — the physical final state. This is saved to the ROOT
tree for possible future jet studies but is not used in the current Python
analysis.

---

# FILE 5 — GenerateDPS.cc (Event Mixing)

## Purpose and physics model

This implements the **simplest DPS approximation**: take a Z event from
`root/Z_N.root` and a W event from `root/W_N.root`, both from completely
separate PYTHIA runs, and combine them into a "DPS event" by simply adding
their four-momenta. The two bosons have **zero** kinematic correlation by
construction.

**Why this is a valid DPS approximation**: In real DPS, the W and Z come from
two independent hard scatters within the same pp collision. If we ignore the
shared underlying event and shared beam remnants, the two bosons behave as
if they came from completely separate events. Event mixing captures this
leading-order approximation exactly: Δφ(W,Z) is exactly uniform on [0,π],
which is the theoretical DPS prediction in the absence of any correlations.

**Limitation**: This model over-decorrelates. In a real DPS event both bosons
share the same pp collision, so they are connected by energy-momentum
conservation and the shared underlying event. These are small but real
correlations absent in event mixing. GenerateDPSReal.cc addresses this.

## Key code

```cpp
TFile* fW = TFile::Open(wPath.c_str(), "READ");
TTree* tW = (TTree*)fW->Get("WTree");
TFile* fZ = TFile::Open(zPath.c_str(), "READ");
TTree* tZ = (TTree*)fZ->Get("ZTree");
```
Opens the ROOT files produced by GenerateW and GenerateZ. `TFile::Open` in
"READ" mode prevents accidental modification. The tree names "WTree" and
"ZTree" must match what was used in GenerateW and GenerateZ — this is the
cross-file contract between the generators.

```cpp
long nEvents = min({(long)nReq, nW, nZ});
```
We can only mix as many pairs as we have events in the smaller of the two
samples. The `min` with an initialiser list is a C++11 feature that takes
the minimum of three values simultaneously.

```cpp
tW->SetBranchAddress("W_pt", &W_pt);
```
Connects the branch "W_pt" in the TTree to the local C++ variable W_pt.
After calling `tW->GetEntry(i)`, the variable W_pt is populated with the
value from event i in the tree. This is how ROOT TTrees are read in C++.

```cpp
for (long i = 0; i < nEvents; ++i) {
    tW->GetEntry(i);
    tZ->GetEntry(i);
```
Reads event i from both trees. The W event and Z event at index i are from
completely separate PYTHIA runs — this is the "mixing". There is no physical
reason event 42 of W and event 42 of Z should be correlated, which is exactly
the DPS-mixing approximation.

```cpp
TLorentzVector pW(W_px, W_py, W_pz, W_E);
TLorentzVector pZ(Z_px, Z_py, Z_pz, Z_E);
TLorentzVector pWZ = pW + pZ;
M_WZ    = pWZ.M();
pT_WZ   = pWZ.Pt();
```
`TLorentzVector` is ROOT's four-vector class (as opposed to PYTHIA's `Vec4`).
It is used here because PYTHIA is not running — we only have ROOT. The four-
momentum components stored in the ROOT file are used directly. `M()` computes
the invariant mass, `Pt()` the transverse momentum.

```cpp
weight  = 1.0;
nJets30 = 0;
```
Since these are mixed events, not real PYTHIA events, there is no event weight
to extract and no jet information. Weight = 1 means all events contribute
equally to histograms, consistent with treating each mixed pair as equally
probable.

---

# FILE 6 — GenerateDPSReal.cc

## Purpose and physics model

This implements "real DPS" by generating BOTH the W and Z within the SAME
pp collision using PYTHIA's `SecondHard` machinery. This is physically more
correct than event mixing: the two bosons share the same underlying event
(MPI activity, beam remnants, ISR), which introduces small correlations.

## PYTHIA SecondHard machinery

```cpp
pythia.readString("WeakSingleBoson:all             = off");
pythia.readString("WeakSingleBoson:ffbar2gmZ       = on");
pythia.readString("WeakBosonAndParton:qqbar2gmZg   = on");
pythia.readString("WeakBosonAndParton:qg2gmZq      = on");
```
This configures the **first hard scatter** (same as GenerateZ). In PYTHIA's
SecondHard framework, the first hard scatter is the primary hard process
defined by these settings.

```cpp
pythia.readString("PhaseSpace:mHatMin = " + to_string(Z_MASS_MIN));
pythia.readString("PhaseSpace:mHatMax = " + to_string(Z_MASS_MAX));
```
The mass window applies to the first hard scatter (the Z production).
This is critical: it ensures the Z in DPSReal has the same mass spectrum
as the Z in DPSMixed (where the mass window is applied in GenerateZ.cc).
Without this, the two DPS models would be comparing different Z mass spectra.

```cpp
pythia.readString("SecondHard:generate = on");
pythia.readString("SecondHard:SingleW  = on");
pythia.readString("SecondHard:WAndJet  = on");
```
**PYTHIA SecondHard mechanism**: This activates PYTHIA's documented two-hard-
scatter framework. It generates a second independent hard scatter (the W
production) within the same pp collision after the first (the Z) is generated.
The second scatter uses its own PDF sampling and its own parton pair, but from
the same proton. The two scatters share:
- Energy-momentum conservation of the full pp system
- The same underlying event (MPI, beam remnants)
- The same ISR structure

`SecondHard:SingleW` = qq̄' → W (bare W)
`SecondHard:WAndJet` = qq̄' → Wg, qg → Wq (W + associated jet)

Using both gives a realistic W pT spectrum from the second scatter.

## Finding W and Z in the same event

The search code is identical to GenerateSPS.cc:
```cpp
int iZ = -1;
for (int i = 0; i < pythia.event.size(); ++i)
    if (pythia.event[i].id() == ID_Z) iZ = i;

int iW = -1;
for (int i = 0; i < pythia.event.size(); ++i) {
    int id = pythia.event[i].id();
    if (id == ID_WPLUS || id == ID_WMINUS) {
        iW = i; W_charge = (id == ID_WPLUS) ? +1 : -1;
    }
}
```
In a SecondHard event, both the Z (from scatter 1) and the W (from scatter 2)
appear in the same event record. The same last-copy search strategy finds both.
This is a key difference from DPS-Mixed where Z and W come from separate
files: here they are in the same PYTHIA event.

## Physics of SecondHard vs event mixing

The DPS cross section formula is:
  σ_DPS(WZ) = σ(W) × σ(Z) / (2 × σ_eff)

where σ_eff ≈ 15 mb encodes the transverse overlap of parton pairs. PYTHIA's
SecondHard machinery samples the second scatter according to the cross section
of the specified process. The DPSReal sample is therefore a physically
motivated simulation, not an approximation.

The comparison DPSMixed vs DPSReal is a physics result in itself: if they
agree, event mixing is validated as a good approximation; if they differ,
the differences quantify shared-event correlations.

---

# FILE 7 — Compare.py

## Purpose

Reads SPS and DPS(Mixed) ROOT files and produces two types of output:
1. Individual PDF plots, one per observable
2. A summary 8-panel figure

## Observable catalogue

```python
ALL_OBSERVABLES = [
    ("dEta_WZ", r"$|\Delta\eta(W,Z)|$", "", (0,8), False, "dEta", True),
    ...
    ("Z_pt",    r"$p_T(Z)$",            "GeV", (0,500), True, "ZpT", False),
    ...
]
```
Each entry is a 7-tuple:
1. Branch name in the ROOT tree — used by uproot to load the array
2. LaTeX label — displayed on plot axes, uses matplotlib's TeX renderer
3. Unit string — appended to axis label in brackets
4. (xlo, xhi) — x-axis range
5. log_y — whether to use logarithmic y-axis (True for pT distributions,
   which span several orders of magnitude)
6. Short tag — used in filenames and table column headers
7. use_kde — True for angular observables, False for pT observables

The last column is the key innovation: angular observables (dEta, dPhi, dY)
and M_WZ use KDE (smooth curves), while pT observables (W_pt, Z_pt, sT,
pT_WZ) use normalised histograms. This separates the two numerics methods
because KDE was failing catastrophically for pT: the DPS pT distribution has
a hard tail from SecondHard producing rare high-pT events, which makes
Silverman's rule choose an enormous bandwidth, collapsing the estimated density
to values of 10^14 in some bins.

## Style configuration

```python
plt.rcParams.update({...})
```
Modifies matplotlib's global style dictionary. `pdf.fonttype: 42` ensures
fonts are embedded as outlines (Type 1 compatible) in PDFs, required by
journals. `xtick.top: True, ytick.right: True` adds axis tick marks on all
four sides, the standard in HEP publications. `xtick.direction: in` puts
ticks inside the plot area. `figure.dpi: 180, savefig.dpi: 300` separates
screen display resolution from file save resolution.

## KDE density estimator

```python
def _kde_density(vals, lo, hi):
    xs = np.linspace(lo, hi, N_KDE_PTS)
    try:
        ys = gaussian_kde(vals)(xs)
    except Exception:
        ys = np.zeros(N_KDE_PTS)
    return xs, ys
```
`scipy.stats.gaussian_kde` constructs a Gaussian kernel density estimator.
Each data point contributes a Gaussian bell curve of width h (bandwidth).
The bandwidth h is chosen by Silverman's rule: h = 0.9 σ N^(-1/5). The
resulting estimate is the sum of all these Gaussians, normalised to integrate
to 1. The `try/except` catches the case where `vals` has only one unique
value (kde would raise a LinAlgError).

## Histogram density estimator

```python
def _hist_density(vals, lo, hi, nbins=N_HIST_BINS):
    edges   = np.linspace(lo, hi, nbins + 1)
    bw      = edges[1] - edges[0]
    counts, _ = np.histogram(vals, bins=edges)
    total   = max(counts.sum(), 1)
    density = counts / (total * bw)
    err     = np.sqrt(counts) / (total * bw)
    centres = 0.5 * (edges[:-1] + edges[1:])
    return centres, density, err
```
`np.histogram` bins the data. Dividing by `total * bw` normalises so the
histogram integrates to 1 (probability density). `np.sqrt(counts)` is the
Poisson uncertainty on each bin count — the standard statistical error for
counting experiments. Dividing by `total * bw` converts it to the density
error. The function returns bin centres (midpoints of edges) for plotting.

## File discovery

```python
def _paired_ns(root_dir):
    sps_files  = {_extract_n(f): f for f in glob.glob(..., "SPS_*.root")}
    dpsmix_files = {_extract_n(f): f for f in glob.glob(..., "DPSMixed_*.root")}
    common = sorted(set(sps_files) & set(dpsmix_files))
```
`glob.glob` finds all matching filenames. A dict comprehension `{n: f}` maps
the extracted N value to the filepath. `set(sps_files) & set(dpsmix_files)`
finds N values that appear in both — the intersection. Only paired samples
can be compared. `sorted()` orders them ascending by N.

## Non-interactive mode detection

```python
if not sys.stdin.isatty():
```
`sys.stdin.isatty()` returns True if stdin is a terminal (interactive shell)
and False if stdin is a pipe or file (e.g., `make all` redirects stdin from
/dev/null). This allows the same script to run interactively (ask the user
what to plot) or non-interactively (auto-pick the most recent pair) from
the Makefile.

## KS test

```python
ks, pv = ks_2samp(s_vals, d_vals)
ax.text(0.03, 0.03, f"KS = {ks:.3f}\n$p = {pv:.1e}$", ...)
```
`scipy.stats.ks_2samp` computes the Kolmogorov-Smirnov two-sample test
statistic: D = max_x |F_SPS(x) − F_DPS(x)| where F is the empirical CDF.
D ranges from 0 (identical distributions) to 1 (no overlap). The p-value
is the probability of observing D as large as this if both samples came from
the same distribution. `p ≈ 0` means the two distributions are statistically
incompatible. The annotation is placed in the lower-left corner of each axis.

---

# FILE 8 — CompareDPSModels.py

## Purpose

Compares DPSMixed (event mixing) vs DPSReal (SecondHard) — two different
physics approximations to Double Parton Scattering.

## Key difference from Compare.py

Compare.py compares SPS vs DPS. CompareDPSModels.py compares two DPS
implementations. If DPSReal ≈ DPSMixed for all observables, the event mixing
approximation is validated. If they differ, the differences quantify the
contribution of shared-event structure to DPS correlations.

## Colour scheme

```python
C_MIXED = "#2ca02c"   # green
C_REAL  = "#ff7f0e"   # orange
```
Different colours from Compare.py (blue/red) to prevent visual confusion when
looking at multiple plots from different scripts.

## Two-panel plots with ratio

```python
def comparison_panel_with_ratio(fig, gs_top, gs_bot, mixed_vals, real_vals, ...):
    ax_top = fig.add_subplot(gs_top)
    ax_bot = fig.add_subplot(gs_bot, sharex=ax_top)
```
Each observable gets TWO subplots vertically: top shows both distributions,
bottom shows their ratio DPSReal/DPSMixed. `sharex=ax_top` links the x-axes
so zooming in the top panel also zooms the bottom. This layout is standard
in HEP: the "ratio panel" immediately shows which x-ranges differ.

## Ratio computation

```python
ratio = np.divide(yr, ym,
                  out=np.full_like(yr, np.nan),
                  where=(ym > 0))
```
`np.divide` with `out` and `where` performs division only where `ym > 0`,
filling other elements with NaN. This avoids division-by-zero warnings.
NaN values are automatically excluded by matplotlib (they appear as gaps
in the plot rather than infinity spikes).

## Output directory

```python
plots_dir = os.path.join(script_dir, "plots", "DPSModelComparison")
os.makedirs(plots_dir, exist_ok=True)
```
`exist_ok=True` prevents an error if the directory already exists (unlike
`os.mkdir`). The separate subdirectory keeps DPS model comparison plots
isolated from the SPS vs DPS comparison plots.

---

# FILE 9 — Convergence.py

## Purpose

Studies how distributions change as a function of sample size N. For each
selected observable, it overlays KDE curves from N = 10k, 20k, 50k ... 800k.
When successive curves overlap, the distribution has "converged" — more
events won't change the shape.

## Alpha fading and Z-ordering

```python
alpha_val = 0.45 + 0.55 * (i / max(1, n_curves - 1))
z_order   = 10 + i
ax.plot(..., alpha=alpha_val, zorder=z_order, ...)
```
Low-N curves (most noisy) are plotted at 45% opacity. High-N curves (most
trustworthy) are at 100% opacity and rendered last (highest zorder = on top).
This makes convergence visually obvious: the converged orange curve (high N)
is bright and solid; the purple low-N curve is faint underneath.

## Plasma colormap

```python
cmap  = plt.get_cmap("plasma")
lo, hi = 0.10, 0.88
return [cmap(lo + (hi - lo) * i / (n_items - 1)) for i in range(n_items)]
```
`plasma` is a perceptually uniform sequential colormap: equal visual steps
correspond to equal value steps. Sampling from 0.10 to 0.88 (not 0 to 1)
avoids the near-black start and near-white yellow end of the colormap, both
of which are hard to see on white paper.

## Change detection cache

```python
def _fp(paths):
    out = {}
    for p in paths:
        st = os.stat(p)
        out[os.path.basename(p)] = [st.st_mtime_ns, st.st_size]
    return out
```
`os.stat` returns file metadata. `st_mtime_ns` is the last modification time
in nanoseconds (more precise than st_mtime in seconds). `st_size` is the file
size in bytes. Storing both guards against the case where a file is modified
but happens to be the same size.

```python
def _needs_replot(out_path, source_paths):
    new_fp = _fp(source_paths)
    if key not in cache: return True, "no previous fingerprint"
    for name, (mt, sz) in new_fp.items():
        if name in old_fp and [mt, sz] != old_fp[name]:
            return True, f"{name} has been modified"
    return False, "no changes detected"
```
Compares current fingerprints against the cached fingerprints from the
previous run. If any file has changed (mtime or size differs), replot.
If nothing changed, ask the user. This prevents accidentally regenerating
expensive plots from 800k events when nothing has changed.

## Layout engine

```python
def _grid(n_obs):
    table = {1:(1,1), 2:(1,2), 3:(2,2), 4:(2,2), 5:(2,3), 6:(2,3), ...}
```
A lookup table maps the number of requested observables to the optimal (rows,
cols) grid. For example, 3 observables uses a 2×2 grid (4 cells) with one
cell hidden — two on top, one on the bottom, which fills the canvas better
than a squished 1×3 layout.

---

# FILE 10 — Analyze.py (Analyse.py)

## Purpose

A comprehensive quantitative analysis suite with 8 distinct analyses
(A through H) for a single SPS+DPS pair.

## Section A: Ratio plots

```python
hs = _hist(sps[key], edges); ns = hs.sum() * bw
hd = _hist(dps[key], edges); nd = hd.sum() * bw
ps = hs / (ns + 1e-12)
pd = hd / (nd + 1e-12)
ratio = np.divide(ps, pd, out=np.full_like(ps, np.nan, dtype=float), where=(pd > 0))
err_r = ratio * np.sqrt((err_s / (ps + 1e-12))**2 + (err_d / (pd + 1e-12))**2)
```
Normalises both histograms to probability densities, then takes the ratio.
The error on the ratio uses standard error propagation:
  σ(A/B) = (A/B) √((σ_A/A)² + (σ_B/B)²)
where σ_A = √N_A / (N_total × bw) and similarly for B. The `+ 1e-12`
prevents division by zero when ps or pd is exactly zero.

## Section B: Correlation matrices

```python
r, _ = pearsonr(data[ki], data[kj])
mat[i, j] = mat[j, i] = r
```
`scipy.stats.pearsonr` computes the Pearson correlation coefficient:
  r = Cov(X,Y) / (σ_X σ_Y)
r=1 means perfect positive linear correlation, r=0 means no linear correlation.
**Physics**: In SPS, the W and Z come from the same hard scatter, so their
kinematics are correlated (e.g., r(η_W, η_Z) should be positive). In DPS,
they are independent, so r ≈ 0 for all pairs. The difference Δr = r_SPS − r_DPS
quantifies this.

## Section C: 2D density plots

```python
h, xedg, yedg = np.histogram2d(x, y, bins=nbins, range=[[xlo,xhi],[ylo,yhi]])
im = ax.pcolormesh(xedg, yedg, h.T, norm=LogNorm(vmin=1), cmap="viridis")
```
`np.histogram2d` bins data in two dimensions simultaneously. `h.T` transposes
because pcolormesh expects the first index to vary along x, but histogram2d
returns the first index varying along the first array (x here). `LogNorm`
uses a logarithmic colour scale, necessary because events cluster strongly at
low pT and there are rare events at high pT.

## Section E: Shape parameters

```python
from scipy.stats import skew, kurtosis
skew(sv)       # Fisher skewness
kurtosis(sv)   # excess kurtosis (0 for Gaussian)
```
**Physics**: Skewness measures asymmetry — all pT distributions have positive
skewness (long tail toward high pT). Excess kurtosis measures whether the
distribution has heavier tails than a Gaussian. Comparing SPS vs DPS for
these parameters quantifies differences beyond just mean and RMS.

## Section F: KS vs N

```python
for n in ns_sorted:
    s = load(sp); d = load(dp)
    for key, *_, tag in OBSERVABLES:
        ks, _ = ks_2samp(s[key], d[key])
        ks_data[tag].append(ks)
```
Loads EACH N sample and computes the KS statistic at each N. If KS is stable
(plateau), the separation between SPS and DPS is robust — it is not a
finite-statistics artefact. If KS is still rising at the largest N, more
events are needed. This is the convergence criterion for the physics result.

## Section H: Chi-squared test

```python
chi2_val = (((hs[mask] - hd_sc[mask])**2) / hd_sc[mask]).sum()
chi2_red = chi2_val / ndof
```
The reduced chi-squared χ²/ndof compares histograms bin-by-bin:
  χ² = Σ_bins (N_SPS − N_DPS_scaled)² / N_DPS_scaled
`χ²/ndof ≫ 1` means the distributions are statistically incompatible (SPS
and DPS differ significantly). `χ²/ndof ≈ 1` means they are consistent
within statistical fluctuations. The `mask` requires ≥5 counts per bin,
which is the validity condition for Poisson statistics to approximate a
Gaussian in the chi-squared formula.

---

# FILE 11 — Makefile

## Purpose

GNU Make automates the build-and-run pipeline. It tracks file dependencies
so that only what needs to be rebuilt (because source files changed) is
recompiled.

## Variable extraction

```makefile
NZ := $(shell grep -E 'constexpr\s+int\s+N_EVENTS_Z\s*=' Settings.h \
          | grep -oE '[0-9]+' | tail -1)
```
`:=` is an immediate assignment (evaluated once at parse time). `$(shell ...)`
runs a shell command and captures its output. The grep regex matches the exact
line in Settings.h, extracts all digit sequences, and takes the last one.
This happens at make-time, before any recipe runs. The result is a number
like "30000" stored in `$(NZ)`.

## Stamped binary names

```makefile
BIN_Z = $(BUILD)/GenerateZ_$(NZ)
```
The binary name encodes N: `build/GenerateZ_30000`. If you change N in
Settings.h, the new binary has a different name and must be recompiled.
The old binary is left in build/ (it will be cleaned by `make clean`).
This prevents accidentally running a binary compiled for N=10000 when you
meant N=800000.

## Dependency rules

```makefile
$(BIN_Z): $(SRC)/GenerateZ.cc Settings.h | $(BUILD)
    $(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $@ $(LDFLAGS)
```
This says: to build `$(BIN_Z)`, the prerequisites are `GenerateZ.cc` AND
`Settings.h`. If either changes, the binary is rebuilt. `$<` is the first
prerequisite (GenerateZ.cc). `$@` is the target (the binary path). The `|`
separator means `$(BUILD)` is an "order-only" prerequisite — the build
directory must exist but its timestamp doesn't trigger a rebuild.

## Parallel execution

```makefile
@$(RUN) $(BIN_Z) > $(LOG_Z) 2>&1 & \
 $(RUN) $(BIN_W) > $(LOG_W) 2>&1 & \
 wait
```
`&` at the end of a shell command runs it in the background. Both GenerateZ
and GenerateW run simultaneously. `wait` blocks until all background jobs
finish. Since GenerateZ and GenerateW are independent (neither reads the
other's output), running them in parallel halves the wall-clock time.

Similarly, SPS, DPS, and DPSReal all run in parallel in step 2:
```makefile
@$(RUN) $(BIN_SPS)     > ... 2>&1 & \
 $(RUN) $(BIN_DPS)     > ... 2>&1 & \
 $(RUN) $(BIN_DPSREAL) > ... 2>&1 & \
 wait
```
DPS (event mixing) can run in parallel with SPS and DPSReal because it only
reads W.root and Z.root (already produced in step 1). DPSReal is self-
contained (runs PYTHIA itself). SPS is also self-contained.

## PHONY targets

```makefile
.PHONY: all build run compare ...
```
A phony target is one that doesn't correspond to a file. Without this, if a
file named "all" existed, Make would think the "all" target was already built
and do nothing. Declaring it `.PHONY` prevents this.

## LD_LIBRARY_PATH

```makefile
RPATH = $(PYTHIA8)/lib:$(shell root-config --libdir)
RUN   = LD_LIBRARY_PATH=$(RPATH)
```
Shared libraries (libpythia8.so, libROOT.so) are searched at runtime in
the directories listed in `LD_LIBRARY_PATH`. Prepending `LD_LIBRARY_PATH=...`
to the command sets it only for that single command, not globally. This
avoids polluting the user's environment.

---

# FILE 12 — settings_editor.py

## Purpose

Provides a safe, menu-driven way to modify Settings.h without knowing C++
or risking syntax errors.

## Regex parsing

```python
m = re.search(rf"constexpr\s+double\s+{token}\s*=\s*([0-9.eE+\-]+)", content)
```
`rf"..."` is a raw f-string: `f` enables variable interpolation (`{token}`),
`r` suppresses backslash interpretation. `\s+` matches one or more whitespace
characters. `([0-9.eE+\-]+)` captures the number including scientific notation
(e.g., "1.36e4"). This regex is specific enough to match the exact constexpr
declaration without accidentally matching a comment.

## In-place modification

```python
new_content, n = re.subn(pattern, repl, content)
if n == 0:
    print(f"  WARNING: token '{token}' not found")
```
`re.subn` returns the modified string AND a count of substitutions made.
If n=0, the pattern didn't match — possibly because the token name changed.
The warning alerts the user without crashing.

## Backup mechanism

```python
shutil.copy2(SETTINGS_FILE, backup)
with open(SETTINGS_FILE, "w") as f:
    f.write(content)
```
`copy2` preserves file metadata (timestamps). It creates a backup before
writing. The Makefile reads timestamps to decide whether to rebuild — the
backup having the old timestamp prevents unnecessary rebuilds.

---

# FILE 13 — run_all.sh

## Purpose

An alternative to `make all` for environments where Make is not available or
where more detailed logging control is needed.

## Pipeline structure

```bash
set -euo pipefail
```
- `-e`: exit immediately on any command failure
- `-u`: treat unset variables as errors
- `-o pipefail`: a pipeline (cmd1 | cmd2) fails if any stage fails, not just
  the last one. Together these three options make the script fail fast and
  loudly rather than silently continuing after an error.

## Parallel execution with explicit PIDs

```bash
"$BUILD/GenerateZ_${NZ}" > "$LOG_Z" 2>&1 &  PID_Z=$!
"$BUILD/GenerateW_${NW}" > "$LOG_W" 2>&1 &  PID_W=$!

wait $PID_Z && echo "done" || { echo "ERROR"; exit 1; }
wait $PID_W && echo "done" || { echo "ERROR"; exit 1; }
```
`$!` captures the process ID (PID) of the most recently backgrounded job.
`wait $PID` waits for that specific process and returns its exit code.
The `&&` chain runs echo on success; the `||` chain runs the error handler
on failure. `{ ...; exit 1; }` is a compound command that exits the entire
script on failure. This is more explicit than Makefile's `& wait` because
it identifies which specific generator failed.

## N extraction

```bash
extract() {
    grep -E "constexpr\s+int\s+${1}\s*=" Settings.h \
        | grep -oE '[0-9]+' | tail -1
}
NZ=$(extract 'N_EVENTS_Z')
```
Bash function using the same grep strategy as the Makefile. `${1}` is the
function's first argument. `tail -1` takes the last matching number — in case
the grep line matches multiple numbers (e.g., in a comment), the last is the
actual constexpr value which appears at the end of the line.

---

## CROSS-FILE CONTRACT SUMMARY

The files communicate through two protocols:

**Protocol 1 — ROOT file naming**:
GenerateZ writes `root/Z_N.root` with tree "ZTree".
GenerateDPS reads `root/W_N.root` ("WTree") and `root/Z_N.root` ("ZTree").
Python scripts read `root/SPS_N.root`, `root/DPS_N.root`, `root/DPSReal_N.root`
all with tree name "Events".

**Protocol 2 — ROOT tree branch names**:
All "Events" trees (SPS, DPS, DPSReal) have IDENTICAL branch names:
W_pt, W_eta, W_phi, W_mass, W_rapidity, W_px, W_py, W_pz, W_E, W_charge,
Z_pt, Z_eta, Z_phi, Z_mass, Z_rapidity, Z_px, Z_py, Z_pz, Z_E,
dEta_WZ, dPhi_WZ, dY_WZ, M_WZ, pT_WZ, sT, weight, nJets30.

This is the engineering reason why GenerateDPS.cc, GenerateDPSReal.cc, and
GenerateSPS.cc all produce the exact same branch list: any Python script that
reads one tree can read any of the three. Compare.py, CompareDPSModels.py,
Convergence.py, and Analyze.py all use `uproot.open(path)["Events"].arrays()`
without branching on which file type they are reading.

**Protocol 3 — N encoding**:
Settings.h defines N. The Makefile reads it at make-time. Binaries are named
with N. ROOT files are named with N. Python scripts discover files by globbing
for `SPS_*.root` and extracting N from the filename. This ensures all pieces
of a given run (binary, log, ROOT file, plots) are labelled with the same N.
