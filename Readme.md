# ZW_DPS — Full Project Code Explanation

This document walks through **every file** in the project, explaining what each
piece of code does, the physics reasoning behind it, why it's written the way
it is, and how the files connect to each other. It is organized file-by-file,
and within each file, block-by-block in the order the code actually runs.

Repetitive boilerplate (e.g. registering 60 near-identical ROOT branches) is
explained once with the pattern, rather than re-explained 60 times — but every
*distinct* line of logic is covered.

---

## 0. The physics question this whole project answers

You are comparing two different mechanisms by which a single proton-proton
collision at the LHC can produce a **W boson and a Z boson simultaneously**:

- **SPS (Single Parton Scattering)**: one quark from proton A and one
  antiquark from proton B annihilate in a single hard interaction,
  `qq̄' → WZ`. Both bosons come from the *same* vertex — the same Feynman
  diagram, the same momentum transfer. They are kinematically tied together
  from the moment they're created: momentum conservation at that one vertex
  forces them to be (nearly) back-to-back in azimuth, correlates their
  rapidities (they share one longitudinal boost), and constrains their
  transverse momenta to nearly cancel.

- **DPS (Double Parton Scattering)**: **two independent hard scatters**
  happen within the *same* proton-proton collision. One parton pair produces
  the Z; a separate, unrelated parton pair produces the W. There is no
  Feynman diagram connecting them. The only thing they share is the
  underlying event: the same two protons, the same beam remnants, the same
  soft QCD activity (multi-parton interactions, colour reconnection).

The project's whole point is to generate both samples with PYTHIA 8, compute
a battery of kinematic observables for each, and quantify how distinguishable
the two production mechanisms are — plus study how much Monte Carlo
statistics (how many events, `N`) is needed before those differences are
trustworthy rather than noise.

---

## 1. `Settings.h` — shared constants and helper functions

This is a C++ header (not a `.cc` file — it has no `main()`), `#include`d by
every generator (`GenerateSPS.cc`, `GenerateDPS.cc`). Anything that must be
*identical* across both generators (beam energy, decay forcing, the `Δφ`
folding convention) lives here exactly once, so the two samples are never
accidentally generated under different physics assumptions.

```cpp
#pragma once
```
Standard include-guard: ensures the header's contents are only inserted once
per translation unit, even if (indirectly) included multiple times.

```cpp
constexpr double ECM_GEV = 13600.;
```
`constexpr` means this is a compile-time constant — the compiler can inline
its value everywhere it's used, and it can't accidentally be reassigned at
runtime. `13600.` GeV = 13.6 TeV, the current LHC Run 3 centre-of-mass
collision energy. This single number is what "how energetic is the collision"
means for every generated event.

```cpp
constexpr double Z_MASS_MIN = 75.;
constexpr double Z_MASS_MAX = 120.;
```
The invariant-mass window PYTHIA is told to restrict the *first hard
process* to, when that process is `qq̄→Z/γ*`. Physically: at low invariant
mass, the intermediate boson is dominated by an off-shell photon (`γ*`)
rather than a real Z — the amplitude blows up as the mass → 0 (photon
propagator `1/q²`). Restricting to 75–120 GeV keeps the sample dominated by
genuine on-shell Z bosons (`M_Z ≈ 91.2` GeV sits comfortably inside this
window) and removes the low-mass Drell-Yan continuum, which isn't the physics
this project is studying.

```cpp
constexpr int N_EVENTS_SPS = 800000;
constexpr int N_EVENTS_DPS = 800000;
```
How many events each generator produces *by default* — overridable from the
command line (`./GenerateSPS 50000`) or via `settings_editor.py`. These are
plain `int`, not `long`, because 800,000 comfortably fits in 32 bits (max
~2.1 billion) and every other piece of code that reads this value (bash
regex, Python parsing) expects a short decimal integer.

```cpp
constexpr bool ENABLE_MPI = true;
```
Multi-Parton Interactions: whether, on top of the one hard process you
explicitly asked PYTHIA to generate, PYTHIA also simulates additional softer
parton-parton scatters happening in the same collision (the "underlying
event"). This is turned on because a real pp collision always has MPI
activity — turning it off gives a cleaner but less realistic sample, useful
only for fast debugging (hence the comment).

```cpp
constexpr int ID_Z = 23;
constexpr int ID_WPLUS = 24;
constexpr int ID_WMINUS = -24;
```
[PDG particle IDs](https://pdg.lbl.gov/) — the standard integer codes every
particle physics tool uses to identify particle species. 23 = Z⁰, 24 = W⁺,
−24 = W⁻ (particle/antiparticle pairs get opposite-sign IDs in the PDG
scheme). These are used throughout the generators to search PYTHIA's event
record for "the Z" or "the W" by matching `particle.id()` against these
constants, instead of hard-coding magic numbers scattered through the code.

```cpp
#include "Pythia8/Pythia.h"
#include <cmath>
#include <string>
```
Brought in *after* the constants (not before) because the helper functions
below need `Pythia8::Pythia`, `std::string`, and `M_PI`/`fabs` from
`<cmath>`.

```cpp
inline void setQuiet(Pythia8::Pythia& p) {
    p.readString("Print:quiet      = on");
    p.readString("Next:numberCount = 0");
}
```
`inline` on a free function defined in a header means: if this header is
included in multiple `.cc` files that both get compiled into the same
program, the linker won't complain about "multiple definition" errors — each
translation unit gets its own copy, and the linker is told that's expected
and fine, since they're identical.

`pythia.readString(...)` is PYTHIA 8's universal configuration interface:
every tunable setting in PYTHIA is set by feeding it a string of the form
`"Setting:Name = value"`, which PYTHIA parses at runtime and looks up against
its internal settings database (this is *why* a typo like `ffbar2Z` instead
of `ffbar2gmZ` causes a runtime "setting not found" error rather than a
compile error — the string is opaque to the C++ compiler).

- `Print:quiet = on` suppresses PYTHIA's normal per-event and per-run
  console spam (parameter listings, cross-section tables mid-run).
- `Next:numberCount = 0` disables the "now processing event N" progress
  print that would otherwise appear periodically — useful for a script but
  noisy for a batch job generating 800,000 events.

```cpp
inline void setBeam(Pythia8::Pythia& p) {
    p.readString("Beams:idA = 2212");
    p.readString("Beams:idB = 2212");
    p.readString("Beams:eCM = " + std::to_string(ECM_GEV));
}
```
2212 is the PDG ID for a proton. Setting both beam A and beam B to 2212
configures a proton-proton collider (as opposed to e.g. proton-antiproton,
which would use −2212 for one beam). `Beams:eCM` sets the total
centre-of-mass energy of the collision — `std::to_string(ECM_GEV)` converts
the `double` constant into the string PYTHIA's settings parser expects
(PYTHIA's config strings are always text, regardless of the underlying
setting's C++ type).

```cpp
inline void setMPI(Pythia8::Pythia& p) {
    p.readString(std::string("PartonLevel:MPI = ") +
                (ENABLE_MPI ? "on" : "off"));
}
```
A ternary expression picks the string `"on"` or `"off"` based on the
compile-time `ENABLE_MPI` flag, concatenated onto the setting name. This is
the single choke point that turns MPI on/off consistently for *every*
generator that calls `setMPI()` — changing `ENABLE_MPI` in one place changes
it everywhere.

```cpp
inline void setZDecayMuMu(Pythia8::Pythia& p) {
    p.readString("23:onMode    = off");
    p.readString("23:onIfMatch = 13 -13");
}
```
This is PYTHIA's particle-decay-table syntax: `"23:onMode = off"` means "for
particle ID 23 (the Z), turn off *all* its default decay channels first."
`"23:onIfMatch = 13 -13"` then turns *on* specifically the channel whose
daughter list matches `{13, -13}` — PDG ID 13 is μ⁻, −13 is μ⁺. So this
forces every Z in the sample to decay to μ⁺μ⁻, and *only* that channel
(never to e⁺e⁻, τ⁺τ⁻, hadrons, invisibly to neutrinos, etc.).

**Why force this at all?** Two reasons. First, muons are the cleanest final
state to reconstruct kinematics from (no missing energy like taus/neutrinos,
no jet-clustering ambiguity like hadronic decays). Second, and just as
important: forcing a *single* decay channel means every event in the sample
has the exact same final-state topology, so every observable computed
downstream (angles between decay leptons, etc.) is directly comparable
event-to-event without having to handle a mix of different decay types.

```cpp
inline void setWDecayMuNu(Pythia8::Pythia& p) {
    p.readString("24:onMode      = off");
    p.readString("24:onIfMatch   = 13 -14");   // W+
    p.readString("-24:onMode     = off");
    p.readString("-24:onIfMatch  = -13 14");   // W-
}
```
Same idea, but note the W needs **two** separate blocks — one for W⁺ (PDG
24) and one for W⁻ (PDG −24) — because they are genuinely different
particles in PYTHIA's decay table (unlike the Z, which is its own
antiparticle). Charge conservation dictates which daughters are physically
allowed: W⁺ (charge +1) can only decay to μ⁺ (id −13, charge +1) plus ν_μ
(id +14, charge 0); W⁻ (charge −1) can only decay to μ⁻ (id +13) plus ν̄_μ
(id −14). `onIfMatch` takes the daughter-ID list for the channel to enable,
matched against that specific particle's own decay table — so `"24:onIfMatch
= 13 -14"` and `"-24:onIfMatch = -13 14"` only take effect if each string
matches a channel that's actually charge-consistent with that W. If you ever
need to edit this block, the safe way to check it is against PYTHIA's own
decay-table listing (`pythia.particleData.list(24)`) rather than re-deriving
the signs by hand — getting the pairing wrong wouldn't cause a compile error,
it would just silently leave that channel disabled, which is the kind of
mistake this project deliberately avoids by using named helper functions
like `setWDecayMuNu()` instead of scattering raw `readString` calls through
every generator file.

```cpp
inline double dPhi(double phi1, double phi2) {
    double d = std::fabs(phi1 - phi2);
    return (d > M_PI) ? 2.0 * M_PI - d : d;
}
```
Azimuthal angle `φ` is defined on a circle, conventionally reported in
`(−π, π]`. A naive `|φ1 − φ2|` can therefore return a value up to `2π`, which
double-counts the same physical separation (e.g. φ1=−179°, φ2=179° are only
2° apart on the circle, but `|φ1−φ2|` naively gives 358°). This function
**folds** the raw difference into `[0, π]`: if the naive difference exceeds
π, the *true* angular separation is `2π` minus that difference (going the
short way around the circle instead of the long way). This is the standard
definition used throughout collider physics for angular separation, and
every `dPhi(...)` observable in every generator and analysis script in this
project calls exactly this one function, so the folding convention is
identical everywhere.

---

## 2. `GenerateSPS.cc` — the SPS generator

### Header comment and includes
```cpp
#include "Settings.h"
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
using namespace Pythia8;
using namespace std;
```
`TFile`/`TTree` are ROOT's file-format and tabular-data classes — ROOT is the
standard HEP data-analysis framework, and its `.root` files (structured,
compressed, column-oriented binary files) are what every downstream Python
script (`Compare.py` etc., via `uproot`) reads. `using namespace Pythia8;`
and `using namespace std;` let the code write `Pythia`, `Vec4`, `string`
instead of `Pythia8::Pythia`, `Pythia8::Vec4`, `std::string` everywhere —
acceptable in a small, single-purpose `.cc` file (would be bad practice in a
shared header, which is why `Settings.h` doesn't do this).

### `main()` setup
```cpp
int nEvents = N_EVENTS_SPS;
if (argc > 1) nEvents = atoi(argv[1]);
```
Default event count comes from the compiled-in constant; but if the binary
is invoked with a command-line argument (`./GenerateSPS_800000 50000`), that
overrides it. `atoi` converts the C-string argument to an `int` (no error
checking — malformed input silently becomes 0, which is acceptable here
since this is an internal tool, not a public-facing one).

```cpp
string outPath = "root/SPS_" + to_string(nEvents) + ".root";
```
The output filename is **stamped with the event count** (`SPS_800000.root`,
not just `SPS.root`). This is deliberate: running the generator again with a
different `N` produces a *new* file instead of silently overwriting the old
one, which is what makes the "convergence study" (comparing the same
observable across many different `N`) possible later — every historical run
survives on disk under its own name.

### PYTHIA configuration — the actual physics process
```cpp
Pythia pythia;
setBeam(pythia);
setMPI(pythia);

pythia.readString("WeakDoubleBoson:all      = off");
pythia.readString("WeakDoubleBoson:ffbar2ZW = on");
pythia.readString("SecondHard:generate      = off");
```
`WeakDoubleBoson` is PYTHIA's process class for **diboson production from a
single hard vertex**. `:all = off` starts from a clean slate (nothing in
this class enabled), then `:ffbar2ZW = on` enables specifically the
`f f̄' → Z W` process — a fermion and antifermion of different flavour (e.g.
u and d̄) annihilating directly into a Z and a W. This is the Feynman diagram
that gives SPS its defining physics: **one vertex, one hard scatter,
momentum-conserving, both bosons correlated by construction**.

`SecondHard:generate = off` explicitly disables PYTHIA's separate machinery
for simulating a *second*, independent hard process in the same event — this
is the DPS-generating switch (see `GenerateDPS.cc` below), and it's
explicitly turned off here to guarantee the SPS sample never accidentally
contains a spurious second hard scatter.

```cpp
setZDecayMuMu(pythia);
setWDecayMuNu(pythia);
setQuiet(pythia);

if (!pythia.init()) { cerr << "ERROR: PYTHIA init failed.\n"; return 1; }
```
`pythia.init()` is where PYTHIA actually validates and locks in all the
settings strings fed to it above, builds its internal process/phase-space
machinery, and returns `false` if anything was inconsistent (e.g. a
misspelled setting name, an invalid mass window). Checking the return value
and bailing out with `return 1` (nonzero exit code) is what lets the
Makefile/shell scripts detect a failed run instead of silently producing a
garbage or empty ROOT file.

### ROOT tree and branch declarations
```cpp
TFile* fOut = TFile::Open(outPath.c_str(), "RECREATE");
TTree* tree = new TTree("Events", "SPS WZ production");
```
`"RECREATE"` mode: if a file already exists at `outPath`, overwrite it
completely (as opposed to `"UPDATE"`, which would append). The `TTree` named
`"Events"` is the actual table — one row per event, one column per branch —
and this exact name (`"Events"`) is what every Python reader
(`uproot.open(path)["Events"]`) expects to find.

The block of `double`/`int` variable declarations that follows —
`W_pt, W_eta, W_phi, ...` — are the **row buffer**: ROOT's `TTree::Branch`
API works by binding a branch name to the *memory address* of a C++
variable (`&W_pt`), and every time `tree->Fill()` is called, ROOT copies
whatever value currently sits at that address into the next row of that
branch's column. This is why the fill loop later always writes to the same
named variable and then calls `Fill()` — the address never changes, only
the value does.

The branches fall into clear groups (each declared with a comment above it
explaining its physics role):
- **W boson 4-vector + derived kinematics**: `W_pt, W_eta, W_phi, W_mass,
  W_rapidity` (the physically meaningful quantities) plus the raw Cartesian
  `W_px, W_py, W_pz, W_E` (needed downstream for any calculation that must
  combine two 4-vectors, like an invariant mass — you can't add two `pT`
  scalars and get a physically meaningful "combined pT", you must add the
  vector components).
- **Z boson**: identical structure, `Z_*`.
- **WZ-system observables**: `dEta_WZ, dPhi_WZ, dY_WZ, M_WZ, pT_WZ, sT` —
  quantities that describe the *relationship* between the two bosons, not
  either boson individually. These are computed once per event in the loop
  below (there's no PYTHIA setting that directly gives you "the angle
  between two already-decayed bosons" — you compute it yourself from their
  stored 4-vectors).
- **W decay products** (`lep_*`, `nu_*`) and **Z decay products**
  (`lep1_*`, `lep2_*`): the actual final-state particles a real detector
  would see, needed for the lepton-level `Δφ` observables requested later
  in the project (§4 below covers the physics reason these were added).
- **Cross-observables**: `dPhi_Z_lep, dPhi_Z_nu, dPhi_lep1_lep2,
  dPhi_lep_nu, dPhi_W_lep1, dPhi_W_lep2, dPhi_Z_lep1, dPhi_Z_lep2` — every
  physically distinct pairwise azimuthal separation among {W, Z, lepton,
  lepton1, lepton2, neutrino} that the project's observable list asks for.
  (Note `Δφ(W,ℓ)` and `Δφ(ℓ,W)` are the same physical angle — the code's
  comment explicitly says only one branch is stored per such pair, avoiding
  redundant/duplicate columns.)
- **`pT_imbalance`**: the vector-based relative-imbalance observable
  (explained fully in the "Physics fix #3" section below).
- **Event-level**: `weight` (PYTHIA's per-event weight — usually 1.0 for
  this kind of unweighted sample, but stored so weighted histograms are
  always possible without regenerating anything) and `nJets30` (a rough jet
  count, explained inline where it's computed).

```cpp
tree->Branch("W_pt", &W_pt);
...
```
Each `tree->Branch("name", &variable)` call is a one-line registration: "the
column called `name` in this tree is filled from whatever is currently in
`variable`." The block of ~60 such calls is mechanical but essential — this
*is* the ROOT file's schema, and every downstream Python script's
`d["W_pt"]`-style access depends on these exact string names matching.

### The event loop — where physics actually happens per event
```cpp
long nFilled = 0, nNoW = 0, nNoZ = 0, nNoLepNu = 0;

for (int iev = 0; iev < nEvents; ++iev) {
    if (!pythia.next()) continue;
```
`pythia.next()` generates one full event (hard process → parton shower →
hadronization → decays) and returns `false` if that particular attempt
failed internally (rare, but possible for numerically difficult phase-space
points) — `continue` simply skips to the next loop iteration without
incrementing any of the counters, silently discarding the failed attempt
(it doesn't count against `nEvents` "requested", but it also isn't logged
as a distinct failure category, since PYTHIA already retries internally
before giving up).

```cpp
int iZ = -1;
for (int i = 0; i < pythia.event.size(); ++i)
    if (pythia.event[i].id() == ID_Z) iZ = i;
if (iZ < 0) { ++nNoZ; continue; }
```
`pythia.event` is PYTHIA's full particle history for this one event — not
just the final-state particles, but every intermediate copy (a particle can
appear multiple times as it radiates, recoils, etc. through the shower).
This loop deliberately does **not** `break` on the first match — it keeps
overwriting `iZ` every time it finds another Z, so after the loop, `iZ`
points at the **last** Z entry in the record. This is intentional and
important: PYTHIA stores intermediate/pre-shower copies of a resonance
*before* its final, physical, post-recoil copy, and the *last* occurrence is
the one whose kinematics you actually want to measure. Searching for the
first match would silently give you a slightly-wrong, pre-recoil Z
four-momentum. If no Z is found at all (`iZ` stays `-1`), the event is
uselessly incomplete for this analysis and is skipped, counted under
`nNoZ`.

```cpp
int iW = -1;
W_charge = 0;
for (int i = 0; i < pythia.event.size(); ++i) {
    int id = pythia.event[i].id();
    if (id == ID_WPLUS || id == ID_WMINUS) {
        iW = i;
        W_charge = (id == ID_WPLUS) ? +1 : -1;
    }
}
if (iW < 0) { ++nNoW; continue; }
```
Same "last occurrence wins" logic, but for *either* W charge — since SPS's
`ffbar2ZW` process produces a W of either sign depending on which specific
quark flavours annihilated, `W_charge` is recorded per event rather than
assumed fixed.

```cpp
int iLep = -1, iNu = -1;
for (int d : pythia.event[iW].daughterList()) {
    int id = pythia.event[d].id();
    if (abs(id) == 13) iLep = d;
    if (abs(id) == 14) iNu  = d;
}
if (iLep < 0 || iNu < 0) { ++nNoLepNu; continue; }
```
`daughterList()` returns the indices of this particle's direct decay
products in the event record. Since `setWDecayMuNu` forced every W to decay
to exactly `{μ, ν_μ}`, this loop just identifies *which* daughter is the
muon (`|id|==13`, `abs()` because it could be μ⁻ or μ⁺ depending on W
charge) and which is the neutrino (`|id|==14`). If somehow neither is found
(shouldn't happen given the forced decay, but defensive coding), the event
is skipped.

```cpp
int iLep1 = -1, iLep2 = -1;
for (int d : pythia.event[iZ].daughterList()) {
    int id = pythia.event[d].id();
    if (id == 13)  iLep1 = d;
    if (id == -13) iLep2 = d;
}
if (iLep1 < 0 || iLep2 < 0) { ++nNoLepNu; continue; }
```
Same idea for the Z's two muon daughters — but here the signs matter and are
used to distinguish them: `lep1` is always specifically the μ⁻ (id=+13,
note PDG convention: the particle μ⁻ has *positive* PDG id 13, μ⁺ has −13 —
this is a common point of confusion, particle/antiparticle sign convention
in the PDG scheme doesn't track electric charge sign directly for leptons),
`lep2` is always the μ⁺. This gives a consistent, unambiguous labelling
(`lep1` = μ⁻ always, not "whichever muon PYTHIA listed first") so that
`dPhi_Z_lep1` means the same physical angle in every single event.

```cpp
nJets30 = 0;
for (int i = 0; i < pythia.event.size(); ++i) {
    const Particle& p = pythia.event[i];
    if (!p.isFinal()) continue;
    if (abs(p.id()) > 8 && p.id() != 21) continue;
    if (p.pT() > 30.) ++nJets30;
}
```
A deliberately crude jet-multiplicity proxy, explicitly documented as such.
`isFinal()` restricts to actual final-state particles (not intermediate
shower copies). `abs(id) > 8 && id != 21` keeps only quarks (PDG IDs 1–8)
and gluons (21) — i.e., only counts partons, not leptons/photons/etc. — with
`pT() > 30.` GeV as an arbitrary-but-standard-ish threshold for "hard
enough to plausibly become a reconstructed jet." This is explicitly *not* a
real jet algorithm (no clustering, no merging of nearby partons into a
single jet) — it's a fast approximation, good enough to distinguish
"basically no extra radiation" from "clearly some hard extra activity" event
by event, which is all this project uses it for.

```cpp
const Particle& pW = pythia.event[iW];
W_pt = pW.pT(); W_eta = pW.eta(); W_phi = pW.phi();
W_mass = pW.m(); W_rapidity = pW.y();
W_px = pW.px(); W_py = pW.py(); W_pz = pW.pz(); W_E = pW.e();
```
Once the right particle index is known, PYTHIA's `Particle` class exposes
all these kinematic quantities as simple accessor methods — `pT()` is
`sqrt(px²+py²)`, `eta()` is pseudorapidity `-ln(tan(θ/2))`, `y()` is true
rapidity `½ln((E+pz)/(E−pz))`, etc. — all computed internally by PYTHIA from
the particle's stored 4-momentum, so this code doesn't need to derive any of
these formulas itself.

```cpp
dEta_WZ = fabs(W_eta - Z_eta);
dPhi_WZ = dPhi(W_phi, Z_phi);
dY_WZ   = fabs(W_rapidity - Z_rapidity);
```
The three "gap" observables between the two bosons: pseudorapidity gap,
folded azimuthal gap (using the shared `dPhi()` helper from `Settings.h`),
rapidity gap. Note these are stored as **absolute values** (`fabs(...)`) —
this branch always holds `|Δη|`/`|Δy|`, never the signed difference. (The
Python analysis layer separately derives *signed* versions from the raw
`W_eta`/`Z_eta` branches when a signed, negative-axis plot is wanted — see
§Compare.py — rather than this generator storing both variants.)

```cpp
Vec4 pWZ = pW.p() + pZ.p();
M_WZ  = pWZ.mCalc();
pT_WZ = pWZ.pT();
sT    = W_pt + Z_pt;
```
`Vec4` is PYTHIA's four-vector class; `pW.p()` returns the W's full 4-vector
`(E, px, py, pz)`, and `+` between two `Vec4`s is genuine 4-vector addition
(component-wise). `pWZ.mCalc()` computes the invariant mass of that *summed*
4-vector, `√(E²−p²)` — the invariant mass of the combined WZ system, which
is not the same thing as `W_mass + Z_mass` (in general, greater than or
equal to it, with equality only in the special case the two bosons are
exactly collinear). `pWZ.pT()` is the transverse momentum of that same
summed vector — physically, "how much net transverse momentum does the WZ
system as a whole carry" (nonzero only from initial-state radiation recoil
in SPS, since two genuinely back-to-back objects sum to zero net pT). `sT`
by contrast is a **scalar** sum of the two individual pT magnitudes — a
different, simpler observable that doesn't care about the *relative
direction* of the two bosons, only their individual hardness added together.

```cpp
const Particle& pLep = pythia.event[iLep];
lep_pt = pLep.pT(); ...
const Particle& pNu = pythia.event[iNu];
nu_pt = pNu.pT(); ...
dPhi_Z_lep = dPhi(Z_phi, lep_phi);
dPhi_Z_nu  = dPhi(Z_phi, nu_phi);
```
Same accessor pattern applied to the W's two decay daughters, then the two
new lepton-level `Δφ` observables computed directly from the just-filled
`phi` values.

```cpp
const Particle& pLep1 = pythia.event[iLep1]; ...
const Particle& pLep2 = pythia.event[iLep2]; ...

dPhi_lep1_lep2 = dPhi(lep1_phi, lep2_phi);
dPhi_lep_nu    = dPhi(lep_phi, nu_phi);
dPhi_W_lep1    = dPhi(W_phi, lep1_phi);
dPhi_W_lep2    = dPhi(W_phi, lep2_phi);
dPhi_Z_lep1    = dPhi(Z_phi, lep1_phi);
dPhi_Z_lep2    = dPhi(Z_phi, lep2_phi);
```
Every remaining pairwise azimuthal separation the project's observable list
asks for, computed the same way, once each event's daughters are known.

```cpp
{
    double sumPx = W_px + Z_px;
    double sumPy = W_py + Z_py;
    double denom = W_pt + Z_pt;
    pT_imbalance = (denom > 0.)
        ? sqrt(sumPx*sumPx + sumPy*sumPy) / denom : 0.;
}
```
This is the corrected, **vector** definition of relative pT imbalance (see
the dedicated physics discussion at the end of this document, §Physics
fix 3). The enclosing `{ }` block is just a local scope so `sumPx`,
`sumPy`, `denom` don't pollute the surrounding function's namespace — they're
genuinely only needed for this one calculation. The ternary guards against
division by zero in the (essentially impossible, but defensively handled)
case both bosons have exactly zero pT.

```cpp
weight = pythia.info.weight();
tree->Fill();
++nFilled;
```
`pythia.info.weight()` retrieves this event's Monte Carlo weight (for this
unweighted leading-order process, effectively always 1.0, but read from
PYTHIA rather than hard-coded so the code is correct even if a weighted
process is ever swapped in later). `tree->Fill()` is the moment all ~60
currently-set branch variables get copied into a new row of the tree.
`nFilled` tracks how many events actually made it all the way through every
selection check above.

### Wrap-up
```cpp
fOut->cd();
tree->Write();
fOut->Close();
pythia.stat();
```
`fOut->cd()` makes sure the current ROOT "directory" context is this output
file (relevant if ROOT has multiple files open, to avoid writing the tree
into the wrong one). `tree->Write()` flushes the tree's data to disk.
`fOut->Close()` closes the file, finalizing it. `pythia.stat()` prints
PYTHIA's own end-of-run statistics — cross-sections, efficiency, error
counts — to the console, independent of anything this project's own code
tracked.

```cpp
cout << "\n--- GenerateSPS summary ---\n";
cout << "  Events requested  : " << nEvents << "\n";
cout << "  Events saved      : " << nFilled << "\n";
...
```
A human-readable summary so you can sanity-check, at a glance, how many of
the requested events actually survived all the selection cuts (a large gap
between requested and saved, or a large `nNoLepNu` count, would be a red
flag worth investigating).

---

## 3. `GenerateDPS.cc` — the DPS generator

Structurally this file is extremely similar to `GenerateSPS.cc` — same
includes, same branch-declaration pattern, same event loop skeleton, same
lepton-finding logic, same summary printout. Rather than re-explain every
identical line, this section focuses on **everything that's genuinely
different**, which is where the actual physics content of "DPS vs SPS" lives.

### Why DPS needs a fundamentally different PYTHIA configuration

```cpp
pythia.readString("WeakSingleBoson:all             = off");
pythia.readString("WeakSingleBoson:ffbar2gmZ       = on");
pythia.readString("WeakBosonAndParton:qqbar2gmZg   = on");
pythia.readString("WeakBosonAndParton:qg2gmZq      = on");
pythia.readString("PhaseSpace:mHatMin = " + to_string(Z_MASS_MIN));
pythia.readString("PhaseSpace:mHatMax = " + to_string(Z_MASS_MAX));
```
This configures the **first hard process** as plain single-boson Z
production (`ffbar2gmZ` = `qq̄→Z/γ*`, plus the two associated
Z+parton channels `qqbar2gmZg`/`qg2gmZq`, which allow the Z to recoil
against a hard extra parton at the matrix-element level rather than only via
the parton shower). This is a completely ordinary, single-vertex process —
by itself indistinguishable from "just generate some Z bosons." The `Z_MASS_MIN`/`Z_MASS_MAX`
window plays the same role as in `GenerateSPS.cc`: keep this on-shell Z, not
low-mass `γ*`.

```cpp
pythia.readString("SecondHard:generate = on");
pythia.readString("SecondHard:SingleW  = on");
pythia.readString("SecondHard:WAndJet  = on");
```
This is the crux of the whole DPS implementation. `SecondHard:generate = on`
tells PYTHIA: *in addition* to the first hard process configured above,
generate a **second, statistically independent hard parton-parton scattering
within the same event** — this is PYTHIA's built-in machinery for exactly
the physics phenomenon this project studies. `SecondHard:SingleW` and
`SecondHard:WAndJet` restrict *what kind* of second process is allowed:
direct `qq̄'→W±` (Drell-Yan-like) and `W`+parton respectively. So the full
picture per event: **one parton pair inside this proton-proton collision
makes a Z, a completely separate parton pair inside the same collision makes
a W** — genuine DPS, not two bosons artificially glued together from
different events.

The header comment in the file explains *why* this particular
first-hard=Z / second-hard=W split was chosen rather than the reverse or
some other combination: PYTHIA's `SecondHard` machinery requires the first
and second hard processes to be drawn from *compatible* process classes, and
there is no built-in PYTHIA mode that mixes "Z" and "W" freely in either
slot — this specific combination (`WeakSingleBoson` first + `SecondHard:
SingleW/WAndJet` second) is the documented, supported way to get one Z and
one W from two independent hard scatters in the same PYTHIA 8 event.

### The "physics fix" block — restoring the expected `Δφ(W,Z)` enhancement

This is the largest and most carefully-commented block in the file, and it's
the direct response to an earlier finding in this project: a *naive*
`SecondHard` DPS configuration produces a **completely flat** `Δφ(W,Z)`
distribution, but the actual physical expectation (per the project's
supervisor) is a **mild upward enhancement as Δφ→π** — not as strong as
SPS's sharp back-to-back peak (which comes from momentum conservation at one
shared vertex, absent here), but a real, non-zero residual correlation.

**Why would a naive DPS configuration be *too* flat?** Because, as
configured with just the settings above, PYTHIA treats the two hard scatters
as essentially fully disconnected sub-events: independent colour flow,
independent local recoil. Real DPS in an actual proton-proton collision is
not quite that independent — both hard scatters are still happening *inside
the same two protons*, so:

1. **Shared beam remnants constrain both scatters.** Whatever momentum,
   colour, and flavour content is pulled out of proton A to make the Z
   changes what's left over in proton A's remnant for the W-producing
   scatter to draw from (and vice versa for proton B). This is a genuine
   physical linkage between the two scatters that a fully-independent
   simulation misses.
2. **Colour reconnection can act across the whole event.** The colour
   strings stretching from the Z-producing system to the beam remnants, and
   from the W-producing system to the (same) beam remnants, are not
   required to be blind to each other during hadronization.
3. **Soft radiation (ISR) recoil bookkeeping matters.** How the shower
   assigns recoil for radiated soft gluons affects each hard system's net
   transverse balance, and doing this locally (dipole-based) rather than
   globally changes how much of that soft activity's momentum leaks between
   the two systems.

None of these three effects create anything like SPS's *hard* momentum-
conservation constraint — they're soft, indirect, "the two systems share a
collision, not a vertex" effects — but they are real and expected to
produce a *mild* large-Δφ enhancement, which is exactly what's missing from
the naive flat baseline.

```cpp
pythia.readString("ColourReconnection:reconnect        = on");
```
Enables cross-system colour reconnection during hadronization: colour lines
from the Z system and the W system, both ultimately anchored in the same
two beam remnants, are allowed to reconnect with each other rather than
being forced to stay entirely separate. The removed/rejected alternative
setting `BeamRemnants:reconnectRange = 10.0` is explicitly documented as
**not valid** in this PYTHIA version (8.316) — the comment records that this
was tried, found to break `pythia.init()`, and removed, so a future editor
doesn't waste time reintroducing the same mistake.

```cpp
pythia.readString("SpaceShower:dipoleRecoil = on");
```
Switches initial-state-radiation recoil from PYTHIA's default global scheme
to a **dipole** scheme, where recoil from a radiated parton is assigned
locally to the specific colour dipole involved in that emission, rather than
smeared uniformly across the whole event. With two independent hard systems
present in the same event, this keeps each system's own hard kinematics
physically sensible while still allowing soft emissions near the shared
beam remnants to correlate the two systems' net transverse balance — which
is exactly the soft-correlation mechanism item 3 above describes.

**What was deliberately left out, and why**: the comment explicitly notes
that `MultipartonInteractions:allowRescatter` (a setting that lets the *same*
outgoing parton from one scatter re-scatter again) is intentionally **not**
enabled — the baseline DPS sample is meant to isolate the soft
beam-remnant/colour-reconnection correlation mechanism specifically, not mix
in a different (and stronger) correlation mechanism at the same time. This
is a controlled, single-variable physics change, documented as such.

### The rest of the file

The branch declarations, event loop, lepton-finding, jet-counting, `Δφ`
computations, and `pT_imbalance` calculation are **structurally identical**
to `GenerateSPS.cc` — same reasoning applies line-for-line (see §2 above).
The only field that's genuinely absent compared to a hypothetical "what if
DPS also needed something SPS doesn't" is nothing — the two generators now
produce ROOT files with an **identical branch schema**, which is exactly
what makes it possible for `Compare.py` to load an `SPS_*.root` and a
`DPS_*.root` and plot them on the same axes: every branch name means the
same physical quantity in both files.

---

## 4. `Compare.py` — SPS vs DPS comparison at a single N

### Purpose and where it sits in the pipeline
Reads exactly one `SPS_<N>.root` and one `DPS_<N>.root` (same `N`, produced
by the two generators above), and produces overlay plots, KS-test
statistics, and summary tables for a single fixed sample size.

### The observable catalogue
```python
ALL_OBSERVABLES = [
    ("dEta_signed", r"$\Delta\eta(W,Z)$", "", (-5, 5), False, "dEta", True),
    ("dPhi_WZ", r"$|\Delta\phi(W,Z)|$", "rad", (0, 3.1416), False, "dPhi", True),
    ...
]
```
This is a single **table-as-data-structure**: every observable the script
knows how to plot is one tuple of `(branch_or_derived_key, LaTeX_label,
axis_unit, (x_min, x_max), use_log_y, short_filename_tag, use_histogram)`.
Every function downstream (menu printing, plotting, file naming) reads this
one list rather than having any observable-specific `if/elif` chains — to
add or remove an observable from every plot type at once, you edit exactly
one line here.

Note `dEta_signed` and `dY_signed` (not `dEta_WZ`/`dY_WZ` directly) — these
are **not** ROOT branches; they're derived in Python (see `load()` below)
from the raw `W_eta`/`Z_eta` branches, specifically so the plotted quantity
can be *signed* (range −5 to 5) rather than the generator's stored `|Δη|`
(which is always ≥0 and would lose the "which boson is more forward"
information). The range is clipped to ±5 rather than the kinematically
possible wider range, matching the finite pseudorapidity acceptance of a
real detector rather than plotting generator-level tails no real experiment
could ever see.

`pT_imbalance`'s label — `r"$|\vec{p}_T^{W}\!+\!\vec{p}_T^{Z}|/(p_T^W\!+\!p_T^Z)$"`
— and its range `(0, 1)` reflect the corrected vector definition (§Physics
fix 3): this quantity is a ratio of a vector magnitude to a scalar sum, so
it's mathematically guaranteed to lie in `[0, 1]`, unlike the old
scalar-difference-based definition which could be negative.

The `dPhi_Z_lep`, `dPhi_Z_nu`, `dPhi_lep1_lep2`, `dPhi_lep_nu`,
`dPhi_W_lep1`, `dPhi_W_lep2`, `dPhi_Z_lep1`, `dPhi_Z_lep2` entries are the
new lepton-level `Δφ` observables (§Physics fix 2) — they read directly from
the correspondingly-named ROOT branches the generators now fill, with no
extra derivation needed in Python (the folding to `[0,π]` already happened
in C++ via the shared `dPhi()` helper).

### `AUTO_MODE` / `safe_input`
```python
AUTO_MODE = "--auto" in sys.argv or os.environ.get("ZW_DPS_AUTO") == "1"

def safe_input(prompt):
    global AUTO_MODE
    try:
        return input(prompt)
    except EOFError:
        print("  [no interactive input available -- switching to auto mode]")
        AUTO_MODE = True
        return ""
```
This solves a genuine, previously-encountered bug: checking
`sys.stdin.isatty()` to decide "should I show an interactive menu" is
unreliable, because several IDE "Run" consoles report `isatty()==False` even
though they *do* forward keystrokes — so the script was silently skipping
its own menu. The fix uses an **explicit** flag instead: only `--auto` on
the command line (used by the Makefile/`run_all.sh` pipeline) or the
`ZW_DPS_AUTO=1` environment variable triggers non-interactive mode. Plain
`python3 Compare.py`, from anywhere, always *tries* to prompt — and
`safe_input()` only falls back to auto-mode if `input()` genuinely raises
`EOFError` (truly no stdin to read from at all), so the script can never
hang forever waiting for input that will never come, but also never
silently skips a menu that a real user could have answered.

### Loading data: `load()`
```python
def load(path, tree="Events"):
    d = uproot.open(path)[tree].arrays(library="np")
    ...
    d["pT_imbalance"] = ... (only if not already a branch, else this is skipped/overridden by the stored one)
    d["dEta_signed"] = d["W_eta"] - d["Z_eta"]
    d["dY_signed"]   = d["W_rapidity"] - d["Z_rapidity"]
    return d
```
`uproot` is a pure-Python library for reading ROOT files *without* needing a
full ROOT installation with Python bindings — it parses ROOT's binary
format directly. `.arrays(library="np")` returns a dictionary mapping every
branch name to a NumPy array of that branch's values across all events —
so `d["W_pt"]` is a NumPy array of every event's `W_pt`, and computing a
derived quantity like `d["W_eta"] - d["Z_eta"]` is a single vectorized NumPy
subtraction across the *entire* dataset at once (not a Python-level loop
over events — this matters for performance with 800,000-row arrays). The
signed Δη/Δy derivation happens here specifically because the generator only
stores the folded `|...|` version — deriving the signed version in Python
avoids needing to regenerate any ROOT files just to get a differently-signed
view of data that's already fully present in the raw branches.

### Binning strategy: `_adaptive_edges` / `_log_edges` / `_choose_edges`
This is one of the more subtle pieces of the codebase, worth explaining
carefully because getting histogram binning wrong produces genuinely
misleading plots.

```python
def _adaptive_edges(pooled_vals, lo, hi, nbins):
    v = np.asarray(pooled_vals)
    v = v[(v >= lo) & (v <= hi)]
    if v.size < nbins * 2:
        return np.linspace(lo, hi, nbins + 1)
    edges = np.quantile(v, np.linspace(0, 1, nbins + 1))
    edges[0], edges[-1] = lo, hi
    edges = np.unique(edges)
    if edges.size < 3:
        return np.linspace(lo, hi, nbins + 1)
    return edges
```
Rather than equal-width bins, this builds bin edges from the **data's own
quantiles** — every bin contains (approximately) the same *number* of
events, which automatically means: many narrow bins where events are dense
(resolving a sharp peak well) and few wide bins where events are sparse
(without those bins being so numerous they're each just noise). `np.unique`
guards against duplicate edges (possible if many events share an identical
value, e.g. from a hard phase-space cut). The `if v.size < nbins*2` and
`if edges.size < 3` checks fall back to plain equal-width bins whenever
there isn't enough data for quantile binning to be meaningful — avoiding
degenerate/garbage bin edges on small samples.

```python
def _log_edges(lo, hi, nbins, floor=1.0):
    floor = min(floor, hi * 0.01) if hi > 0 else floor
    if lo <= 0:
        tail_edges = np.geomspace(max(floor, 1e-6), hi, nbins)
        return np.concatenate(([0.0], tail_edges))
    return np.geomspace(lo, hi, nbins + 1)
```
For a **steeply falling spectrum** (like `pT(W)`, `pT(Z)`, `M_WZ`) plotted
on a log-y axis, quantile binning is actually the *wrong* tool: since most
events sit at low pT, quantile bins get almost entirely used up describing
the first small slice of the range, leaving the long tail covered by only a
handful of enormous bins (visually blocky, chunky steps). The correct tool
here is **log-spaced** bins — `np.geomspace` produces edges that grow by a
constant *ratio* rather than a constant absolute step, giving smoothly
growing bin widths all the way out through the tail, resolving both the
peak and the tail properly. The `floor` handles the fact that log-spacing
is undefined starting exactly at 0 (a real, physical lower limit for a pT
spectrum) — a `[0, floor)` bin catches everything below the smallest
meaningful geometric step, and geometric spacing takes over from `floor`
onward.

```python
def _choose_edges(pooled_vals, lo, hi, nbins, logy):
    if logy:
        return _log_edges(lo, hi, nbins)
    return _adaptive_edges(pooled_vals, lo, hi, nbins)
```
The dispatcher: every observable's own `use_log_y` flag (from the catalogue)
determines which binning strategy it gets, automatically — falling
pT/mass spectra get log bins, bounded/peaked angular observables get
quantile bins. Nothing else in the codebase needs to know or care which
strategy was used; it's decided once, here, from data already in the
catalogue.

### `_norm_hist` / `_hist`
```python
def _norm_hist(vals, edges):
    counts, _ = np.histogram(vals, bins=edges)
    bw = np.diff(edges)
    total = counts.sum()
    density = counts / (total * bw + 1e-12)
    err = np.sqrt(counts) / (total * bw + 1e-12)
    cx = 0.5 * (edges[:-1] + edges[1:])
    return cx, density, err
```
Turns raw per-bin counts into a **normalized probability density** — dividing
by `total * bin_width` rather than just `total` is what makes the histogram
correctly integrate to 1 over its range *even with variable-width bins*
(each bin's height reflects density, not raw count, so a wide bin with many
events and a narrow bin with few events can still be visually compared on
the same y-axis meaningfully). `np.sqrt(counts)` is the **Poisson**
statistical uncertainty on a raw count (the standard assumption for
"how many independent random events landed in this bin"), propagated through
the same normalization to get an uncertainty on the *density*. The `1e-12`
in the denominator is a numerical-safety epsilon, preventing division by
exactly zero if a bin somehow has zero total width or zero events summed
(defensive, essentially never triggered in practice).

### Plot styling and structure
The `plt.rcParams.update({...})` block sets a single, consistent visual
style (serif font, specific tick direction/size, DPI, PDF font embedding via
`pdf.fonttype: 42` so text stays selectable/searchable in the output PDF
rather than being rasterized) applied automatically to *every* figure this
script creates — no per-plot styling code needed elsewhere.

`plot_panel(...)`, `make_process_single(...)`, `make_summary(...)`,
`make_eta_heatmap(...)` are the actual figure-drawing functions: given
already-loaded data and an observable's catalogue entry, draw the correct
axis (histogram + error bars via `ax.stairs`/`ax.errorbar`, or the
`(η_W,η_Z)` 2D joint-density heatmap via `ax.hist2d`), compute and annotate
the KS statistic (`scipy.stats.ks_2samp`, testing whether the SPS and DPS
samples for this observable could plausibly be drawn from the same
underlying distribution — a low p-value means "no, they're statistically
distinguishable"), and save to the correct path under `plots/`.

`main()` ties it together: discover available `(SPS,DPS)` file pairs for
matching `N` values, ask (or auto-pick, per `AUTO_MODE`) which pair(s), which
observables, and which plot types to produce, then loop over the selection
calling the drawing functions above.

---

## 5. `Convergence.py` — how observables change with sample size N

### Purpose
Where `Compare.py` looks at *one* fixed `N`, `Convergence.py` overlays the
**same observable at several different N values** (e.g. 10,000 through
800,000) on one set of axes, to answer: "has this distribution's shape
stopped changing as I add more Monte Carlo statistics, or is what I'm seeing
still just noise from too few events?"

### Shared machinery with `Compare.py`
`ALL_OBSERVABLES`, `_adaptive_edges`, `_log_edges`, `_choose_edges`,
`_hist`, `_kde`, `safe_input` all follow the **same reasoning** as the
identically-named (or near-identical) pieces in `Compare.py` explained
above — deliberately kept consistent so a reader who understands one script
already understands the corresponding piece of the other. The differences
worth calling out specifically:

```python
def make_palette(n_items):
def _alpha(idx, total):
```
Since this script can overlay many curves (one per N value) on a single
axis, it needs a colour/transparency scheme that stays readable with, say,
8 different N values plotted at once — `make_palette` assigns each N a
distinct colour, `_alpha` scales transparency so overlapping curves (which
will happen a lot, since higher-N curves should converge toward the same
shape as lower-N ones) don't turn into an illegible solid block.

```python
def _dir_single(base, proc):
def _dir_multi(base, proc):
def _nsuf(n_list):
def _ind_path(base, proc, tag, n_list):
def _sum_path(base, proc, tags, n_list):
```
Output-path builders: since a convergence plot's *filename* now depends on
*which set of N values* were plotted (not just one N), `_nsuf` builds a
filename suffix like `_N10000_50000_800000` from whatever subset was
actually chosen, so different N-selections never silently overwrite each
other's output files.

```python
def _cache_key(out_path):
def _load_cache():
def _save_cache(cache):
def _needs_replot(out_path, source_paths):
def _record(out_path, source_paths):
def _ask_replot(out_path, reason):
```
A small on-disk cache (keyed by output path, storing the modification time
and size of every source ROOT file that went into producing it) that lets
the script detect: "has anything actually changed since I last produced this
exact plot?" If nothing changed, it asks before silently re-doing (and
overwriting) work that would produce byte-identical output — useful when
convergence plots can be slow to regenerate across many N values and many
observables.

```python
def draw_panel(...):
def make_individual(...):
def make_summary(...):
def print_stats_table(...):
```
The actual per-observable and multi-observable-grid drawing functions,
following the same histogram/error-bar drawing logic as `Compare.py`'s
`plot_panel`, but looping over multiple `(N, dataset)` pairs instead of a
fixed SPS/DPS pair.

```python
def ask_mode_menu(sps_avail, dps_avail):
def ask_n_selection(pairs, label):
```
Interactive menus specific to this script's job: which *process* (SPS, DPS,
or both) and which *subset of available N values* to include in the
convergence overlay.

---

## 6. `Analyze.py` — the full quantitative SPS-vs-DPS study

### Purpose
The most analysis-heavy script: rather than just overlaying distributions,
it runs a battery of **statistical tests and derived-observable studies**
meant to answer "how distinguishable are SPS and DPS, quantitatively, and
is that distinguishability robust or a statistics artifact."

### `OBSERVABLES` and `CORR_PAIRS`
`OBSERVABLES` is the same kind of catalogue as the other two scripts (kept
independently per-file rather than imported from a shared module, per the
project's stated goal of not restructuring the existing architecture).
`CORR_PAIRS` is a separate list specifically of *pairs* of raw kinematic
branches (e.g. `(W_eta, Z_eta)`) used for the 2D correlation/joint-density
analyses (§ `analysis_2d`, `analysis_corr_matrix` below) — this needs to be
a distinct list from `OBSERVABLES` because those functions inherently work
on *two* variables at once, not one.

### `_paired_ns` / `discover` / `load`
Same file-discovery and branch-loading logic as `Compare.py`, adapted to
this script's own needs (e.g. `_paired_ns` specifically finds every `N` for
which *both* an SPS and a DPS file exist, since every analysis here needs a
matched pair).

### The `analysis_*` functions — one per statistical test
Each of these is a self-contained "run this one test/plot, save its output,
print its numbers" unit:

- **`analysis_ratio`**: for each observable, plots the SPS and DPS
  normalized histograms plus their bin-by-bin **ratio** (SPS/DPS) — a
  ratio flat at 1 means the two shapes agree in that region; deviations
  show exactly where and how strongly they differ.
- **`analysis_corr_matrix`**: computes and displays the full pairwise
  Pearson correlation matrix among the kinematic variables in
  `CORR_PAIRS`, separately for SPS and DPS — directly testing the core
  physics claim of this whole project (SPS bosons should be correlated by
  sharing one hard vertex; DPS bosons should be much closer to
  uncorrelated).
- **`analysis_2d`**: 2D joint-density heatmaps (like the `η_W`⊗`η_Z` one in
  `Compare.py`, but generalized to every pair in `CORR_PAIRS`) — the direct
  visual counterpart to the correlation matrix's single numbers.
- **`analysis_deltaR`**: `ΔR(W,Z) = √(Δη² + Δφ²)`, the standard collider
  angular-separation observable combining both angles into one number.
- **`analysis_shape`**: skewness, kurtosis, and other distribution-shape
  statistics for each observable — a more quantitative alternative to
  eyeballing "does this histogram look different."
- **`analysis_ks_vs_n`**: runs the KS test (SPS vs DPS) separately at
  *every available N*, then plots KS statistic vs N — directly answering
  "does the observed SPS/DPS separation hold up as statistics increase, or
  does it wash out (meaning it was a low-statistics fluctuation all
  along)."
- **`analysis_fb`**: forward-backward asymmetry, `(N_forward − N_backward) /
  (N_forward + N_backward)`, computed on **signed** rapidity-type
  observables (this is exactly where the earlier `dY_WZ`→`dY_signed` bug
  fix mattered — an asymmetry computed on an always-non-negative branch is
  meaningless by construction, since "backward" would never register any
  events).
- **`analysis_chi2`**: bin-by-bin χ² goodness-of-fit between the SPS and
  DPS histograms, on **fixed**-width bins deliberately (not the adaptive
  binning used for visual plots elsewhere) — because a χ² test's meaning
  depends on its bin definition being fixed *independently* of the data
  being tested; adaptive bins built from the data itself would make the
  test circular.
- **`analysis_uniformity`**: the χ²-vs-flat test built specifically to give
  a rigorous, numeric answer to "is DPS's `Δφ` distribution actually flat,
  or does it have real (if subtle) structure" — again using fixed bins and
  the *folded* `|Δφ|` branch (not the signed version) for the same
  "test needs a range where 'flat' is a meaningful null hypothesis" reason.

### `ANALYSES`, `ask_analyses`, `ask_n_selection`, `main`
`ANALYSES` is the master menu list (letter code, description, dispatch key)
tying every `analysis_*` function to a human-readable menu entry;
`ask_analyses`/`ask_n_selection` are the interactive selection prompts;
`main()` loads the requested N's paired data and dispatches to whichever
`analysis_*` functions were selected.

---

## 7. `Makefile` — the build system

```make
MAKEFLAGS += --no-builtin-rules
.SUFFIXES:
```
Disables Make's built-in implicit rules (e.g. its default `%.o: %.cc` C++
compilation rule). This matters because without it, a command like
`make GenerateSPS.cc` would get silently intercepted by Make's own built-in
pattern-matching *before* this Makefile's own explicit rule for that exact
target name gets a chance to run — a subtle, easy-to-hit trap when target
names happen to look like source filenames.

```make
$(BUILD)/.stamp:
	@mkdir -p $(BUILD)
	@touch $@
```
The binaries' compile rules depend on `$(BUILD)/.stamp` (an empty marker
file) as an order-only prerequisite for "the build directory must exist,"
rather than depending on the literal directory name `$(BUILD)` itself. This
avoids a real bug that was hit earlier in the project: if `BUILD=build`
(the directory name) and there's *also* a phony target literally named
`build`, Make can't tell the difference between "depend on the directory
existing" and "depend on the phony target build," causing a circular
dependency warning. The stamp-file indirection sidesteps the naming
collision entirely.

```make
all:
	@chmod +x interactive_all.sh
	@./interactive_all.sh

all_auto: build
	...
	@python3 Compare.py --auto
```
`make all` (the default, interactive entry point) doesn't build or run
anything itself — it delegates entirely to `interactive_all.sh`, which owns
the actual interactive loop (see §8). `make all_auto` is the
**non-interactive** equivalent used by `run_all.sh` and any automated/CI
context: build, run both generators, then call `Compare.py --auto` so
nothing ever blocks waiting on a prompt.

```make
NSPS := $(shell grep -E 'constexpr\s+int\s+N_EVENTS_SPS\s*=' Settings.h | grep -oE '[0-9]+' | tail -1)
NDPS := $(shell grep -E 'constexpr\s+int\s+N_EVENTS_DPS\s*=' Settings.h | grep -oE '[0-9]+' | tail -1)
```
Make reads the *current* event count directly out of `Settings.h` at
Makefile-parse time (via a `grep`/regex extraction, not by trying to
actually parse C++), so every binary/output filename automatically stays in
sync with whatever `settings_editor.py` last wrote there — no separate
"remember to update N here too" step.

The remaining targets (`build`, individual `GenerateSPS`/`GenerateDPS`
rules, `compare`/`convergence`/`analyze` Python-script shortcuts, `info`,
`clean*`) follow standard Make conventions: each binary's rule depends on
its source `.cc` file and `Settings.h`, so editing either automatically
triggers a rebuild next time; `clean_root`/`clean_plots`/`clean_all` give
explicit, separate ways to wipe generated ROOT files vs. plots vs.
everything, rather than one all-or-nothing clean.

---

## 8. `interactive_all.sh` — the `make all` driver

```bash
set -uo pipefail   # NOTE: no -e
```
Deliberately **omits** `-e` (which would abort the whole script on any
command's non-zero exit code) — explained in the comment: this script must
survive a failed submenu action (e.g. a failed compile) and return to its
own menu loop rather than dying outright. `-u` (error on unset variables)
and `pipefail` (a pipeline's exit status is its last *failing* command, not
just its last command) are kept, since those failure modes genuinely
indicate bugs worth catching.

```bash
do_build() { ... }
do_run() { ... }
```
Two shell functions wrapping the compile step and the run step
respectively, each independently callable — `do_build` reports how many
binaries were actually (re)compiled this round vs. already up to date (via
the same up-to-date timestamp check the Makefile itself uses:
`[ "$src" -nt "$bin" ]`, i.e. "is the source newer than the binary").

```bash
while true; do
    echo "  What do you want to do?"
    ...
    read -r -p "  Choice [1-4]: " choice
    case "$choice" in
        1) do_build; do_run ;;
        2) python3 settings_editor.py; do_build; do_run ;;
        3) python3 Compare.py ;;
        4) break ;;
        *) echo "  Enter 1, 2, 3, or 4." ;;
    esac
done
```
The main interactive loop, deliberately **never** auto-building on launch
(the menu is shown first, before anything happens) — per an earlier explicit
design change in this project, so that simply typing `make all` and then
immediately exiting never has an unwanted side effect of silently
compiling/running anything you didn't ask for. Option 2 (change event
counts) chains directly into a rebuild afterward, since a changed `N` is
meaningless until the binaries are recompiled against it.

---

## 9. `run_all.sh` — non-interactive pipeline

Structurally a simpler, linear version of the same compile→run→analyze
pipeline, with no menu at all — every step just runs unconditionally in
sequence: compile both generators if needed (same up-to-date check as
`interactive_all.sh`), run both in parallel (`&` backgrounding + `wait`, so
SPS and DPS generation happen concurrently rather than one after another,
since they're fully independent of each other), then call
`python3 Compare.py --auto` at the end. `set -euo pipefail` (note: **with**
`-e` this time) is appropriate here, unlike the interactive script — this is
meant to be a script you can run unattended (e.g. from a cron job or a CI
pipeline) where a failure *should* immediately stop everything rather than
silently continuing.

```bash
if [ $# -ge 1 ]; then NSPS=$1; NDPS=$1; fi
if [ $# -ge 2 ]; then NSPS=$1; NDPS=$2; fi
```
Command-line overrides: no arguments uses whatever `Settings.h` currently
says; one argument sets *both* SPS and DPS event counts to that value;
two arguments sets them independently.

---

## 10. `settings_editor.py` — interactive `Settings.h` editor

### Why this exists
Editing `Settings.h` by hand is error-prone (typos in the `constexpr`
syntax could silently fail to update the intended value, or worse, break
compilation). This script provides a numbered menu, reads and rewrites
`Settings.h` via regex, and keeps a `.bak` backup before every save.

```python
PARAMS = [
    ("Centre-of-mass energy", "ECM_GEV", float, "Beam energy", "GeV"),
    ...
]
```
Same "one table drives everything" pattern as the plotting scripts'
observable catalogues — each entry names a display label, the exact C++
token to search for in `Settings.h`, its Python type (`float`/`int`/`bool`,
used to pick the right regex and the right input-parsing/validation logic),
a description, and a unit string for display.

```python
def get_value(content, token, typ):
    if typ == bool:
        m = re.search(rf"constexpr\s+bool\s+{token}\s*=\s*(true|false)", content)
        ...
```
Regex-based extraction rather than an actual C++ parser — deliberately
simple, matching the exact `constexpr <type> <TOKEN> = <value>;` pattern
this project's `Settings.h` always uses. `\s*` allows for whitespace
variation (tabs vs spaces, extra alignment spaces) without requiring exact
formatting.

```python
def set_value(content, token, typ, new_val):
    ...
    new_content, n = re.subn(pattern, repl, content)
    if n == 0:
        print(f"  WARNING: token '{token}' not found in Settings.h")
    return new_content
```
`re.subn` returns both the modified string *and* a count of how many
substitutions were made — checking `n == 0` catches the case where the
regex simply didn't match anything (e.g. someone manually reformatted
`Settings.h` in a way this script's pattern no longer recognizes), warning
rather than silently doing nothing.

```python
def save(content):
    backup = SETTINGS_FILE + ".bak"
    shutil.copy2(SETTINGS_FILE, backup)
    with open(SETTINGS_FILE, "w") as f:
        f.write(content)
```
Always copies the *current on-disk* file to `.bak` immediately before
overwriting it — a simple, single-level undo mechanism (the `R` menu
command restores from this backup) in case a save turns out to be a
mistake.

```python
def main():
    ...
    if "--show" in sys.argv or "-s" in sys.argv:
        print_current(read_settings())
        return
    ...
    while True:
        content = read_settings()
        print_current(content)
        print("  Enter a parameter number to change it, or: A / S / Q / R")
        cmd = input("\n  Command: ").strip().upper()
        ...
```
`--show` is a fast, non-interactive path (just print current values and
exit — used by other scripts/humans who just want to check current settings
without entering the full menu). The main loop re-reads `Settings.h` from
disk at the top of every iteration (not just once at startup) — meaning if
you use the `R` (restore backup) command mid-session, or if the file changed
on disk for any other reason, the displayed values are always accurate
rather than a stale in-memory copy.

The per-parameter edit flow (`cmd` is a number → `prompt_value(...)` →
`set_value(...)` → immediate write to disk) writes the change to disk
**immediately** after each individual edit (not batched until you choose
"S" to save) — so if you change three parameters in one session, each one
is durably saved as soon as you confirm it, and only the very first change
in a session goes through the explicit backup-then-overwrite `save()` path
via the initial `S`/`A` command (subsequent individual-parameter edits
bypass `save()` and write directly, which is a deliberate simplification:
the single backup made at the very first `S` covers "what it looked like
before this session started editing," which is the meaningful restore
point).

---

## Summary: how all the pieces connect

```
Settings.h  ──────────────┬──────────────────────┐
   (shared constants,       │                      │
    shared helper fns)      ▼                      ▼
                      GenerateSPS.cc         GenerateDPS.cc
                       (SecondHard OFF)       (SecondHard ON +
                                                colour-reconnection/
                                                dipole-recoil fix)
                             │                      │
                             ▼                      ▼
                     root/SPS_<N>.root      root/DPS_<N>.root
                             │                      │
              ┌──────────────┼──────────────────────┤
              ▼              ▼                      ▼
        Compare.py    Convergence.py           Analyze.py
        (one N,       (many N, one            (full statistical
         SPS vs DPS    process at a            battery: ratio,
         overlay)      time, shape vs N)        correlation, KS
                                                 vs N, chi2, etc.)
              │              │                      │
              └──────────────┴──────────────────────┘
                             ▼
                        plots/*.pdf

Makefile / run_all.sh / interactive_all.sh
   orchestrate compiling the generators, running them, and invoking
   the three Python scripts above, using settings_editor.py to change
   Settings.h (and therefore N, MPI, energy, mass window) in between runs.
```

Every observable added anywhere in this project follows the same three-step
pattern: (1) compute it in the C++ generators and store it as a ROOT branch
(or derive it from existing branches in Python, if it's just a
transformation of data already stored), (2) add one entry to the
`ALL_OBSERVABLES`/`OBSERVABLES` catalogue in each Python script that should
plot it, (3) nothing else — every plotting function, menu, filename builder,
and statistical test reads from that one catalogue, so no observable-specific
code needs to be written anywhere else.
