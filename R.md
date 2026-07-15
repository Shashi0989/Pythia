# ZW_DPS — Full Project Documentation

This document explains every file in the project: what each block of code
does, the physics reasoning behind it, why it's written the way it is, and
how the files connect to each other. It reflects the **current** state of
the project — a two-generator pipeline (SPS and DPS only; there are no
standalone `GenerateZ`/`GenerateW` generators and no event-mixing DPS
model). Repetitive boilerplate (e.g. registering ROOT branches) is
explained once as a pattern rather than re-explained line by line, but
every distinct piece of logic is covered.

---

## 0. The physics question this project answers

You are comparing two different mechanisms by which a single
proton-proton collision at the LHC can produce a **W boson and a Z boson
simultaneously**:

- **SPS (Single Parton Scattering)** — one quark from proton A and one
  antiquark from proton B annihilate in a single hard interaction,
  `qq̄' → WZ`. Both bosons come from the same vertex — the same Feynman
  diagram, the same momentum transfer. They are kinematically tied
  together from the moment they're created: momentum conservation at that
  one vertex forces them to be (nearly) back-to-back in azimuth, correlates
  their rapidities (they share one longitudinal boost), and constrains
  their transverse momenta to nearly cancel.

- **DPS (Double Parton Scattering)** — two independent hard scatters
  happen within the *same* proton-proton collision. One parton pair
  produces the Z; a separate, unrelated parton pair produces the W. There
  is no Feynman diagram connecting them. The only thing they share is the
  underlying event: the same two protons, the same beam remnants, the same
  soft QCD activity (multi-parton interactions, colour reconnection).

The project generates both samples with PYTHIA 8, computes a battery of
kinematic and angular observables for each, and quantifies how
distinguishable the two production mechanisms are — plus studies how much
Monte Carlo statistics (how many events, N) is needed before those
differences are trustworthy rather than noise.

### Data flow

```
Settings.h  ───────────────────┬────────────────────┐
 (shared constants,             │                     │
  shared helper functions)      ▼                     ▼
                         GenerateSPS.cc          GenerateDPS.cc
                        (SecondHard OFF)     (SecondHard ON + colour-
                                              reconnection / dipole-
                                              recoil correlation fix)
                                │                     │
                    (raw per-particle kinematics only)
                                ▼                     ▼
                        root/SPS_<N>.root      root/DPS_<N>.root
                                │                     │
                 ┌──────────────┼─────────────────────┤
                 ▼              ▼                     ▼
           Compare.py    Convergence.py           Analyze.py
          (one N, SPS    (many N, one process     (full statistical
           vs DPS         at a time, shape vs N)   battery: ratio,
           overlay)                                correlation, KS
                                                     vs N, chi², etc.)
                 │              │                     │
                 └──────────────┴─────────────────────┘
                                ▼
                           plots/*.pdf
```

`Makefile` / `run_all.sh` / `interactive_all.sh` orchestrate compiling the
generators, running them, and invoking the three Python scripts above.
`settings_editor.py` changes `Settings.h` (and therefore N, MPI, beam
energy, mass window) in between runs.

Both generators write a tree named **`Events`** with an **identical
branch schema** — every branch name means the same physical quantity in
both files. Critically, **both generators now write only the raw,
per-particle kinematics PYTHIA/ROOT hand back directly** (pt, eta, phi,
mass/rapidity, and the raw px, py, pz, E four-vector, for every relevant
particle), plus the handful of event-record-only quantities that
genuinely cannot be recomputed outside PYTHIA (`nJets30`, `W_charge`,
`weight`). Every *derived*, cross-particle observable — Δη/Δφ/Δy between
the two bosons, the WZ invariant mass and system pT, the scalar pT sum,
the relative pT imbalance, and every lepton-level Δφ — is **no longer
computed or stored by the generators at all**. It is computed entirely in
the Python analysis layer instead, once per file load, inside each
script's own `load()` function (`Compare.py` / `Convergence.py` /
`Analyze.py`), via vectorized NumPy operating on the raw branches. This
is a deliberate design split: the C++ generators are pure data producers
("generate events, save exactly what PYTHIA measured"), and every
downstream calculation lives in the analysis scripts that already own
that responsibility. This is what lets every Python script load an
`SPS_*.root` and a `DPS_*.root` and plot them directly against each other
with no file-type-specific logic anywhere — and it also means the three
Python scripts, not the two generators, are now the single source of
truth for how every derived observable is defined.

---

## 1. `Settings.h` — shared constants and helper functions

A C++ header (no `main()`), `#include`d by both generators. Anything that
must be identical across both samples (beam energy, decay forcing, the
Δφ folding convention) lives here exactly once, so the two samples are
never accidentally generated under different physics assumptions.

```cpp
#pragma once
```
Include-guard: ensures the header's contents are inserted only once per
translation unit. `#pragma once` is a non-standard but universally
supported shorthand for the classic `#ifndef`/`#define`/`#endif` guard.

```cpp
constexpr double ECM_GEV = 13600.;
```
`constexpr` is a compile-time constant — the compiler inlines its value
everywhere it's used and it can't be reassigned at runtime. `13600.` GeV =
13.6 TeV, the LHC Run 3 centre-of-mass collision energy. This single
number defines "how energetic is the collision" for every generated event.

```cpp
constexpr double Z_MASS_MIN = 75.;
constexpr double Z_MASS_MAX = 120.;
```
The invariant-mass window PYTHIA restricts the Z-producing hard process
to. At low invariant mass, the intermediate boson is dominated by an
off-shell photon (γ*) rather than a real Z — the amplitude blows up as
mass → 0 (photon propagator ∝ 1/q²). Restricting to 75–120 GeV keeps the
sample dominated by genuine on-shell Z bosons (M_Z ≈ 91.2 GeV sits
comfortably inside this window) and removes the low-mass Drell-Yan
continuum, which isn't the physics this project studies.

```cpp
constexpr int N_EVENTS_SPS = 10000;
constexpr int N_EVENTS_DPS = 10000;
```
Default event counts, overridable from the command line
(`./GenerateSPS_10000 50000`) or via `settings_editor.py`. Plain `int`
(not `long`) because these values comfortably fit in 32 bits and every
other piece of code that reads them (bash regex, Python parsing) expects a
short decimal integer.

```cpp
constexpr bool ENABLE_MPI = true;
```
Multi-Parton Interactions: whether, on top of the one hard process
explicitly requested, PYTHIA also simulates additional softer
parton-parton scatters in the same collision (the "underlying event"). ON
because a real pp collision always has MPI activity; OFF is only useful
for fast debugging.

```cpp
constexpr int ID_Z      =  23;
constexpr int ID_WPLUS  =  24;
constexpr int ID_WMINUS = -24;
```
PDG particle IDs — the standard integer codes every particle-physics tool
uses to identify particle species (23 = Z⁰, 24 = W⁺, −24 = W⁻; a
particle/antiparticle pair gets opposite-sign IDs in the PDG scheme).
Used throughout the generators to search PYTHIA's event record by
matching `particle.id()` against these constants, instead of scattering
magic numbers through the code.

```cpp
#include "Pythia8/Pythia.h"
#include <cmath>
#include <string>
```
Brought in after the constants because the helper functions below need
`Pythia8::Pythia`, `std::string`, and `M_PI`/`fabs` from `<cmath>`.

```cpp
inline void setQuiet(Pythia8::Pythia& p) {
    p.readString("Print:quiet      = on");
    p.readString("Next:numberCount = 0");
}
```
`inline` on a free function defined in a header avoids "multiple
definition" linker errors if the header ends up compiled into more than
one translation unit that link together — each gets its own identical
copy, which is fine.

`pythia.readString(...)` is PYTHIA 8's universal configuration interface:
every tunable setting is set by feeding a string of the form
`"Setting:Name = value"`, parsed at runtime against PYTHIA's internal
settings database. (This is why a typo causes a *runtime* "setting not
found" error rather than a compile error — the string is opaque to the
C++ compiler.) `Print:quiet = on` suppresses PYTHIA's normal console
spam; `Next:numberCount = 0` disables the periodic "now processing event
N" progress print — useful for a script but noisy for a batch job
generating hundreds of thousands of events.

```cpp
inline void setBeam(Pythia8::Pythia& p) {
    p.readString("Beams:idA = 2212");
    p.readString("Beams:idB = 2212");
    p.readString("Beams:eCM = " + std::to_string(ECM_GEV));
}
```
`2212` is the PDG ID for a proton; setting both beams to `2212` configures
a proton-proton collider. `std::to_string(ECM_GEV)` converts the `double`
constant into the string PYTHIA's settings parser expects. Centralising
this ensures every generator uses identical beam conditions — if the
beams differed between SPS and DPS, the comparison would be meaningless.

```cpp
inline void setMPI(Pythia8::Pythia& p) {
    p.readString(std::string("PartonLevel:MPI = ") + (ENABLE_MPI ? "on" : "off"));
}
```
A ternary picks `"on"`/`"off"` from the compile-time `ENABLE_MPI` flag.
This is the single choke point that turns MPI on/off consistently for
every generator that calls `setMPI()`. Note `PartonLevel:MPI` (many soft
extra scatters) is a distinct mechanism from `SecondHard` (one additional
*hard* scatter, used by the DPS generator — see §3).

```cpp
inline void setZDecayMuMu(Pythia8::Pythia& p) {
    p.readString("23:onMode    = off");
    p.readString("23:onIfMatch = 13 -13");
}
```
PYTHIA's particle-decay-table syntax: `"23:onMode = off"` turns off all
default Z decay channels; `"23:onIfMatch = 13 -13"` re-enables
specifically the channel whose daughters are `{13, -13}` — PDG 13 = μ⁻,
−13 = μ⁺. Every Z in the sample is forced to decay to μ⁺μ⁻, and only that
channel. Two reasons: muons are the cleanest final state to reconstruct
(no missing energy like taus/neutrinos, no jet-clustering ambiguity like
hadronic decays); and forcing a single decay channel means every event has
the exact same final-state topology, so every lepton-level observable
computed downstream is directly comparable event-to-event.

```cpp
inline void setWDecayMuNu(Pythia8::Pythia& p) {
    p.readString("24:onMode      = off");
    p.readString("24:onIfMatch   = 13 -14");   // W+
    p.readString("-24:onMode     = off");
    p.readString("-24:onIfMatch  = -13 14");   // W-
}
```
The W needs two separate blocks — one for W⁺ (PDG 24), one for W⁻ (PDG
−24) — because, unlike the Z, they are genuinely different particles in
PYTHIA's decay table. Charge conservation dictates the allowed daughters:
W⁺ (charge +1) → μ⁺ (id −13) + ν_μ (id +14); W⁻ (charge −1) → μ⁻ (id +13)
+ ν̄_μ (id −14). If you ever need to edit this block, check it against
PYTHIA's own decay-table listing (`pythia.particleData.list(24)`) rather
than re-deriving the signs by hand — getting the pairing wrong wouldn't
cause a compile error, it would just silently leave that channel
disabled.

```cpp
inline double dPhi(double phi1, double phi2) {
    double d = std::fabs(phi1 - phi2);
    return (d > M_PI) ? 2.0 * M_PI - d : d;
}
```
Azimuthal angle φ lives on a circle, conventionally reported in (−π, π].
A naive `|φ1 − φ2|` can return up to 2π, double-counting the same
physical separation (e.g. φ1 = −179°, φ2 = 179° are only 2° apart on the
circle, but the naive difference gives 358°). This function folds the raw
difference into `[0, π]`: if it exceeds π, the true separation is `2π`
minus that difference (the short way around the circle). This is the
standard definition used throughout collider physics.

Note: as of the current generators, this `dPhi()` helper is declared in
`Settings.h` but is **not actually called from the C++ generators
anymore** — neither `GenerateSPS.cc` nor `GenerateDPS.cc` computes any
Δφ observable itself. Its folding convention lives on only as the
reference definition that the Python analysis layer's own `_dphi_np()`-
style helpers reimplement (see §4) — every `dPhi_*` observable in every
plot and every test in the project still ends up folded to `[0, π]`, just
computed in NumPy on the Python side instead of in C++ per-event.

---

## 2. `GenerateSPS.cc` — the SPS generator

### Includes and setup

```cpp
#include "Settings.h"
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
```
`TFile`/`TTree` are ROOT's file-format and tabular-data classes — ROOT is
the standard HEP data-analysis framework, and its `.root` files
(structured, compressed, column-oriented binary files) are what every
downstream Python script reads via `uproot`. `using namespace Pythia8;`
and `using namespace std;` are acceptable in a small, single-purpose
`.cc` file (would be bad practice in a shared header, which is why
`Settings.h` avoids it).

```cpp
int nEvents = N_EVENTS_SPS;
if (argc > 1) nEvents = atoi(argv[1]);
string outPath = "root/SPS_" + to_string(nEvents) + ".root";
```
Default event count comes from the compiled-in constant, but a
command-line argument overrides it. The output filename is stamped with
the event count (`SPS_10000.root`, not just `SPS.root`) — deliberate, so
re-running with a different N produces a new file instead of overwriting
the old one, which is exactly what makes the convergence study (comparing
the same observable across many N) possible: every historical run
survives on disk under its own name.

### PYTHIA configuration — the SPS physics process

```cpp
pythia.readString("WeakDoubleBoson:all      = off");
pythia.readString("WeakDoubleBoson:ffbar2ZW = on");
pythia.readString("SecondHard:generate      = off");
```
`WeakDoubleBoson` is PYTHIA's process class for diboson production from a
*single* hard vertex. `:all = off` starts from a clean slate, then
`:ffbar2ZW = on` enables specifically `f f̄' → Z W` — a fermion and
antifermion of different flavour annihilating directly into a Z and a W.
This is the Feynman diagram that gives SPS its defining physics: one
vertex, one hard scatter, momentum-conserving, both bosons correlated by
construction. `SecondHard:generate = off` explicitly disables PYTHIA's
machinery for a second, independent hard process in the same event (the
DPS-generating switch, see §3) — turned off here to guarantee the SPS
sample never accidentally contains a spurious second hard scatter.

```cpp
setZDecayMuMu(pythia);
setWDecayMuNu(pythia);
setQuiet(pythia);
if (!pythia.init()) { cerr << "ERROR: PYTHIA init failed.\n"; return 1; }
```
`pythia.init()` validates and locks in every setting fed to it above,
builds PYTHIA's internal process/phase-space machinery, and returns
`false` if anything was inconsistent (e.g. a misspelled setting name, an
invalid mass window, or — historically in this project — an invalid
setting string; see §3). Checking the return value and bailing out with a
nonzero exit code is what lets the Makefile/shell scripts detect a failed
run instead of silently producing a garbage or empty ROOT file.

### ROOT tree and branch declarations — raw kinematics only

```cpp
TFile* fOut = TFile::Open(outPath.c_str(), "RECREATE");
TTree* tree = new TTree("Events", "SPS WZ production (raw kinematics)");
```
`"RECREATE"` overwrites any existing file at `outPath`. The tree is named
`"Events"` — exactly the name every Python reader
(`uproot.open(path)["Events"]`) expects. The tree's title has been
updated to say "raw kinematics" — a direct reflection of the current
design: **this file stores nothing but what PYTHIA/ROOT hand back for
each individual particle.**

Each declared `double`/`int` variable (`W_pt`, `W_eta`, …) is a **row
buffer**: ROOT's `TTree::Branch` API binds a branch name to the memory
address of a C++ variable, and every `tree->Fill()` call copies whatever
value currently sits at that address into the next row of that branch's
column. This is why the fill loop always writes to the same named
variable and then calls `Fill()` — the address never changes, only the
value does.

The branches now fall into exactly these groups, and **no others**:

- **W boson raw kinematics**: `W_pt, W_eta, W_phi, W_mass, W_rapidity`
  plus the raw Cartesian `W_px, W_py, W_pz, W_E` (needed for any
  downstream calculation that must combine two 4-vectors, like an
  invariant mass — you can't add two pT scalars and get a physically
  meaningful "combined pT", you must add the vector components).
  `W_charge` records ±1 — an event-record-only quantity (a PDG-id-sign
  lookup) that genuinely can't be recomputed from the stored kinematics.
- **Z boson raw kinematics**: identical structure, `Z_*`.
- **W decay products**: `lep_*` (the W's charged lepton, μ⁻ or μ⁺
  depending on `W_charge`) and `nu_*` (its neutrino) — the actual
  final-state particles a real detector would see.
- **Z decay products**: `lep1_*` (always the μ⁻) and `lep2_*` (always the
  μ⁺) — the Z's own two decay leptons, labelled unambiguously by charge
  sign rather than "whichever muon PYTHIA listed first," so that e.g.
  a downstream Δφ(Z, lep1) means the same physical angle in every event.
- **Event-level**: `weight` (PYTHIA's per-event weight — effectively
  always 1.0 for this unweighted LO sample, but read from PYTHIA rather
  than hard-coded so the code stays correct if a weighted process is ever
  swapped in) and `nJets30` (a rough jet-multiplicity proxy, explained
  where it's computed below).

**What is *not* here anymore:** there are no `dEta_WZ`, `dPhi_WZ`,
`dY_WZ`, `M_WZ`, `pT_WZ`, `sT`, `pT_imbalance`, or any `dPhi_*`
lepton-pair branches. Every one of those cross-particle, derived
quantities used to be computed in this file's event loop and written as
its own branch; all of that logic has been deleted from the generator
entirely and now lives in each Python analysis script's `load()`
function instead (see §4), computed from the raw branches listed above.
The file's own header comment states this design choice explicitly:
*"generators only GENERATE, they do not ANALYSE."*

```cpp
tree->Branch("W_pt", &W_pt);
...
```
Each `tree->Branch("name", &variable)` call registers one column: "the
column called `name` is filled from whatever is currently in `variable`."
This block of calls is mechanical but essential — it *is* the ROOT file's
schema, and every downstream Python script's `d["W_pt"]`-style access
depends on these exact string names matching.

### The event loop

```cpp
long nFilled = 0, nNoW = 0, nNoZ = 0, nNoLepNu = 0;
for (int iev = 0; iev < nEvents; ++iev) {
    if (!pythia.next()) continue;
```
`pythia.next()` generates one full event (hard process → parton shower →
hadronization → decays) and returns `false` if that particular attempt
failed internally (rare, but possible for numerically difficult
phase-space points); `continue` skips to the next iteration without
incrementing any counter — the failed attempt is silently discarded (it
doesn't count against `nEvents` requested, and isn't logged as a distinct
failure category, since PYTHIA already retries internally before giving
up).

```cpp
int iZ = -1;
for (int i = 0; i < pythia.event.size(); ++i)
    if (pythia.event[i].id() == ID_Z) iZ = i;
if (iZ < 0) { ++nNoZ; continue; }
```
`pythia.event` is PYTHIA's full particle history for this event — not
just final-state particles, but every intermediate copy (a particle can
appear multiple times as it radiates, recoils, etc. through the shower).
This loop deliberately does **not** break on the first match — it keeps
overwriting `iZ` every time it finds another Z, so after the loop `iZ`
points at the *last* Z entry in the record. This matters: PYTHIA stores
intermediate/pre-shower copies of a resonance before its final,
post-recoil copy, and the last occurrence is the one whose kinematics you
actually want. If no Z is found at all, the event is skipped and counted
under `nNoZ`.

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
Same "last occurrence wins" logic for either W charge — since SPS's
`ffbar2ZW` process produces a W of either sign depending on which quark
flavours annihilated, `W_charge` is recorded per event rather than
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
`daughterList()` returns the indices of a particle's direct decay
products. Since `setWDecayMuNu` forced every W to decay to exactly
`{μ, ν_μ}`, this identifies which daughter is the muon (`|id|==13`,
`abs()` because it could be μ⁻ or μ⁺) and which is the neutrino
(`|id|==14`).

```cpp
int iLep1 = -1, iLep2 = -1;
for (int d : pythia.event[iZ].daughterList()) {
    int id = pythia.event[d].id();
    if (id == 13)  iLep1 = d;
    if (id == -13) iLep2 = d;
}
if (iLep1 < 0 || iLep2 < 0) { ++nNoLepNu; continue; }
```
Same idea for the Z's two muon daughters — but here the *signs* matter and
are used to distinguish them: `lep1` is always specifically the μ⁻
(id = +13; note the PDG convention: the μ⁻ *particle* has *positive* PDG
id 13, μ⁺ has −13 — a common point of confusion, since particle/
antiparticle sign convention for leptons doesn't track electric charge
sign directly), `lep2` is always the μ⁺. This gives a consistent,
unambiguous labelling so that any downstream `dPhi_Z_lep1`-type
observable means the same physical angle in every single event.

```cpp
nJets30 = 0;
for (int i = 0; i < pythia.event.size(); ++i) {
    const Particle& p = pythia.event[i];
    if (!p.isFinal()) continue;
    if (abs(p.id()) > 8 && p.id() != 21) continue;
    if (p.pT() > 30.) ++nJets30;
}
```
A deliberately crude jet-multiplicity proxy. `isFinal()` restricts to
actual final-state particles (not intermediate shower copies).
`abs(id) > 8 && id != 21` keeps only quarks (PDG 1–8) and gluons (21) —
i.e. only partons, not leptons/photons — with `pT() > 30.` GeV as an
arbitrary-but-standard-ish "hard enough to plausibly become a
reconstructed jet" threshold. This is explicitly *not* a real jet
algorithm (no clustering) — a fast approximation, good enough to
distinguish "basically no extra radiation" from "clearly some hard extra
activity" event by event. Like `W_charge`, this stays in the generator
because it requires walking the full PYTHIA particle list — it is not a
"derived observable" recomputable from any single stored branch.

```cpp
const Particle& pW = pythia.event[iW];
W_pt = pW.pT(); W_eta = pW.eta(); W_phi = pW.phi();
W_mass = pW.m(); W_rapidity = pW.y();
W_px = pW.px(); W_py = pW.py(); W_pz = pW.pz(); W_E = pW.e();
```
Once the right particle index is known, PYTHIA's `Particle` class exposes
all these kinematic quantities as simple accessor methods —
`pT() = √(px²+py²)`, `eta() = -ln(tan(θ/2))` (pseudorapidity),
`y() = ½ln((E+pz)/(E−pz))` (true rapidity), etc. — all computed internally
by PYTHIA from the particle's stored 4-momentum. The same pattern fills
the `Z_*` block from `pZ`, and then the `lep_*`/`nu_*`/`lep1_*`/`lep2_*`
blocks from the corresponding daughter particles.

```cpp
weight = pythia.info.weight();
tree->Fill();
++nFilled;
```
`pythia.info.weight()` retrieves this event's Monte Carlo weight.
`tree->Fill()` is the moment every currently-set branch variable is
copied into a new row of the tree. Immediately before this, the source
carries an explicit comment recording the design decision: every
cross-particle derived quantity that used to be computed here (`dEta_WZ`,
`dPhi_WZ`, `dY_WZ`, `M_WZ`, `pT_WZ`, `sT`, every lepton-level Δφ,
`pT_imbalance`, …) is now computed in the Python analysis layer's
`load()` functions instead, directly from the raw branches above — this
file only ever writes what PYTHIA measured.

### Wrap-up

```cpp
fOut->cd(); tree->Write(); fOut->Close(); pythia.stat();
```
`fOut->cd()` ensures the current ROOT directory context is this output
file; `tree->Write()` flushes the tree to disk; `fOut->Close()`
finalizes the file; `pythia.stat()` prints PYTHIA's own end-of-run
statistics (cross-sections, efficiency, error counts).

A final human-readable summary (`Events requested / saved / Skipped (no
W) / (no Z) / (no lep/nu)`) lets you sanity-check at a glance how many
requested events actually survived every selection cut — a large gap, or
a large `nNoLepNu` count, is a red flag worth investigating.

---

## 3. `GenerateDPS.cc` — the DPS generator

Structurally this file is extremely similar to `GenerateSPS.cc` — same
includes, same raw-kinematics-only branch-declaration pattern (identical
schema, including the same `lep1_*`/`lep2_*` branches described in §2),
same event loop skeleton, same lepton-finding logic, same jet-counting,
same summary printout, and the **same absence** of any derived,
cross-particle branch (no `dEta_WZ`/`dPhi_WZ`/`dY_WZ`/`M_WZ`/`pT_WZ`/
`sT`/`pT_imbalance`/`dPhi_*` — all of it lives in Python now, see §4).
This section focuses on everything that's genuinely *different* — which
is where the actual "DPS vs SPS" physics content lives.

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
production: `ffbar2gmZ` (`qq̄ → Z/γ*`), plus the two associated
Z+parton channels `qqbar2gmZg`/`qg2gmZq`, which let the Z recoil against
a hard extra parton at the matrix-element level rather than only via the
parton shower. By itself this is a completely ordinary, single-vertex
process — indistinguishable from "just generate some Z bosons." The
`Z_MASS_MIN`/`Z_MASS_MAX` window plays the same role as in the SPS
generator: keep this an on-shell Z, not low-mass γ*.

```cpp
pythia.readString("SecondHard:generate = on");
pythia.readString("SecondHard:SingleW  = on");
pythia.readString("SecondHard:WAndJet  = on");
```
This is the crux of the DPS implementation. `SecondHard:generate = on`
tells PYTHIA: in addition to the first hard process configured above,
generate a **second, statistically independent hard parton-parton
scattering within the same event** — PYTHIA's built-in machinery for
exactly the physics phenomenon this project studies. `SecondHard:SingleW`
and `SecondHard:WAndJet` restrict the second process to direct
`qq̄' → W±` (Drell-Yan-like) and W+parton respectively. So the full
picture per event: one parton pair inside this pp collision makes a Z, a
completely separate parton pair inside the *same* collision makes a W —
genuine DPS, not two bosons artificially glued together from different
events. This particular first-hard=Z / second-hard=W split was chosen
because PYTHIA's `SecondHard` machinery requires the first and second
hard processes to be drawn from compatible process classes, and this
combination (`WeakSingleBoson` first + `SecondHard:SingleW`/`WAndJet`
second) is the documented, supported way to get one Z and one W from two
independent hard scatters in a single PYTHIA 8 event.

### The correlation fix — restoring the expected Δφ(W,Z) enhancement

This block is the direct response to an earlier finding in this project:
a *naive* `SecondHard` DPS configuration produces a completely **flat**
Δφ(W,Z) distribution, but the physically expected behaviour (per the
project's supervisor) is a mild upward enhancement as Δφ → π — nowhere
near as strong as SPS's sharp back-to-back peak (which comes from
momentum conservation at one shared vertex, absent here), but a real,
non-zero residual correlation.

**Why a naive configuration is too flat:** as configured with just the
settings above, PYTHIA treats the two hard scatters as essentially fully
disconnected sub-events — independent colour flow, independent local
recoil. Real DPS in an actual pp collision is not quite that independent,
because both hard scatters still happen inside the same two protons:

1. **Shared beam remnants constrain both scatters** — whatever momentum,
   colour, and flavour content is pulled out of proton A to make the Z
   changes what's left over in proton A's remnant for the W-producing
   scatter to draw from (and symmetrically for proton B).
2. **Colour reconnection can act across the whole event** — colour
   strings from the Z-producing and the W-producing systems, both
   ultimately anchored in the same beam remnants, are not required to be
   colour-blind to one another during hadronization.
3. **Soft-radiation (ISR) recoil bookkeeping matters** — how the shower
   assigns recoil for radiated soft gluons affects each hard system's net
   transverse balance; doing this locally (dipole-based) rather than
   globally changes how much of that soft activity's momentum leaks
   between the two systems.

None of these three effects create anything like SPS's hard
momentum-conservation constraint — they're soft, indirect, "the two
systems share a collision, not a vertex" effects — but they are real and
are expected to produce a mild large-Δφ enhancement, which is exactly
what's missing from the naive flat baseline.

```cpp
pythia.readString("ColourReconnection:reconnect        = on");
```
Enables cross-system colour reconnection during hadronization: colour
lines from the Z system and the W system, both anchored in the same two
beam remnants, are allowed to reconnect with each other rather than being
forced to stay entirely separate.

```cpp
pythia.readString("SpaceShower:dipoleRecoil = on");
```
Switches initial-state-radiation recoil from PYTHIA's default global
scheme to a dipole scheme, where recoil from a radiated parton is
assigned locally to the specific colour dipole involved, rather than
smeared uniformly across the whole event. With two independent hard
systems present, this keeps each system's own hard kinematics physically
sensible while still allowing soft emissions near the shared beam
remnants to correlate the two systems' net transverse balance — exactly
the soft-correlation mechanism described above.

**What is deliberately *not* used, and why:** `MultipartonInteractions:
allowRescatter` (letting the same outgoing parton from one scatter
re-scatter again) is intentionally **not enabled**. The baseline DPS
sample is meant to isolate the soft beam-remnant/colour-reconnection
correlation mechanism specifically, not mix in a different (and
stronger) correlation mechanism at the same time — this is a controlled,
single-variable physics change, and the setting has been removed from the
file entirely (not merely commented out) rather than left as dead code.

Also removed entirely: `BeamRemnants:reconnectRange = 10.0`. This setting
is **not valid** in PYTHIA 8.316 — it would cause `pythia.init()` to
fail — and has been deleted rather than kept around commented-out, so a
future editor doesn't waste time reintroducing the same mistake.
`ColourReconnection:reconnect = on` above already enables cross-system
colour reconnection using PYTHIA's default reconnection range, which is
sufficient for the soft correlation mechanism described; no directly
equivalent replacement setting is needed.

### The rest of the file — raw kinematics only, same as §2

Branch declarations, the event loop, lepton-finding (for both the W's
`lep`/`nu` and the Z's `lep1`/`lep2`), and jet-counting are structurally
identical to `GenerateSPS.cc` — the same reasoning applies line-for-line
(see §2). As in the SPS generator, the fill block ends with an explicit
comment recording that every cross-particle derived quantity that used to
live here (`dEta_WZ`, `dPhi_WZ`, `dY_WZ`, `M_WZ`, `pT_WZ`, `sT`, every
lepton-level Δφ, `pT_imbalance`, …) has been moved out to the Python
analysis layer's `load()` functions, computed from the raw branches. The
two generators produce ROOT files with an **identical, raw-kinematics-only
branch schema**, which is exactly what makes it possible for `Compare.py`
(and every other analysis script) to load an `SPS_*.root` and a
`DPS_*.root` and plot them on the same axes.

---

## 4. `Compare.py` — SPS vs DPS comparison at a single N

### Purpose

Reads exactly one `SPS_<N>.root` and one `DPS_<N>.root` (same N), and
produces overlay plots, KS-test statistics, and summary tables for a
single fixed sample size.

### The observable catalogue

```python
ALL_OBSERVABLES = [
    ("dEta_signed",    r"$\Delta\eta(W,Z)$",              "",    (-5, 5),     False, "dEta",     True),
    ("dPhi_WZ",        r"$|\Delta\phi(W,Z)|$",             "rad", (0, 3.1416), False, "dPhi",     True),
    ("dY_signed",      r"$\Delta y(W,Z)$",                 "",    (-5, 5),     False, "dY",       True),
    ("M_WZ",           r"$M_{WZ}$",                        "GeV", (0, 1500),   True,  "MWZ",      True),
    ("Z_pt",           r"$p_T(Z)$",                        "GeV", (0, 500),    True,  "ZpT",      True),
    ("W_pt",           r"$p_T(W)$",                        "GeV", (0, 500),    True,  "WpT",      True),
    ("sT",             r"$S_T$",                           "GeV", (0, 800),    True,  "sT",       True),
    ("pT_WZ",          r"$p_T^{WZ}$",                      "GeV", (0, 400),    True,  "pTWZ",     True),
    ("pT_imbalance",   r"$|\vec{p}_T^{W}+\vec{p}_T^{Z}|/(p_T^W+p_T^Z)$", "", (0, 1), False, "pTimb", True),
    ("dPhi_Z_lep",     r"$|\Delta\phi(Z,\ell_W)|$",        "rad", (0, 3.1416), False, "dPhiZLep", True),
    ("dPhi_Z_nu",      r"$|\Delta\phi(Z,\nu_W)|$",         "rad", (0, 3.1416), False, "dPhiZNu",  True),
    ("dPhi_lep1_lep2", r"$|\Delta\phi(\ell_1,\ell_2)|$",   "rad", (0, 3.1416), False, "dPhiLL",   True),
    ("dPhi_lep_nu",    r"$|\Delta\phi(\ell,\nu)|$",        "rad", (0, 3.1416), False, "dPhiLNu",  True),
    ("dPhi_W_lep1",    r"$|\Delta\phi(W,\ell_1)|$",        "rad", (0, 3.1416), False, "dPhiWL1",  True),
    ("dPhi_W_lep2",    r"$|\Delta\phi(W,\ell_2)|$",        "rad", (0, 3.1416), False, "dPhiWL2",  True),
    ("dPhi_Z_lep1",    r"$|\Delta\phi(Z,\ell_1)|$",        "rad", (0, 3.1416), False, "dPhiZL1",  True),
    ("dPhi_Z_lep2",    r"$|\Delta\phi(Z,\ell_2)|$",        "rad", (0, 3.1416), False, "dPhiZL2",  True),
]
```
This is a single table-as-data-structure: every observable the script
knows how to plot is one tuple of `(branch_or_derived_key, LaTeX_label,
axis_unit, (x_min, x_max), use_log_y, short_filename_tag, use_histogram)`.
Every function downstream (menu printing, plotting, file naming) reads
this one list rather than any observable-specific `if`/`elif` chain — to
add or remove an observable from every plot type at once, you edit
exactly one line here. (`Analyze.py` and `Convergence.py` each keep their
own independent copy of this same catalogue, per the project's design
choice of not sharing a common module — see §9.) Note that **every
non-raw key in this table** — `dEta_signed`, `dPhi_WZ`, `dY_signed`,
`M_WZ`, `sT`, `pT_WZ`, `pT_imbalance`, and every `dPhi_*` lepton-pair
entry — is now a *derived* key computed inside `load()` (below), not a
ROOT branch; only `Z_pt` and `W_pt` are read straight off the tree.

`dEta_signed`/`dY_signed` are **not** ROOT branches; they're derived in
Python (see `load()` below) from the raw `W_eta`/`Z_eta` branches,
specifically so the plotted quantity can be signed (range −5 to 5) rather
than an always-non-negative `|Δη|`/`|Δy|`, which would lose the "which
boson is more forward" information. The range is clipped to ±5, matching
the finite pseudorapidity acceptance of a real detector rather than
plotting generator-level tails no real experiment could see.

`pT_imbalance`'s range `(0, 1)` reflects the vector definition described
in §2/§3: a ratio of a vector magnitude to a scalar sum, mathematically
guaranteed to lie in `[0, 1]`.

The eight `dPhi_*` entries after `pT_imbalance` are every physically
distinct pairwise azimuthal separation the project requires: `dPhi_Z_lep`
and `dPhi_Z_nu` (Z vs. the W's own lepton/neutrino), and
`dPhi_lep1_lep2`, `dPhi_lep_nu`, `dPhi_W_lep1`, `dPhi_W_lep2`,
`dPhi_Z_lep1`, `dPhi_Z_lep2` completing the full set among
`{W, Z, lepton1, lepton2, neutrino}`. Since the generators no longer store
any Δφ branch at all, every one of these is now computed in `load()` from
the raw `*_phi` branches, via the same `[0, π]`-folding convention as the
C++ `dPhi()` helper in `Settings.h` (reimplemented as a NumPy-vectorized
helper on the Python side — see below).

### AUTO_MODE / safe_input

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
unreliable, because several IDE "Run" consoles report `isatty()==False`
even though they do forward keystrokes — the script was silently
skipping its own menu. The fix uses an explicit flag instead: only
`--auto` on the command line (used by the Makefile/`run_all.sh` pipeline)
or the `ZW_DPS_AUTO=1` environment variable triggers non-interactive
mode. Plain `python3 Compare.py` always tries to prompt, and
`safe_input()` only falls back to auto-mode if `input()` genuinely raises
`EOFError` — so the script can never hang forever, but also never
silently skips a menu a real user could have answered.

### Loading data: `load()`

```python
def load(path, tree="Events"):
    """Load the RAW per-particle branches from a ROOT TTree via uproot,
    then compute every cross-particle derived observable here in Python.

    The generators (GenerateSPS.cc / GenerateDPS.cc) only ever write raw
    kinematics -- pt, eta, phi, mass/rapidity, and the raw px,py,pz,E
    four-vector -- for each relevant particle (W, Z, the W's lepton and
    neutrino, and the Z's two leptons), plus a couple of event-record-only
    quantities (nJets30, W_charge, weight) that can't be recomputed
    outside PYTHIA. Every *derived*, multi-particle observable -- angular
    gaps, the WZ invariant mass and system pT, the scalar pT sum, every
    lepton-level Delta-phi, and the relative pT imbalance -- is computed
    right here instead, once per file load, via vectorized NumPy across
    the whole sample. This keeps all analysis logic in the analysis
    scripts, not the generators."""
    print(f"  Loading {path} ...")
    d = uproot.open(path)[tree].arrays(library="np")
    print(f"    {len(d['W_pt']):,} events")

    # ---- WZ-system observables ----
    d["dEta_WZ"] = np.abs(d["W_eta"] - d["Z_eta"])
    d["dPhi_WZ"] = _dphi_np(d["W_phi"], d["Z_phi"])
    d["dY_WZ"]   = np.abs(d["W_rapidity"] - d["Z_rapidity"])

    sumPx_WZ = d["W_px"] + d["Z_px"]
    sumPy_WZ = d["W_py"] + d["Z_py"]
    sumPz_WZ = d["W_pz"] + d["Z_pz"]
    sumE_WZ  = d["W_E"]  + d["Z_E"]
    d["M_WZ"]  = np.sqrt(np.clip(sumE_WZ**2 - sumPx_WZ**2 - sumPy_WZ**2 - sumPz_WZ**2, 0, None))
    d["pT_WZ"] = np.hypot(sumPx_WZ, sumPy_WZ)
    d["sT"]    = d["W_pt"] + d["Z_pt"]

    # ---- pT imbalance: vector definition ----
    denom = d["W_pt"] + d["Z_pt"]
    d["pT_imbalance"] = np.where(denom > 0,
                                 np.hypot(sumPx_WZ, sumPy_WZ) / denom, 0.0)

    # ---- Lepton-level Delta-phi observables ----
    d["dPhi_Z_lep"]     = _dphi_np(d["Z_phi"], d["lep_phi"])
    d["dPhi_Z_nu"]      = _dphi_np(d["Z_phi"], d["nu_phi"])
    d["dPhi_lep1_lep2"] = _dphi_np(d["lep1_phi"], d["lep2_phi"])
    d["dPhi_lep_nu"]    = _dphi_np(d["lep_phi"], d["nu_phi"])
    d["dPhi_W_lep1"]    = _dphi_np(d["W_phi"], d["lep1_phi"])
    d["dPhi_W_lep2"]    = _dphi_np(d["W_phi"], d["lep2_phi"])
    d["dPhi_Z_lep1"]    = _dphi_np(d["Z_phi"], d["lep1_phi"])
    d["dPhi_Z_lep2"]    = _dphi_np(d["Z_phi"], d["lep2_phi"])

    # ---- Derived: SIGNED Δη and Δy (not folded to |...|) ----
    d["dEta_signed"] = d["W_eta"]      - d["Z_eta"]
    d["dY_signed"]   = d["W_rapidity"] - d["Z_rapidity"]

    return d
```
`uproot` is a pure-Python library for reading ROOT files without needing a
full ROOT installation with Python bindings. `.arrays(library="np")`
returns a dict mapping every branch name to a NumPy array of that
branch's values across all events — so computing a derived quantity is a
single vectorized NumPy operation across the entire dataset (not a
Python-level per-event loop, which matters for performance with
hundreds-of-thousands-of-rows arrays). `load()` now does **all** of the
derivation work that the C++ generators used to do: `_dphi_np()` is the
NumPy-vectorized counterpart of `Settings.h`'s C++ `dPhi()` helper (fold
the raw angular difference to `[0, π]`), and every WZ-system quantity
(`dEta_WZ`, `dPhi_WZ`, `dY_WZ`, `M_WZ`, `pT_WZ`, `sT`), the vector
`pT_imbalance`, all eight lepton-level `dPhi_*` observables, and the
signed `dEta_signed`/`dY_signed` pair are computed here from the raw
branches, exactly once per file load, rather than being read pre-computed
off the tree. `M_WZ`'s `np.clip(..., 0, None)` guards against a tiny
negative value under the square root from floating-point round-off right
at the (near-)massless-system edge case — the 4-vector-sum analogue of
`Vec4::mCalc()`'s own internal safety in PYTHIA.

### Binning strategy: `_adaptive_edges` / `_log_edges` / `_choose_edges`

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
Rather than equal-width bins, this builds bin edges from the data's own
quantiles — every bin contains (approximately) the same number of
events, giving many narrow bins where events are dense (resolving a sharp
peak well) and few wide bins where events are sparse (without those bins
being individually noisy). `np.unique` guards against duplicate edges
(possible with many identical values, e.g. from a hard phase-space cut);
the two size checks fall back to plain equal-width bins whenever there
isn't enough data for quantile binning to be meaningful. This is used for
the pT/mass-type observables' histograms (see below).

```python
def _log_edges(lo, hi, nbins, floor=1.0):
    floor = min(floor, hi * 0.01) if hi > 0 else floor
    if lo <= 0:
        tail_edges = np.geomspace(max(floor, 1e-6), hi, nbins)
        return np.concatenate(([0.0], tail_edges))
    return np.geomspace(lo, hi, nbins + 1)
```
For a steeply falling spectrum (pT(W), pT(Z), M_WZ) plotted on a log-y
axis, quantile binning is the *wrong* tool — since most events sit at low
pT, quantile bins get almost entirely used up describing the first small
slice of the range, leaving the long tail covered by a handful of
enormous, blocky bins. Log-spaced bins (`np.geomspace`, equal-ratio
rather than equal-width) instead grow smoothly all the way through the
tail, resolving both the peak and the tail properly. `floor` handles the
fact that log-spacing is undefined starting exactly at 0 — everything
below the smallest meaningful geometric step lands in one small `[0,
floor)` bin, and geometric spacing takes over from there.

### Plotting strategy: KDE for angular observables, histograms for pT/mass

The current script's own header now documents an explicit split in how
different observable families are drawn:

- **Angular observables (`dEta`, `dPhi`, `dY`) and any Δ*R*-type
  quantity**: drawn with a **Gaussian KDE** (`scipy.stats.gaussian_kde`),
  since these distributions are smooth and bounded, and a KDE gives a
  cleaner visual curve than a histogram for that shape.
- **pT observables (`Z_pt`, `W_pt`, `sT`, `pT_WZ`) and `M_WZ`**: drawn
  with **normalized histograms and Poisson error bars** instead, using
  the `_adaptive_edges`/`_log_edges` machinery above. A KDE is
  deliberately *not* used here because it fails on the DPS sample's
  heavy-tailed pT distributions — rare hard `SecondHard:WAndJet` events
  inflate the sample variance, which causes the KDE's Silverman
  bandwidth rule to pick an absurdly wide or narrow kernel. Fixed/adaptive-
  bin histograms are the numerically robust tool for these observables,
  so the script keeps them for that family rather than trying to force a
  single drawing method onto every observable.

### `_choose_edges`, `_norm_hist`, and the drawing functions

```python
def _choose_edges(pooled_vals, lo, hi, nbins, logy):
    return _log_edges(lo, hi, nbins) if logy else _adaptive_edges(pooled_vals, lo, hi, nbins)
```
The dispatcher for the histogram family above: each observable's own
`use_log_y` flag (from the catalogue) automatically picks the right
binning strategy — falling pT/mass spectra get log bins, bounded/peaked
histogram-drawn observables get quantile bins.

```python
def _norm_hist(vals, edges):
    counts, _ = np.histogram(vals, bins=edges)
    bw        = np.diff(edges)
    total     = counts.sum()
    density   = counts / (total * bw + 1e-12)
    err       = np.sqrt(counts) / (total * bw + 1e-12)
    cx        = 0.5 * (edges[:-1] + edges[1:])
    return cx, density, err
```
Turns raw per-bin counts into a normalized probability density — dividing
by `total * bin_width` (not just `total`) is what makes the histogram
correctly integrate to 1 even with variable-width bins. `np.sqrt(counts)`
is the standard Poisson uncertainty on a raw count, propagated through
the same normalization. The `1e-12` is a numerical-safety epsilon
against division by exactly zero.

`plt.rcParams.update({...})` sets a single, consistent visual style
(serif font, tick direction/size, DPI, `pdf.fonttype: 42` so PDF text
stays selectable/searchable rather than rasterized) applied automatically
to every figure — no per-plot styling code needed elsewhere. This global
style block, the colour choices (`C_SPS`/`C_DPS`), figure sizes, panel
layouts, legends, and output filenames/directory structure are all
unchanged from the project's original design and must stay that way.

`plot_panel(...)`, `make_process_single(...)`, `make_summary(...)`,
`make_eta_heatmap(...)` are the figure-drawing functions: given
already-loaded data and an observable's catalogue entry, they draw the
correct axis (KDE curve via `_kde_curve`/`gaussian_kde`, or normalized
histogram + error bars via `ax.stairs`/`ax.errorbar`, or the (η_W, η_Z)
2D joint-density heatmap via `ax.hist2d`), compute and annotate the KS
statistic (`scipy.stats.ks_2samp` — tests whether the SPS and DPS samples
for this observable could plausibly be drawn from the same underlying
distribution; a low p-value means "no, they're statistically
distinguishable"), and save to the correct path under `plots/`.

`main()` ties it together: discover available (SPS, DPS) file pairs for
matching N values, ask (or auto-pick, per `AUTO_MODE`) which pair(s),
which observables, and which plot types to produce, then loop over the
selection calling the drawing functions above.

---

## 5. `Convergence.py` — how observables change with sample size N

### Purpose

Where `Compare.py` looks at one fixed N, `Convergence.py` overlays the
**same** observable at several different N values (e.g. 10,000 through
several hundred thousand) on one set of axes, to answer: "has this
distribution's shape stopped changing as I add more Monte Carlo
statistics, or is what I'm seeing still just noise from too few events?"

### Shared machinery with `Compare.py`

`ALL_OBSERVABLES` (extended with the same eight `dPhi_*` entries
described in §4), the shared `load()`-style derivation of every
cross-particle observable from raw branches, `_adaptive_edges`,
`_log_edges`, `_choose_edges`, histogram helpers, and `safe_input` follow
the same reasoning as the identically-named pieces in `Compare.py` —
deliberately kept consistent so a reader who understands one script
already understands the corresponding piece of the other. The differences
worth calling out:

```python
def make_palette(n_items): ...
def _alpha(idx, total): ...
```
Since this script can overlay many curves (one per N value) on a single
axis, it needs a colour/transparency scheme that stays readable with
several N values plotted at once. `make_palette` assigns each N a
distinct colour from a perceptually uniform colormap; `_alpha` scales
transparency so overlapping curves — which happen a lot, since higher-N
curves should converge toward the same shape as lower-N ones — don't turn
into an illegible solid block. Low-N (noisiest) curves are the most
transparent; high-N (most trustworthy) curves are drawn last and fully
opaque, so convergence is visually obvious.

```python
def _dir_single(base, proc): ...
def _dir_multi(base, proc): ...
def _nsuf(n_list): ...
def _ind_path(base, proc, tag, n_list): ...
def _sum_path(base, proc, tags, n_list): ...
```
Output-path builders: since a convergence plot's filename depends on
*which set* of N values were plotted (not just one N), `_nsuf` builds a
filename suffix (e.g. `_N10000_50000_800000`) from whatever subset was
actually chosen, so different N-selections never silently overwrite each
other's output files.

```python
def _cache_key(out_path): ...
def _load_cache(): ...
def _save_cache(cache): ...
def _needs_replot(out_path, source_paths): ...
def _record(out_path, source_paths): ...
def _ask_replot(out_path, reason): ...
```
A small on-disk cache (keyed by output path, storing the modification
time and size of every source ROOT file that went into producing it) that
detects "has anything actually changed since I last produced this exact
plot?" If nothing changed, it asks before silently redoing (and
overwriting) work that would produce byte-identical output — useful when
convergence plots can be slow to regenerate across many N values and many
observables.

`draw_panel(...)`, `make_individual(...)`, `make_summary(...)`,
`print_stats_table(...)` are the drawing functions, following the same
histogram/error-bar/KDE logic as `Compare.py`'s drawing functions, but
looping over multiple `(N, dataset)` pairs instead of a fixed SPS/DPS
pair. `ask_mode_menu(...)` and `ask_n_selection(...)` are the interactive
menus specific to this script: which process (SPS, DPS, or both) and
which subset of available N values to include in the overlay.

---

## 6. `Analyze.py` — the full quantitative SPS-vs-DPS study

### Purpose

The most analysis-heavy script: rather than just overlaying
distributions, it runs a battery of statistical tests and derived-
observable studies meant to answer "how distinguishable are SPS and DPS,
quantitatively, and is that distinguishability robust or a statistics
artifact."

### `OBSERVABLES` and `CORR_PAIRS`

`OBSERVABLES` is the same kind of catalogue as the other two scripts
(kept independently per-file rather than imported from a shared module,
per the project's stated goal of not restructuring the existing
architecture), extended with the same eight lepton-level `dPhi_*` entries
described in §4 — and, like the other scripts, every non-raw entry in it
is derived in Python from raw branches rather than read pre-computed off
the tree. `CORR_PAIRS` is a separate list of raw kinematic branch pairs
(e.g. `(W_eta, Z_eta)`) used specifically by the 2D
correlation/joint-density analyses — a distinct list because those
functions inherently work on two variables at once, not one.

### `_paired_ns` / `discover` / `load`

Same file-discovery and raw-branch-loading-plus-derivation logic as
`Compare.py`, adapted to this script's needs (e.g. `_paired_ns`
specifically finds every N for which both an SPS and a DPS file exist,
since every analysis here needs a matched pair).

### The `analysis_*` functions — one per statistical test

Each is a self-contained "run this one test/plot, save its output, print
its numbers" unit:

- **`analysis_ratio`** — for each observable, plots the SPS and DPS
  normalized histograms plus their bin-by-bin ratio (SPS/DPS); a ratio
  flat at 1 means the two shapes agree in that region, deviations show
  exactly where and how strongly they differ.
- **`analysis_corr_matrix`** — computes and displays the full pairwise
  Pearson correlation matrix among the kinematic variables in
  `CORR_PAIRS`, separately for SPS and DPS — directly testing the core
  physics claim of the whole project (SPS bosons should be correlated by
  sharing one hard vertex; DPS bosons should be much closer to
  uncorrelated).
- **`analysis_2d`** — 2D joint-density heatmaps (the same idea as the
  η_W⊗η_Z heatmap in `Compare.py`, generalized to every pair in
  `CORR_PAIRS`) — the direct visual counterpart to the correlation
  matrix's single numbers.
- **`analysis_deltaR`** — `ΔR(W,Z) = √(Δη² + Δφ²)`, the standard collider
  angular-separation observable combining both angles into one number,
  built from the Python-derived `dEta_WZ`/`dPhi_WZ` (no longer ROOT
  branches — see §4).
- **`analysis_shape`** — skewness, kurtosis, and other distribution-shape
  statistics for each observable, a more quantitative alternative to
  eyeballing "does this histogram look different."
- **`analysis_ks_vs_n`** — runs the KS test (SPS vs DPS) separately at
  every available N, then plots the KS statistic vs. N — directly
  answering "does the observed SPS/DPS separation hold up as statistics
  increase, or does it wash out (meaning it was a low-statistics
  fluctuation all along)."
- **`analysis_fb`** — forward-backward asymmetry,
  `(N_forward − N_backward) / (N_forward + N_backward)`, computed on
  *signed* rapidity-type observables (this is exactly where the signed
  `dY_signed`/`dEta_signed` derivation matters — an asymmetry computed on
  an always-non-negative quantity would be meaningless by construction,
  since "backward" would never register any events).
- **`analysis_chi2`** — bin-by-bin χ² goodness-of-fit between the SPS and
  DPS histograms, on fixed-width bins deliberately (not the adaptive
  binning used for visual plots elsewhere), because a χ² test's meaning
  depends on its bin definition being fixed independently of the data
  being tested — adaptive bins built from the data itself would make the
  test circular.
- **`analysis_uniformity`** — the χ²-vs-flat test built to give a
  rigorous, numeric answer to "is DPS's Δφ distribution actually flat, or
  does it have real (if subtle) structure," using fixed bins and the
  folded `|Δφ|` observables. Its range dictionary now also covers every
  one of the lepton-level Δφ observables (`dPhi_lep1_lep2`, `dPhi_lep_nu`,
  `dPhi_W_lep1`, `dPhi_W_lep2`, `dPhi_Z_lep1`, `dPhi_Z_lep2`), not just
  `dPhi_WZ`/`dPhi_Z_lep`/`dPhi_Z_nu`, so the same rigorous flatness test
  applies to every angular observable the project tracks — all of them
  now computed in this script's own `load()`, per the current design.

### `ANALYSES`, menus, `main`

`ANALYSES` is the master menu list (letter code, description, dispatch
key) tying every `analysis_*` function to a human-readable menu entry;
`ask_analyses`/`ask_n_selection` are the interactive selection prompts;
`main()` loads the requested N's paired data and dispatches to whichever
`analysis_*` functions were selected.

---

## 7. `Makefile` — the build system

```makefile
MAKEFLAGS += --no-builtin-rules
.SUFFIXES:
```
Disables Make's built-in implicit rules (e.g. its default `%.o: %.cc`
compilation rule). Without this, a target name that happens to look like
a source filename could get silently intercepted by Make's own
pattern-matching before this Makefile's own explicit rule gets a chance
to run.

```makefile
$(BUILD)/.stamp:
	@mkdir -p $(BUILD)
	@touch $@
```
Binaries' compile rules depend on `$(BUILD)/.stamp` (an empty marker
file) as an order-only prerequisite for "the build directory must
exist," rather than depending on the literal directory name. This
sidesteps a real trap: if `BUILD=build` and there's also a phony target
literally named `build`, Make can't tell "depend on the directory
existing" from "depend on the phony target `build`," causing a circular
dependency warning.

```makefile
all:
	@chmod +x interactive_all.sh
	@./interactive_all.sh

all_auto: build
	...
	@python3 Compare.py --auto
```
`make all` (the default, interactive entry point) delegates entirely to
`interactive_all.sh`, which owns the actual interactive loop (§8).
`make all_auto` is the non-interactive equivalent used by `run_all.sh`
and any automated context: build, run both generators, then call
`Compare.py --auto` so nothing ever blocks waiting on a prompt.

```makefile
NSPS := $(shell grep -E 'constexpr\s+int\s+N_EVENTS_SPS\s*=' Settings.h | grep -oE '[0-9]+' | tail -1)
NDPS := $(shell grep -E 'constexpr\s+int\s+N_EVENTS_DPS\s*=' Settings.h | grep -oE '[0-9]+' | tail -1)
```
Make reads the current event count directly out of `Settings.h` at
parse time (via regex extraction, not by parsing C++), so every
binary/output filename automatically stays in sync with whatever
`settings_editor.py` last wrote there — no separate "remember to update N
here too" step.

The remaining targets (`build`, individual `GenerateSPS`/`GenerateDPS`
rules, `compare`/`convergence`/`analyze` Python-script shortcuts, `info`,
`clean*`) follow standard Make conventions: each binary's rule depends on
its source `.cc` file *and* `Settings.h`, so editing either automatically
triggers a rebuild next time; `clean_root`/`clean_plots`/`clean_all` give
explicit, separate ways to wipe generated ROOT files vs. plots vs.
everything.

---

## 8. `interactive_all.sh` — the `make all` driver

```bash
set -uo pipefail   # NOTE: no -e
```
Deliberately omits `-e` (which would abort the whole script on any
command's non-zero exit code): this script must survive a failed
submenu action (e.g. a failed compile) and return to its own menu loop
rather than dying outright. `-u` (error on unset variables) and
`pipefail` (a pipeline's exit status is its last *failing* command) are
kept, since those failure modes genuinely indicate bugs worth catching.

`do_build()` and `do_run()` are two shell functions wrapping the compile
step and the run step respectively, each independently callable —
`do_build` reports how many binaries were actually (re)compiled this
round vs. already up to date, using the same up-to-date timestamp check
the Makefile itself uses (`[ "$src" -nt "$bin" ]`, i.e. "is the source
newer than the binary").

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
The main interactive loop, deliberately never auto-building on launch (the
menu is shown first, before anything happens) — so simply typing
`make all` and immediately exiting never has an unwanted side effect of
silently compiling/running anything you didn't ask for. Option 2 (change
event counts) chains directly into a rebuild afterward, since a changed N
is meaningless until the binaries are recompiled against it.

---

## 9. `run_all.sh` — non-interactive pipeline

A simpler, linear version of the same compile → run → analyze pipeline,
with no menu at all: compile both generators if needed (same up-to-date
check as `interactive_all.sh`), run both **in parallel** (`&`
backgrounding + `wait`, since SPS and DPS generation are fully
independent of each other), then call `python3 Compare.py --auto` at the
end.

```bash
set -euo pipefail
```
Unlike the interactive script, this one keeps `-e` — appropriate here,
since this is meant to be run unattended (e.g. from a cron job or CI
pipeline) where a failure should immediately stop everything rather than
silently continuing.

```bash
if [ $# -ge 1 ]; then NSPS=$1; NDPS=$1; fi
if [ $# -ge 2 ]; then NSPS=$1; NDPS=$2; fi
```
Command-line overrides: no arguments uses whatever `Settings.h` currently
says; one argument sets both SPS and DPS event counts to that value; two
arguments sets them independently.

```bash
"$BUILD/GenerateSPS_${NSPS}" > "$LOG_SPS" 2>&1 &  PID_SPS=$!
"$BUILD/GenerateDPS_${NDPS}" > "$LOG_DPS" 2>&1 &  PID_DPS=$!

wait $PID_SPS && echo "..." || { echo "ERROR"; exit 1; }
wait $PID_DPS && echo "..." || { echo "ERROR"; exit 1; }
```
`$!` captures the PID of the most recently backgrounded job; `wait $PID`
waits for that specific process and returns its exit code. This is more
explicit than a bare `& wait` because it identifies exactly which
generator failed, if either does.

---

## 10. `settings_editor.py` — interactive `Settings.h` editor

### Why this exists

Editing `Settings.h` by hand is error-prone (a typo in the `constexpr`
syntax could silently fail to update the intended value, or worse, break
compilation). This script provides a numbered menu, reads and rewrites
`Settings.h` via regex, and keeps a `.bak` backup before every save.

```python
PARAMS = [
    ("Centre-of-mass energy",    "ECM_GEV",      float, "Beam energy",              "GeV"),
    ("Z mass window min",        "Z_MASS_MIN",   float, "Low edge of Z pole window","GeV"),
    ("Z mass window max",        "Z_MASS_MAX",   float, "High edge of Z pole window","GeV"),
    ("N events  (SPS)",          "N_EVENTS_SPS", int,   "Events for GenerateSPS",    ""),
    ("N events  (DPS)",          "N_EVENTS_DPS", int,   "Events for GenerateDPS",    ""),
    ("Enable MPI",               "ENABLE_MPI",   bool,  "Multi-parton interactions", ""),
]
```
The same "one table drives everything" pattern as the plotting scripts'
observable catalogues — each entry names a display label, the exact C++
token to search for in `Settings.h`, its Python type (used to pick the
right regex and the right input-parsing/validation logic), a description,
and a unit string for display.

```python
def get_value(content, token, typ):
    if typ == bool:
        m = re.search(rf"constexpr\s+bool\s+{token}\s*=\s*(true|false)", content)
        ...
```
Regex-based extraction rather than an actual C++ parser — deliberately
simple, matching the exact `constexpr <type> <TOKEN> = <value>;` pattern
this project's `Settings.h` always uses. `\s*` allows whitespace variation
(tabs vs. spaces, alignment) without requiring exact formatting.

```python
def set_value(content, token, typ, new_val):
    ...
    new_content, n = re.subn(pattern, repl, content)
    if n == 0:
        print(f"  WARNING: token '{token}' not found in Settings.h")
    return new_content
```
`re.subn` returns both the modified string and a count of substitutions
made; checking `n == 0` catches the case where the regex simply didn't
match anything (e.g. someone manually reformatted `Settings.h` in a way
this script's pattern no longer recognizes), warning rather than
silently doing nothing.

```python
def save(content):
    backup = SETTINGS_FILE + ".bak"
    shutil.copy2(SETTINGS_FILE, backup)
    with open(SETTINGS_FILE, "w") as f:
        f.write(content)
```
Always copies the current on-disk file to `.bak` immediately before
overwriting it — a simple, single-level undo mechanism (the `R` menu
command restores from this backup) in case a save turns out to be a
mistake. `shutil.copy2` preserves file metadata (timestamps), which is
relevant because the Makefile's rebuild logic is timestamp-based.

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
        cmd = input("\n  Command: ").strip().upper()
        ...
```
`--show` is a fast, non-interactive path (just print current values and
exit). The main loop re-reads `Settings.h` from disk at the top of
*every* iteration (not just once at startup), so if you use the `R`
(restore backup) command mid-session, the displayed values are always
accurate rather than a stale in-memory copy. The per-parameter edit flow
(`cmd` is a number → `prompt_value(...)` → `set_value(...)` → immediate
write to disk) writes each individual change to disk as soon as it's
confirmed — the single backup made at the very first `S`/`A` command
covers "what it looked like before this session started editing," which
is the meaningful restore point.

---

## Cross-file contract summary

**Protocol 1 — ROOT tree naming.** Both `SPS_<N>.root` and
`DPS_<N>.root` contain a single tree named `"Events"`.

**Protocol 2 — ROOT tree branch names (raw kinematics only).** Both
`Events` trees have an **identical** branch schema, and — unlike earlier
versions of this project — that schema now contains **only raw,
per-particle kinematics and the handful of event-record-only quantities
that must be computed inside PYTHIA**:

```
W_pt, W_eta, W_phi, W_mass, W_rapidity, W_px, W_py, W_pz, W_E, W_charge,
Z_pt, Z_eta, Z_phi, Z_mass, Z_rapidity, Z_px, Z_py, Z_pz, Z_E,
lep_pt..lep_E, nu_pt..nu_E,               (W's own lepton + neutrino)
lep1_pt..lep1_E, lep2_pt..lep2_E,         (Z's two decay leptons, μ-/μ+)
weight, nJets30.
```

There are no `dEta_WZ`, `dPhi_WZ`, `dY_WZ`, `M_WZ`, `pT_WZ`, `sT`,
`pT_imbalance`, or any `dPhi_*` lepton-pair branches in the ROOT file
anymore — every one of those is now a **derived Python quantity**,
computed inside each analysis script's own `load()` function from the
raw branches above (see §4). This is the engineering reason why
`GenerateSPS.cc` and `GenerateDPS.cc` produce the exact same, purely-raw
branch list: any Python script that reads one tree can read the other
with no file-type-specific logic, and all analysis-defining logic lives
in exactly one layer (Python), not split across C++ and Python.

**Protocol 3 — N encoding.** `Settings.h` defines the default N. The
Makefile reads it at make-time via regex. Binaries, ROOT files, logs, and
plots are all stamped with N in their filename, so every artifact of a
given run stays traceable and no two different-N runs ever collide on
disk.

**Protocol 4 — the Δφ convention.** Every angular separation used
anywhere in the project — `dPhi_WZ`, `dPhi_Z_lep`, `dPhi_Z_nu`,
`dPhi_lep1_lep2`, `dPhi_lep_nu`, `dPhi_W_lep1`, `dPhi_W_lep2`,
`dPhi_Z_lep1`, `dPhi_Z_lep2` — is folded to `[0, π]` using the same
convention everywhere: defined once in C++ as `dPhi()` in `Settings.h`
(now retained there only as the canonical reference definition, since the
generators no longer call it), and reimplemented as a NumPy-vectorized
helper (`_dphi_np()`) in each Python analysis script, which is where
every one of these observables is actually computed today. No script
recomputes an angular separation with different folding logic, so every
Δφ observable is directly comparable across every plot and every test in
the project.

**Adding a new observable — the three-step pattern.** Every observable
added anywhere in this project follows the same pattern: (1) if it needs
data not already stored, add the necessary raw per-particle kinematics as
a ROOT branch in both C++ generators (only genuinely PYTHIA-only
quantities, like a new event-record scan, belong in C++ — anything
derivable from existing raw branches should be derived in Python
instead); (2) add its derivation (if any) plus one entry to the
`ALL_OBSERVABLES`/`OBSERVABLES` catalogue in each Python script that
should plot it; (3) nothing else — every plotting function, menu,
filename builder, and statistical test reads from that one catalogue, so
no observable-specific code needs to be written anywhere else.
