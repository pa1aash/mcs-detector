# `sim/` — Geant4 application (built S1)

Single-particle gun fired through a target slab (a topology geometry from
`geom/`) into downstream scoring planes, recording the **projected MCS kink
angle** per event. This is the transport engine behind E0, the full campaign,
and the MSC systematic (E5).

> **Status: stubs only (S0).** The files below are empty placeholders describing
> intended structure. S1 implements them and gets the app building + running a
> single track through a slab. Nothing here compiles yet.

## Intended structure

| File | Role |
|---|---|
| `CMakeLists.txt` | Standard Geant4 CMake project; finds Geant4, builds the `mcs_sim` executable, installs macros. |
| `src/mcs_sim.cc` | `main()` — run manager, UI/vis, macro execution. |
| `include/DetectorConstruction.hh` / `src/DetectorConstruction.cc` | World, target volume (loads STL/GDML topology), upstream + downstream scoring planes, optional telescope/air to mimic the measurement chain. |
| `include/PhysicsList.hh` / `src/PhysicsList.cc` | **FTFP_BERT + EmStandard_opt4 + WentzelVI + single scattering**, tightened production cuts, small max step — the tail-trustworthy configuration (`ROADMAP_SOURCE.md` §2.4). Switchable single-scattering reference for E5. |
| `include/PrimaryGeneratorAction.hh` / `src/PrimaryGeneratorAction.cc` | Particle gun; momentum set from macro (the 3 campaign energies). |
| `include/SteppingAction.hh` / `src/SteppingAction.cc` + scoring | Record in/out directions → projected kink angle; traversed line-integral `t` per event where available. ROOT or flat-file output to `data/runs/<tag>/`. |
| `macros/vis.mac` | Interactive visualization. |
| `macros/run.mac` | Batch template (energy, geometry, N events, output tag). |

## Physics configuration (locked, do not drift)

`FTFP_BERT` hadronic + `G4EmStandardPhysics_option4` EM, with `WentzelVI` +
single-scattering MSC, tightened cuts, small `G4UserLimits` max step. This is
what lets a referee trust the tails (`ROADMAP_SOURCE.md` §2.4, §12). The
single-scattering-only reference is a **one-point** systematic (E5) — expensive,
do not grid.

## Output contract

Per run, write to `data/runs/<tag>/`: one record per event with the projected
kink angles (x, y) and, where scored, the traversed material line-integral `t`.
The analysis library (`analysis/lib/`) consumes this.
