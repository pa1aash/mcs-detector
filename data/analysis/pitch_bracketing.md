# Δκ4-vs-pitch transport bracketing (CONV-PITCH, transport level; Fable 4.4)

Full Geant4 at f=0.40, 200 MeV, campaign statistics (rect 3e6, gyroid 6e6), locked
physics/windowing/all-order floor. Regenerated on the pod 2026-07-16 (Geant4 11.3.2).
Bracketing the operating pitch (rect spc=48/52µm; gyroid spc=32/78µm) with a coarser and a
finer voxel field. `analysis`: pod_task3.py (build field at res -> mcs_sim -> deconvolved Δκ4).

| config | pitch | Δκ4 [rad^4] (95% CI) | γ2(raw,in-window) | k_eff |
|---|---|---|--:|--:|
| rect spc40 | 62 µm | 5.204e-10 [4.937,5.470]e-10 | 7.98 | 1.015 |
| **rect spc48 (operating)** | **52 µm** | **5.170e-10 [4.907,5.433]e-10** | 8.07 | 1.009 |
| rect spc80 | 31 µm | 5.428e-10 [5.161,5.695]e-10 | 7.97 | 1.055 |
| gyro spc24 | 104 µm | 2.670e-10 [2.481,2.859]e-10 | 7.05 | 0.947 |
| **gyro spc32 (operating)** | **78 µm** | **2.779e-10 [2.583,2.975]e-10** | 7.09 | 1.005 |
| gyro spc48 | 52 µm | 2.707e-10 [2.529,2.884]e-10 | 7.08 | 0.975 |

**Verdict: PASS.** For both topologies the deconvolved Δκ4 at the three pitches have mutually
overlapping 95% bootstrap CIs (rect spread ≤5%, gyroid ≤4%); γ2 is flat across pitch. The
Geant4-transport Δκ4, not only the geometric input, is pitch-converged at the campaign
operating resolution. Closes CONV-PITCH at the transport level.
