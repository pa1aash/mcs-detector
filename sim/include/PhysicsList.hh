// PhysicsList --- LOCKED configuration (S3; ROADMAP_SOURCE.md s2.4):
//   FTFP_BERT hadronic  +  G4EmStandardPhysics_option4 EM  (which provides the
//   WentzelVI multiple-scattering + G4eCoulombScattering single-scattering split
//   for protons/hadrons)  +  G4StepLimiterPhysics (so the in-target G4UserLimits
//   max step is enforced).
// Default production range cut 0.05 mm (overridable by /run/setCut).
#ifndef MCS_PHYSICSLIST_HH
#define MCS_PHYSICSLIST_HH

#include "FTFP_BERT.hh"
#include "globals.hh"

class PhysicsList : public FTFP_BERT {
 public:
  explicit PhysicsList(G4int verbose = 1);
  ~PhysicsList() override = default;

  void SetCuts() override;
};

#endif  // MCS_PHYSICSLIST_HH
