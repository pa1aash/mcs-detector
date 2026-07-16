// PrimaryGeneratorAction --- single-particle pencil beam along +z (S3).
// Particle and kinetic energy come from the standard /gun/ macro commands
// (default: proton, 200 MeV). The beam starts upstream of the target in vacuum,
// so no scattering occurs before the target boundary. Optional Gaussian angular
// divergence via /mcs/gun/divergence (default 0 = pencil).
#ifndef MCS_PRIMARYGENERATORACTION_HH
#define MCS_PRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class G4GenericMessenger;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event* event) override;

  const G4ParticleGun* GetParticleGun() const { return fGun; }

 private:
  void DefineCommands();

  G4ParticleGun* fGun = nullptr;
  G4GenericMessenger* fMessenger = nullptr;
  G4double fZStart;       // beam origin z (upstream, in vacuum)
  G4double fDivergence;   // Gaussian sigma of the projected entry angle [rad]
  G4double fSpotXY = 0.;  // full width of uniform transverse sampling [mm]
                          // (0 = on-axis pencil; set = cell to sample one cell)
};

#endif  // MCS_PRIMARYGENERATORACTION_HH
