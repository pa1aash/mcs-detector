// PrimaryGeneratorAction implementation (S3).
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4GenericMessenger.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : fZStart(-50. * mm), fDivergence(0.) {
  fGun = new G4ParticleGun(1);
  auto proton = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  fGun->SetParticleDefinition(proton);
  fGun->SetParticleEnergy(200. * MeV);  // kinetic; override via /gun/energy
  fGun->SetParticlePosition({0., 0., fZStart});
  fGun->SetParticleMomentumDirection({0., 0., 1.});
  DefineCommands();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete fGun;
  delete fMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  // Beam on the z axis (or uniformly across one transverse cell if fSpotXY>0,
  // to sample the lattice). Direction +z (+ optional small Gaussian divergence).
  G4double bx = 0., by = 0.;
  if (fSpotXY > 0.) {
    bx = (G4UniformRand() - 0.5) * fSpotXY;
    by = (G4UniformRand() - 0.5) * fSpotXY;
  }
  fGun->SetParticlePosition({bx, by, fZStart});
  if (fDivergence > 0.) {
    const G4double ax = G4RandGauss::shoot(0., fDivergence);
    const G4double ay = G4RandGauss::shoot(0., fDivergence);
    G4ThreeVector dir(std::tan(ax), std::tan(ay), 1.0);
    fGun->SetParticleMomentumDirection(dir.unit());
  } else {
    fGun->SetParticleMomentumDirection({0., 0., 1.});
  }
  fGun->GeneratePrimaryVertex(event);
}

void PrimaryGeneratorAction::DefineCommands() {
  fMessenger =
      new G4GenericMessenger(this, "/mcs/gun/", "primary beam controls");
  fMessenger->DeclarePropertyWithUnit("zStart", "mm", fZStart,
                                      "beam origin z (upstream, in vacuum)");
  fMessenger->DeclareProperty(
      "divergence", fDivergence,
      "Gaussian sigma of the projected entry angle [rad] (0 = pencil)");
  fMessenger->DeclarePropertyWithUnit(
      "spotXY", "mm", fSpotXY,
      "full width of uniform transverse beam sampling (0 = on-axis pencil; "
      "set = cell to sample one lattice cell)");
}
