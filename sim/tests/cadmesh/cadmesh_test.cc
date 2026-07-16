// cadmesh_test.cc --- S3.5 dependency gate: load an STL via CADMesh into a
// G4TessellatedSolid, check volume, fire 100 primaries with no navigation error.
#include <cstdlib>

#include "CADMesh.hh"
#include "FTFP_BERT.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4TessellatedSolid.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "Randomize.hh"
#include "globals.hh"

static G4String gStl;

class Det : public G4VUserDetectorConstruction {
 public:
  G4VPhysicalVolume* Construct() override {
    auto nist = G4NistManager::Instance();
    auto vac = nist->FindOrBuildMaterial("G4_Galactic");
    auto water = nist->FindOrBuildMaterial("G4_WATER");
    auto worldS = new G4Box("World", 0.5 * m, 0.5 * m, 0.5 * m);
    auto worldL = new G4LogicalVolume(worldS, vac, "World");
    auto worldP =
        new G4PVPlacement(nullptr, {}, worldL, "World", nullptr, false, 0, true);

    auto mesh = CADMesh::TessellatedMesh::FromSTL(gStl);
    auto solid = mesh->GetTessellatedSolid();  // SetSolidClosed(true) internally
    G4cout << "[cadmesh] facets = " << solid->GetNumberOfFacets()
           << "  GetCubicVolume = " << solid->GetCubicVolume() / mm3 << " mm3"
           << G4endl;
    auto cubeL = new G4LogicalVolume(solid, water, "Cube");
    new G4PVPlacement(nullptr, {}, cubeL, "Cube", worldL, false, 0, true);
    return worldP;
  }
};

class Gun : public G4VUserPrimaryGeneratorAction {
 public:
  Gun() {
    fGun = new G4ParticleGun(1);
    fGun->SetParticleDefinition(
        G4ParticleTable::GetParticleTable()->FindParticle("proton"));
    fGun->SetParticleEnergy(200 * MeV);
    fGun->SetParticlePosition({0., 0., -50. * mm});
    fGun->SetParticleMomentumDirection({0., 0., 1.});
  }
  ~Gun() override { delete fGun; }
  void GeneratePrimaries(G4Event* e) override {
    // small lateral spread so rays sample faces/edges of the cube
    fGun->SetParticlePosition({(G4UniformRand() - 0.5) * 8. * mm,
                               (G4UniformRand() - 0.5) * 8. * mm, -50. * mm});
    fGun->GeneratePrimaryVertex(e);
  }

 private:
  G4ParticleGun* fGun;
};

class AI : public G4VUserActionInitialization {
 public:
  void Build() const override { SetUserAction(new Gun); }
};

int main(int argc, char** argv) {
  if (argc < 2) {
    G4cerr << "usage: cadmesh_test <file.stl>" << G4endl;
    return 2;
  }
  gStl = argv[1];
  auto rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  rm->SetUserInitialization(new Det);
  rm->SetUserInitialization(new FTFP_BERT(0));
  rm->SetUserInitialization(new AI);
  rm->Initialize();
  rm->BeamOn(100);
  delete rm;
  G4cout << "[cadmesh_test] completed 100 primaries" << G4endl;
  return 0;
}
