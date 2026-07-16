// DetectorConstruction implementation (S5.5: macro-selectable slab/stl/voxel).
#include "DetectorConstruction.hh"

#include <cmath>

#include "CADMesh.hh"
#include "G4Box.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4TessellatedSolid.hh"
#include "G4ThreeVector.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"
#include "VoxelParam.hh"

DetectorConstruction::DetectorConstruction()
    : fThickness(16. * mm),
      fMaxStep(0.1 * mm),
      fMaterialName("PLA") {
  fCellSize = 2.5 * mm;
  fLatticeDepth = 10. * mm;
  DefineMaterials();
  DefineCommands();
}

DetectorConstruction::~DetectorConstruction() { delete fMessenger; }

void DetectorConstruction::DefineMaterials() {
  auto nist = G4NistManager::Instance();
  if (!G4Material::GetMaterial("PLA", false)) {
    auto pla = new G4Material("PLA", 1.24 * g / cm3, 3);
    pla->AddElement(nist->FindOrBuildElement("C"), 3);
    pla->AddElement(nist->FindOrBuildElement("H"), 4);
    pla->AddElement(nist->FindOrBuildElement("O"), 2);
  }
  nist->FindOrBuildMaterial("G4_Galactic");
  nist->FindOrBuildMaterial("G4_AIR");
  nist->FindOrBuildMaterial("G4_Si");   // M1 telescope planes (optional)
  // M4: amorphous carbon at the detector-foam ligament density (rho = 1.7 g/cm3).
  if (!G4Material::GetMaterial("C_amorph", false)) {
    auto c = new G4Material("C_amorph", 1.7 * g / cm3, 1);
    c->AddElement(nist->FindOrBuildElement("C"), 1);
  }
}

G4double DetectorConstruction::GetPlaneZ(G4int i) const {
  static const G4double z[kNPlanes] = {-300. * mm, -200. * mm, -100. * mm,
                                       100. * mm, 200. * mm, 300. * mm};
  return (i >= 0 && i < kNPlanes) ? z[i] : 0.;
}

G4RotationMatrix* DetectorConstruction::MakeTiltRot() const {
  if (std::abs(fTilt) < 1e-9) return nullptr;      // no tilt -> unrotated placement
  auto rot = new G4RotationMatrix();
  rot->rotateX(fTilt);                             // tilt the target about the x axis
  return rot;
}

G4double DetectorConstruction::GetTiltDeg() const { return fTilt / deg; }

G4Material* DetectorConstruction::ResolveMaterial(const G4String& name) {
  if (name == "PLA") return G4Material::GetMaterial("PLA");
  auto mat = G4NistManager::Instance()->FindOrBuildMaterial(name);
  if (!mat) {
    G4Exception("DetectorConstruction::ResolveMaterial", "mat001",
                FatalException, ("unknown material: " + name).c_str());
  }
  return mat;
}

G4VPhysicalVolume* DetectorConstruction::Construct() {
  auto vacuum = G4Material::GetMaterial("G4_Galactic");
  auto target = ResolveMaterial(fMaterialName);
  fTargetX0 = target->GetRadlen();
  fTargetDensity = target->GetDensity() / (g / cm3);

  const G4double worldXY = 300. * mm, worldZ = 500. * mm;
  auto worldSolid = new G4Box("World", worldXY, worldXY, worldZ);
  auto worldLV = new G4LogicalVolume(worldSolid, vacuum, "World");
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr,
                                   false, 0, true);

  if (fGeomType == "slab") {
    BuildSlab(worldLV, target);
  } else if (fGeomType == "stl") {
    BuildSTL(worldLV, target);
  } else if (fGeomType == "voxel") {
    BuildVoxel(worldLV, target);
  } else {
    G4Exception("DetectorConstruction::Construct", "geom001", FatalException,
                ("unknown /mcs/det/geom: " + fGeomType).c_str());
  }

  // M1 telescope: optionally place 50 um Si scoring planes (else the crossings are
  // recorded geometrically in vacuum -- massless planes that isolate target physics).
  if (fTelescope && fPlaneSi) {
    auto si = G4Material::GetMaterial("G4_Si");
    auto planeS = new G4Box("SiPlane", 50. * mm, 50. * mm, 0.025 * mm);  // 50 um thick
    auto planeLV = new G4LogicalVolume(planeS, si, "SiPlane");
    auto pvis = new G4VisAttributes(G4Colour(0.9, 0.4, 0.2, 0.4));
    pvis->SetForceSolid(true);
    planeLV->SetVisAttributes(pvis);
    for (G4int i = 0; i < kNPlanes; ++i)
      new G4PVPlacement(nullptr, {0., 0., GetPlaneZ(i)}, planeLV, "SiPlane",
                        worldLV, false, i, true);
  }

  G4cout << "[DetectorConstruction] geom=" << fGeomType
         << "  material=" << fMaterialName << "  X0=" << fTargetX0 / mm << " mm"
         << "  depth=" << fTargetDepth / mm << " mm"
         << "  cell=" << fCellSize / mm << " mm  ntile=" << fNTile
         << "  reflective=" << (fReflective ? "on" : "off")
         << "  maxStep=" << fMaxStep / mm << " mm" << G4endl;
  return worldPV;
}

// ---- slab: the S3-validated homogeneous G4Box (UNCHANGED) -----------------
G4VPhysicalVolume* DetectorConstruction::BuildSlab(G4LogicalVolume* worldLV,
                                                   G4Material* target) {
  fTargetDepth = fThickness;
  const G4double tgtXY = 200. * mm;
  auto tgtSolid = new G4Box("Target", tgtXY, tgtXY, 0.5 * fThickness);
  auto tgtLV = new G4LogicalVolume(tgtSolid, target, "Target");
  tgtLV->SetUserLimits(new G4UserLimits(fMaxStep));
  auto vis = new G4VisAttributes(G4Colour(0.2, 0.6, 0.9, 0.3));
  vis->SetForceSolid(true);
  tgtLV->SetVisAttributes(vis);
  fTargetPV = new G4PVPlacement(MakeTiltRot(), {}, tgtLV, "Target", worldLV, false,
                                0, true);
  return fTargetPV;
}

// ---- stl: CADMesh unit cell, tiled transversely (periodic) + along z -------
G4VPhysicalVolume* DetectorConstruction::BuildSTL(G4LogicalVolume* worldLV,
                                                  G4Material* target) {
  if (fStlPath.empty()) {
    G4Exception("DetectorConstruction::BuildSTL", "geom002", FatalException,
                "/mcs/det/geom stl needs /mcs/det/stl <file>");
  }
  fTargetDepth = fLatticeDepth;
  const G4int ntxy = fReflective ? 1 : fNTile;     // transverse cells per side
  const G4int nz = std::max(1, (int)std::lround(fLatticeDepth / fCellSize));
  auto vacuum = G4Material::GetMaterial("G4_Galactic");

  // Container box exactly enclosing the tiled lattice (transverse = ntxy cells,
  // z = nz cells). Reflective walls (S5.5) are the transverse faces of THIS box.
  const G4double halfXY = 0.5 * ntxy * fCellSize;
  const G4double halfZ = 0.5 * nz * fCellSize;
  auto contSolid = new G4Box("LatContainer", halfXY, halfXY, halfZ);
  auto contLV = new G4LogicalVolume(contSolid, vacuum, "LatContainer");
  contLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  fTargetPV = new G4PVPlacement(MakeTiltRot(), {}, contLV, "LatContainer", worldLV,
                                false, 0, true);

  // One shared tessellated solid from the unit-cell STL (facet memory = 1 cell).
  auto mesh = CADMesh::TessellatedMesh::FromSTL(fStlPath);
  auto cellSolid = mesh->GetTessellatedSolid();
  auto cellLV = new G4LogicalVolume(cellSolid, target, "Cell");
  cellLV->SetUserLimits(new G4UserLimits(fMaxStep));
  auto vis = new G4VisAttributes(G4Colour(0.2, 0.6, 0.9, 0.35));
  vis->SetForceSolid(true);
  cellLV->SetVisAttributes(vis);

  // Place copies on the cell grid, block centred at the origin.
  const G4double x0 = -0.5 * (ntxy - 1) * fCellSize;
  const G4double z0 = -0.5 * (nz - 1) * fCellSize;
  G4int copy = 0;
  for (G4int iz = 0; iz < nz; ++iz)
    for (G4int ix = 0; ix < ntxy; ++ix)
      for (G4int iy = 0; iy < ntxy; ++iy) {
        G4ThreeVector pos(x0 + ix * fCellSize, x0 + iy * fCellSize,
                          z0 + iz * fCellSize);
        new G4PVPlacement(nullptr, pos, cellLV, "Cell", contLV, false, copy++,
                          false);
      }
  G4cout << "[BuildSTL] " << fStlPath << "  facets="
         << cellSolid->GetNumberOfFacets() << "  placements=" << copy
         << " (" << ntxy << "x" << ntxy << "x" << nz << ")" << G4endl;
  return fTargetPV;
}

// ---- voxel: binary field (Voronoi) as a flat G4 parameterisation -----------
G4VPhysicalVolume* DetectorConstruction::BuildVoxel(G4LogicalVolume* worldLV,
                                                    G4Material* target) {
  if (fVoxelPath.empty()) {
    G4Exception("DetectorConstruction::BuildVoxel", "geom003", FatalException,
                "/mcs/det/geom voxel needs /mcs/det/voxel <file.raw>");
  }
  auto vacuum = G4Material::GetMaterial("G4_Galactic");
  auto param = new VoxelParam(fVoxelPath, target, vacuum);  // reads Nx,Ny,Nz+field
  fVoxelN = param->GetNz();
  const G4double bx = param->GetBoxX(), by = param->GetBoxY(), bz = param->GetBoxZ();
  fTargetDepth = bz;
  fCellSize = param->GetCellSize();

  // Mother box holding the Nx*Ny*Nz parameterised voxels (centred at origin).
  auto motherS = new G4Box("VoxMother", 0.5 * bx, 0.5 * by, 0.5 * bz);
  auto motherLV = new G4LogicalVolume(motherS, vacuum, "VoxMother");
  motherLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  fTargetPV = new G4PVPlacement(MakeTiltRot(), {}, motherLV, "VoxMother", worldLV,
                                false, 0, true);

  const G4double h = 0.5 * param->GetVoxelSize();
  auto voxS = new G4Box("Voxel", h, h, h);
  auto voxLV = new G4LogicalVolume(voxS, target, "Voxel");  // material per-copy
  voxLV->SetUserLimits(new G4UserLimits(fMaxStep));
  new G4PVParameterised("Voxels", voxLV, motherLV, kUndefined,
                        param->GetNtotal(), param);
  G4cout << "[BuildVoxel] " << fVoxelPath << "  " << param->GetNx() << "x"
         << param->GetNy() << "x" << param->GetNz() << "  voxel="
         << param->GetVoxelSize() / mm << " mm  box=" << bx / mm << "x" << by / mm
         << "x" << bz / mm << " mm" << G4endl;
  return fTargetPV;
}

void DetectorConstruction::DefineCommands() {
  fMessenger = new G4GenericMessenger(this, "/mcs/det/", "target geometry");
  fMessenger->DeclarePropertyWithUnit("thickness", "mm", fThickness,
                                      "slab depth along the beam (z)");
  fMessenger->DeclarePropertyWithUnit("maxStep", "mm", fMaxStep,
                                      "max step length inside the target");
  fMessenger->DeclareProperty("material", fMaterialName,
                              "target material: PLA or any NIST name");
  fMessenger->DeclareProperty("geom", fGeomType, "slab | stl | voxel");
  fMessenger->DeclareProperty("stl", fStlPath, "unit-cell ASCII STL path");
  fMessenger->DeclareProperty("voxel", fVoxelPath, "binary voxel .raw path");
  fMessenger->DeclarePropertyWithUnit("cell", "mm", fCellSize,
                                      "lattice unit-cell edge");
  fMessenger->DeclarePropertyWithUnit("depth", "mm", fLatticeDepth,
                                      "lattice target depth along z");
  fMessenger->DeclareProperty("ntile", fNTile,
                              "transverse cells per side (periodic tiling)");
  fMessenger->DeclareProperty("reflective", fReflective,
                              "single transverse cell + specular walls");
  fMessenger->DeclareProperty("topology", fTopology, "provenance tag");
  fMessenger->DeclareProperty("infill", fInfill, "provenance tag");
  fMessenger->DeclareProperty("telescope", fTelescope,
                              "M1: enable 6-plane telescope scoring (default off)");
  fMessenger->DeclareProperty("planeSi", fPlaneSi,
                              "M1: place 50 um Si planes (else massless recording)");
  fMessenger->DeclarePropertyWithUnit("tilt", "deg", fTilt,
                                      "M2: target incidence-angle tilt about x");
}
