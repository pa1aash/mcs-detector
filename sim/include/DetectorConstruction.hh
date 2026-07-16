// DetectorConstruction --- vacuum world + a macro-selectable target (S5.5).
//
//   /mcs/det/geom slab   : homogeneous PLA G4Box (S3 validated path; DEFAULT,
//                          unchanged -- the Highland-regression baseline).
//   /mcs/det/geom stl    : an ASCII-STL unit cell loaded via vendored CADMesh
//                          into a G4TessellatedSolid, tiled transversely
//                          (ntile x ntile, periodic) and along z (depth/cell)
//                          to build the L-deep lattice the ray-tracer scored.
//   /mcs/det/geom voxel  : a binary voxel field (Voronoi) as a G4 nested
//                          parameterisation (PLA / vacuum per voxel).
//   /mcs/det/reflective 1: single transverse cell with specular transverse
//                          walls (SteppingAction fold) -- the cheap deep-cell
//                          option, validated against true multi-cell in S5.5.
//
// Cuts / maxStep / material stay exactly as S3 (macro params). Topology / cell /
// infill tags are carried into the metadata sidecar for provenance.
#ifndef MCS_DETECTORCONSTRUCTION_HH
#define MCS_DETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"  // G4RotationMatrix is a CLHEP typedef, not a class
#include "globals.hh"

class G4VPhysicalVolume;
class G4Material;
class G4GenericMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction {
 public:
  DetectorConstruction();
  ~DetectorConstruction() override;

  G4VPhysicalVolume* Construct() override;

  // Accessors used by RunAction (sidecar) and SteppingAction (scoring/reflection).
  const G4VPhysicalVolume* GetTargetPV() const { return fTargetPV; }
  G4double GetThickness() const { return fThickness; }
  G4double GetMaxStep() const { return fMaxStep; }
  G4String GetMaterialName() const { return fMaterialName; }
  G4double GetTargetX0() const { return fTargetX0; }            // mm, GetRadlen()
  G4double GetTargetDensity() const { return fTargetDensity; }  // g/cm3

  // Geometry-agnostic kink scoring: planes at z = -/+ 0.5*depth.
  G4double GetTargetDepth() const { return fTargetDepth; }      // mm
  // Reflective single-cell support (SteppingAction).
  G4bool GetReflective() const { return fReflective; }
  G4double GetCellSize() const { return fCellSize; }

  // M1 telescope: 6 scoring planes at z = {-300,-200,-100,+100,+200,+300} mm,
  // gated OFF by default (campaign path unchanged). SteppingAction records the
  // primary's (x,y) crossing at each plane; EventAction fills the "telescope"
  // ntuple. Optional 50 um Si planes add realistic telescope scattering.
  static constexpr G4int kNPlanes = 6;
  G4bool GetTelescope() const { return fTelescope; }
  G4int GetNPlanes() const { return kNPlanes; }
  G4double GetPlaneZ(G4int i) const;                            // mm (G4 units)

  // M2 incidence-angle tilt: rotate the target about x by fTilt (rad). Provenance getter
  // returns degrees for the meta sidecar.
  G4double GetTiltDeg() const;
  G4double GetTilt() const { return fTilt; }   // rad (used by SteppingAction scoring)

  // Provenance for the sidecar.
  G4String GetGeomType() const { return fGeomType; }
  G4String GetTopology() const { return fTopology; }
  G4double GetInfill() const { return fInfill; }
  G4int GetNTile() const { return fNTile; }

 private:
  void DefineMaterials();
  void DefineCommands();
  G4Material* ResolveMaterial(const G4String& name);

  G4VPhysicalVolume* BuildSlab(G4LogicalVolume* world, G4Material* target);
  G4VPhysicalVolume* BuildSTL(G4LogicalVolume* world, G4Material* target);
  G4VPhysicalVolume* BuildVoxel(G4LogicalVolume* world, G4Material* target);
  G4RotationMatrix* MakeTiltRot() const;   // M2: target rotation about x (nullptr if 0)

  G4GenericMessenger* fMessenger = nullptr;
  G4VPhysicalVolume* fTargetPV = nullptr;

  // S3 macro params (defaults = headline Highland-check point).
  G4double fThickness;     // slab depth along beam (z)
  G4double fMaxStep;       // G4UserLimits max step in the target
  G4String fMaterialName;  // "PLA" or any NIST name

  // S5.5 lattice params.
  G4String fGeomType = "slab";   // slab | stl | voxel
  G4String fStlPath = "";        // unit-cell ASCII STL
  G4String fVoxelPath = "";      // .raw binary voxel field (uint8) for voxel mode
  G4double fCellSize = 2.5 * 1.0;  // mm (set in ctor with units)
  G4double fLatticeDepth = 10.0 * 1.0;  // mm, z extent of the lattice target
  G4int fNTile = 3;              // transverse cells per side (periodic tiling)
  G4bool fReflective = false;    // single transverse cell + specular walls
  G4String fTopology = "";       // provenance tag
  G4double fInfill = 0.;         // provenance tag
  G4int fVoxelN = 0;             // voxel grid edge count (cubic), voxel mode

  // M1 telescope.
  G4bool fTelescope = false;     // enable 6-plane telescope scoring
  G4bool fPlaneSi = false;       // place 50 um Si planes (else massless/virtual)

  // M2 incidence-angle tilt (radians; set via /mcs/det/tilt in deg).
  G4double fTilt = 0.;

  // Filled at Construct() time.
  G4double fTargetX0 = 0.;
  G4double fTargetDensity = 0.;
  G4double fTargetDepth = 0.;    // mm; = fThickness (slab) or fLatticeDepth
};

#endif  // MCS_DETECTORCONSTRUCTION_HH
