// VoxelParam --- flat G4 parameterisation of a binary (PLA/vacuum) voxel field
// (S5.5; generalised to non-cubic Nx*Ny*Nz in S5-rebuilt for L-deep validation
// blocks). Reads the field from <path> (uint8, 0=void 1=solid) plus a
// <path>.meta text sidecar: either "Nx Ny Nz voxel_mm cell_mm" (5 tokens) or the
// legacy cubic "N voxel_mm cell_mm" (3 tokens). Each parameterised copy is one
// cubic voxel; material = PLA where solid, vacuum otherwise.
#ifndef MCS_VOXELPARAM_HH
#define MCS_VOXELPARAM_HH

#include <vector>

#include "G4VPVParameterisation.hh"
#include "globals.hh"

class G4Material;
class G4VPhysicalVolume;
class G4Box;

class VoxelParam : public G4VPVParameterisation {
 public:
  VoxelParam(const G4String& path, G4Material* solid, G4Material* voidMat);
  ~VoxelParam() override = default;

  void ComputeTransformation(const G4int copyNo,
                             G4VPhysicalVolume* physVol) const override;
  G4Material* ComputeMaterial(const G4int copyNo, G4VPhysicalVolume* physVol,
                              const G4VTouchable* parentTouch = nullptr) override;

  G4int GetNx() const { return fNx; }
  G4int GetNy() const { return fNy; }
  G4int GetNz() const { return fNz; }
  G4int GetNtotal() const { return fNx * fNy * fNz; }
  G4double GetVoxelSize() const { return fVoxel; }
  G4double GetCellSize() const { return fCell; }
  G4double GetBoxX() const { return fNx * fVoxel; }
  G4double GetBoxY() const { return fNy * fVoxel; }
  G4double GetBoxZ() const { return fNz * fVoxel; }

 private:
  G4int fNx = 0, fNy = 0, fNz = 0;
  G4double fVoxel = 0.;   // voxel edge (mm)
  G4double fCell = 0.;    // provenance: lattice cell (mm)
  G4double fOx = 0., fOy = 0., fOz = 0.;  // block centred at origin
  G4Material* fSolid = nullptr;
  G4Material* fVoid = nullptr;
  std::vector<unsigned char> fField;  // ordering ix*Ny*Nz + iy*Nz + iz
};

#endif  // MCS_VOXELPARAM_HH
