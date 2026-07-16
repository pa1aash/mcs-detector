// VoxelParam implementation (S5.5; non-cubic in S5-rebuilt).
#include "VoxelParam.hh"

#include <algorithm>
#include <fstream>
#include <sstream>

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ios.hh"

VoxelParam::VoxelParam(const G4String& path, G4Material* solid,
                       G4Material* voidMat)
    : fSolid(solid), fVoid(voidMat) {
  std::ifstream meta((path + ".meta").c_str());
  if (!meta) {
    G4Exception("VoxelParam", "vox001", FatalException,
                ("missing voxel meta sidecar: " + path + ".meta").c_str());
  }
  // Read the whole line and parse 5 tokens (Nx Ny Nz voxel cell) or legacy 3
  // tokens (N voxel cell -> cubic).
  std::string line;
  std::getline(meta, line);
  std::istringstream is(line);
  std::vector<double> tok;
  double v;
  while (is >> v) tok.push_back(v);
  G4double voxel_mm = 0., cell_mm = 0.;
  if (tok.size() >= 5) {
    fNx = (int)tok[0]; fNy = (int)tok[1]; fNz = (int)tok[2];
    voxel_mm = tok[3]; cell_mm = tok[4];
  } else if (tok.size() >= 3) {
    fNx = fNy = fNz = (int)tok[0];
    voxel_mm = tok[1]; cell_mm = tok[2];
  } else {
    G4Exception("VoxelParam", "vox004", FatalException,
                "voxel meta needs 'Nx Ny Nz voxel cell' or 'N voxel cell'");
  }
  fVoxel = voxel_mm * mm;
  fCell = cell_mm * mm;
  fOx = -0.5 * fNx * fVoxel;
  fOy = -0.5 * fNy * fVoxel;
  fOz = -0.5 * fNz * fVoxel;

  const size_t ntot = (size_t)fNx * fNy * fNz;
  fField.resize(ntot);
  std::ifstream raw(path.c_str(), std::ios::binary);
  if (!raw) {
    G4Exception("VoxelParam", "vox002", FatalException,
                ("missing voxel field: " + path).c_str());
  }
  raw.read(reinterpret_cast<char*>(fField.data()), (std::streamsize)ntot);
  if ((size_t)raw.gcount() != ntot) {
    G4Exception("VoxelParam", "vox003", FatalException,
                "voxel field shorter than Nx*Ny*Nz");
  }
  G4cout << "[VoxelParam] " << fNx << "x" << fNy << "x" << fNz
         << " voxels, solid frac="
         << (double)std::count(fField.begin(), fField.end(), (unsigned char)1) /
                (double)ntot
         << G4endl;
}

void VoxelParam::ComputeTransformation(const G4int copyNo,
                                       G4VPhysicalVolume* physVol) const {
  const G4int ix = copyNo / (fNy * fNz);
  const G4int iy = (copyNo / fNz) % fNy;
  const G4int iz = copyNo % fNz;
  physVol->SetTranslation(G4ThreeVector(fOx + (ix + 0.5) * fVoxel,
                                        fOy + (iy + 0.5) * fVoxel,
                                        fOz + (iz + 0.5) * fVoxel));
}

G4Material* VoxelParam::ComputeMaterial(const G4int copyNo,
                                        G4VPhysicalVolume* /*physVol*/,
                                        const G4VTouchable* /*parentTouch*/) {
  return fField[(size_t)copyNo] ? fSolid : fVoid;
}
