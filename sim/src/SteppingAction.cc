// SteppingAction implementation (S5.5).
//
// Geometry-agnostic kink scoring by z-plane crossing (planes at z = -/+ 0.5*depth):
// identical to the S3 slab-face scoring for the slab (Highland-regression
// preserved) and well-defined for the tiled/voxel lattice (where the primary
// enters/leaves material many times). Also accumulates the primary's PLA path
// length (t = integral chi dz, the geometry probe) and, in reflective single-cell
// mode, specularly folds the primary at the transverse container walls.
#include "SteppingAction.hh"

#include <cmath>

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "G4Material.hh"
#include "G4Step.hh"
#include "G4Track.hh"

SteppingAction::SteppingAction(EventAction* evt,
                               const DetectorConstruction* det)
    : fEventAction(evt), fDetector(det) {}

void SteppingAction::UserSteppingAction(const G4Step* step) {
  auto track = step->GetTrack();
  if (track->GetParentID() != 0) return;  // primary only

  // Kink-scoring planes. For a target tilted by theta (M2), the beam is inside the target
  // over a z-span of depth/cos(theta), which EXCEEDS +/-0.5*depth. Push the scoring planes
  // just outside the tilted z-extent so the recorded entry/exit directions bracket the FULL
  // traversal (they sit in vacuum, so the direction is unchanged there). tilt=0 -> +/-0.5*depth,
  // exactly the normal-incidence scoring (campaign unchanged).
  const G4double tilt = fDetector->GetTilt();
  const G4double halfDepth = (std::abs(tilt) < 1e-9)
      ? 0.5 * fDetector->GetTargetDepth()
      : 0.5 * fDetector->GetTargetDepth() / std::cos(tilt) + 2.0;  // +2 mm margin (mm units)
  const auto pre = step->GetPreStepPoint();
  const auto post = step->GetPostStepPoint();
  const G4double z0 = pre->GetPosition().z();
  const G4double z1 = post->GetPosition().z();

  // Entry plane z = -halfDepth: record the incident direction.
  if (z0 < -halfDepth && z1 >= -halfDepth)
    fEventAction->SetEntryDir(post->GetMomentumDirection());
  // Exit plane z = +halfDepth: record the outgoing direction.
  if (z0 < halfDepth && z1 >= halfDepth)
    fEventAction->SetExitDir(post->GetMomentumDirection());

  // PLA path length: add the step length when the step is INSIDE target material.
  auto preMat = pre->GetMaterial();
  if (preMat && preMat->GetName() == fDetector->GetMaterialName())
    fEventAction->AddPathLen(step->GetStepLength());

  // M1 telescope: record the primary's (x,y) where this step crosses each scoring
  // plane z (linear interpolation between pre/post points). Gated OFF by default.
  if (fDetector->GetTelescope()) {
    const G4ThreeVector p0 = pre->GetPosition();
    const G4ThreeVector p1 = post->GetPosition();
    for (G4int i = 0; i < fDetector->GetNPlanes(); ++i) {
      const G4double zp = fDetector->GetPlaneZ(i);
      if ((z0 < zp && z1 >= zp) || (z0 > zp && z1 <= zp)) {
        const G4double dz = z1 - z0;
        const G4double frac = (std::abs(dz) > 1e-12) ? (zp - z0) / dz : 0.0;
        fEventAction->SetPlaneHit(i, p0.x() + frac * (p1.x() - p0.x()),
                                  p0.y() + frac * (p1.y() - p0.y()));
      }
    }
  }

  // Reflective single transverse cell: specularly fold the primary at the cell
  // walls x,y = +/- cell/2 so a one-cell block mimics an extended lattice.
  if (fDetector->GetReflective()) {
    const G4double half = 0.5 * fDetector->GetCellSize();
    const G4double eps = 1e-7;
    G4ThreeVector dir = track->GetMomentumDirection();
    const G4ThreeVector p = post->GetPosition();
    G4bool fold = false;
    if ((p.x() >= half - eps && dir.x() > 0) ||
        (p.x() <= -half + eps && dir.x() < 0)) { dir.setX(-dir.x()); fold = true; }
    if ((p.y() >= half - eps && dir.y() > 0) ||
        (p.y() <= -half + eps && dir.y() < 0)) { dir.setY(-dir.y()); fold = true; }
    if (fold) track->SetMomentumDirection(dir.unit());
  }
}
