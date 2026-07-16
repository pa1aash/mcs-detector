// EventAction implementation (S3 + M1 telescope).
#include "EventAction.hh"

#include <cmath>

#include "G4AnalysisManager.hh"
#include "G4Event.hh"

void EventAction::BeginOfEventAction(const G4Event*) {
  fHasEntry = false;
  fHasExit = false;
  fPathLenPLA = 0.;
  for (G4int i = 0; i < DetectorConstruction::kNPlanes; ++i) fPlaneHit[i] = false;
}

void EventAction::EndOfEventAction(const G4Event*) {
  if (!(fHasEntry && fHasExit)) return;  // primary must traverse the target

  // Projected kink angle = exit minus entry projected angle, per plane.
  const G4double thetax = std::atan2(fExitDir.x(), fExitDir.z()) -
                          std::atan2(fEntryDir.x(), fEntryDir.z());
  const G4double thetay = std::atan2(fExitDir.y(), fExitDir.z()) -
                          std::atan2(fEntryDir.y(), fEntryDir.z());

  auto am = G4AnalysisManager::Instance();
  am->FillNtupleDColumn(0, thetax);      // "kinks" ntuple (id 0) -- unchanged
  am->FillNtupleDColumn(1, thetay);
  am->FillNtupleDColumn(2, fPathLenPLA);
  am->AddNtupleRow();

  // M1 telescope: fill the "telescope" ntuple (id 1) when enabled and all 6
  // planes were crossed. Columns: x0..x5 (0..5), y0..y5 (6..11), tpla (12),
  // thetax (13), thetay (14). Positions in mm (G4 internal units -> /mm below).
  if (fDet && fDet->GetTelescope()) {
    G4bool all = true;
    for (G4int i = 0; i < DetectorConstruction::kNPlanes; ++i)
      if (!fPlaneHit[i]) { all = false; break; }
    if (all) {
      const G4int nP = DetectorConstruction::kNPlanes;
      // Positions are in Geant4 base length units (mm), stored directly.
      for (G4int i = 0; i < nP; ++i) {
        am->FillNtupleDColumn(1, i, fPlaneX[i]);
        am->FillNtupleDColumn(1, nP + i, fPlaneY[i]);
      }
      am->FillNtupleDColumn(1, 2 * nP + 0, fPathLenPLA);
      am->FillNtupleDColumn(1, 2 * nP + 1, thetax);
      am->FillNtupleDColumn(1, 2 * nP + 2, thetay);
      am->AddNtupleRow(1);
    }
  }
}
