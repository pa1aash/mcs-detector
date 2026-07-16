// EventAction --- per-primary projected kink angle (S3) + optional M1 telescope.
// SteppingAction sets the primary's momentum direction at target entry and at
// target exit; at end of event the two projected kink angles (x-z and y-z
// planes) are computed and filled into the "kinks" ntuple. When the telescope is
// enabled (M1), SteppingAction also records the primary's (x,y) crossing at each
// of the 6 scoring planes, filled into the "telescope" ntuple.
#ifndef MCS_EVENTACTION_HH
#define MCS_EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"  // for kNPlanes + telescope query

class EventAction : public G4UserEventAction {
 public:
  explicit EventAction(const DetectorConstruction* det = nullptr) : fDet(det) {}
  ~EventAction() override = default;

  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

  void SetEntryDir(const G4ThreeVector& d) { fEntryDir = d; fHasEntry = true; }
  void SetExitDir(const G4ThreeVector& d) { fExitDir = d; fHasExit = true; }
  // PLA path length per primary (geometry probe: t = integral chi dz).
  void AddPathLen(G4double dl) { fPathLenPLA += dl; }
  // M1: (x,y) of the primary at telescope plane i.
  void SetPlaneHit(G4int i, G4double x, G4double y) {
    if (i >= 0 && i < DetectorConstruction::kNPlanes) {
      fPlaneX[i] = x; fPlaneY[i] = y; fPlaneHit[i] = true;
    }
  }

 private:
  const DetectorConstruction* fDet = nullptr;
  G4ThreeVector fEntryDir;
  G4ThreeVector fExitDir;
  G4bool fHasEntry = false;
  G4bool fHasExit = false;
  G4double fPathLenPLA = 0.;  // mm of target material traversed by the primary

  // M1 telescope plane hits.
  G4double fPlaneX[DetectorConstruction::kNPlanes] = {0.};
  G4double fPlaneY[DetectorConstruction::kNPlanes] = {0.};
  G4bool fPlaneHit[DetectorConstruction::kNPlanes] = {false};
};

#endif  // MCS_EVENTACTION_HH
