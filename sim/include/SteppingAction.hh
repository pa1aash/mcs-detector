// SteppingAction --- record the primary's momentum direction at target entry
// and exit boundaries (S3). Primary track only (parentID == 0).
#ifndef MCS_STEPPINGACTION_HH
#define MCS_STEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class EventAction;
class DetectorConstruction;
class G4VPhysicalVolume;

class SteppingAction : public G4UserSteppingAction {
 public:
  SteppingAction(EventAction* evt, const DetectorConstruction* det);
  ~SteppingAction() override = default;

  void UserSteppingAction(const G4Step* step) override;

 private:
  EventAction* fEventAction = nullptr;
  const DetectorConstruction* fDetector = nullptr;
};

#endif  // MCS_STEPPINGACTION_HH
