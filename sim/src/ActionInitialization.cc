// ActionInitialization implementation (S3).
#include "ActionInitialization.hh"

#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization(const DetectorConstruction* det)
    : fDetector(det) {}

void ActionInitialization::BuildForMaster() const {
  SetUserAction(new RunAction);  // master writes the merged file + sidecar
}

void ActionInitialization::Build() const {
  SetUserAction(new PrimaryGeneratorAction);
  SetUserAction(new RunAction);
  auto eventAction = new EventAction(fDetector);
  SetUserAction(eventAction);
  SetUserAction(new SteppingAction(eventAction, fDetector));
}
