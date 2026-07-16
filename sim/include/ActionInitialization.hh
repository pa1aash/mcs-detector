// ActionInitialization --- wires user actions for master + worker threads (S3).
#ifndef MCS_ACTIONINITIALIZATION_HH
#define MCS_ACTIONINITIALIZATION_HH

#include "G4VUserActionInitialization.hh"

class DetectorConstruction;

class ActionInitialization : public G4VUserActionInitialization {
 public:
  explicit ActionInitialization(const DetectorConstruction* det);
  ~ActionInitialization() override = default;

  void Build() const override;
  void BuildForMaster() const override;

 private:
  const DetectorConstruction* fDetector = nullptr;
};

#endif  // MCS_ACTIONINITIALIZATION_HH
