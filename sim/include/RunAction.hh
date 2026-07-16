// RunAction --- owns the G4AnalysisManager ntuple of per-primary kink angles
// (S3). ROOT backend, MT-merged into a single file. On the master thread at end
// of run it also writes a JSON sidecar of run metadata (X0 from GetRadlen,
// thickness, beam kinetic energy, cut, max step, N events) consumed by the
// Highland-validation script.
#ifndef MCS_RUNACTION_HH
#define MCS_RUNACTION_HH

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class RunAction : public G4UserRunAction {
 public:
  RunAction();
  ~RunAction() override = default;

  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

 private:
  void WriteMetadataSidecar(const G4Run* run) const;
};

#endif  // MCS_RUNACTION_HH
