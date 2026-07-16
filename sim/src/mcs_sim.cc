// mcs_sim --- main() for the MCS kink-angle transport engine (S3).
// Usage:
//   mcs_sim run.mac [nThreads]    batch (default threads = ncores-2)
//   mcs_sim                       interactive UI + vis (vis.mac)
#include <cstdlib>

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4String.hh"
#include "G4Types.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "PhysicsList.hh"
#include "Randomize.hh"

int main(int argc, char** argv) {
  G4UIExecutive* ui = nullptr;
  if (argc == 1) ui = new G4UIExecutive(argc, argv);

  auto runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
  if (argc > 2) runManager->SetNumberOfThreads(std::atoi(argv[2]));

  auto detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new PhysicsList);
  runManager->SetUserInitialization(new ActionInitialization(detector));

  auto visManager = new G4VisExecutive;
  visManager->Initialize();
  auto uiManager = G4UImanager::GetUIpointer();

  if (ui) {
    uiManager->ApplyCommand("/control/execute macros/vis.mac");
    ui->SessionStart();
    delete ui;
  } else {
    uiManager->ApplyCommand(G4String("/control/execute ") + argv[1]);
  }

  delete visManager;
  delete runManager;
  return 0;
}
