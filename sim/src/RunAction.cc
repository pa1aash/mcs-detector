// RunAction implementation (S3).
#include "RunAction.hh"

#include <filesystem>
#include <fstream>
#include <string>
#include <system_error>

#include "DetectorConstruction.hh"
#include "G4AnalysisManager.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "Randomize.hh"

RunAction::RunAction() {
  // ROOT ntuple, MT-merged into a single file. One row per traversing primary.
  auto am = G4AnalysisManager::Instance();
  am->SetDefaultFileType("root");
  am->SetNtupleMerging(true);
  am->SetVerboseLevel(0);
  am->SetFileName("data/runs/run");  // default; macro overrides via /analysis/setFileName
  am->CreateNtuple("kinks", "projected MCS kink angles [rad]");
  am->CreateNtupleDColumn("thetax");
  am->CreateNtupleDColumn("thetay");
  am->CreateNtupleDColumn("tpla");   // PLA path length [mm] (geometry probe)
  am->FinishNtuple();

  // M1 telescope ntuple (id 1): filled only when /mcs/det/telescope 1 (else empty).
  am->CreateNtuple("telescope", "6-plane telescope hits [mm] + truth kink [rad]");
  for (G4int i = 0; i < DetectorConstruction::kNPlanes; ++i)
    am->CreateNtupleDColumn("x" + std::to_string(i));
  for (G4int i = 0; i < DetectorConstruction::kNPlanes; ++i)
    am->CreateNtupleDColumn("y" + std::to_string(i));
  am->CreateNtupleDColumn("tpla");
  am->CreateNtupleDColumn("thetax");
  am->CreateNtupleDColumn("thetay");
  am->FinishNtuple();
}

void RunAction::BeginOfRunAction(const G4Run*) {
  auto am = G4AnalysisManager::Instance();
  std::error_code ec;
  const auto parent = std::filesystem::path(am->GetFileName()).parent_path();
  if (!parent.empty()) std::filesystem::create_directories(parent, ec);
  am->OpenFile();
}

void RunAction::EndOfRunAction(const G4Run* run) {
  auto am = G4AnalysisManager::Instance();
  am->Write();
  am->CloseFile();
  if (IsMaster()) WriteMetadataSidecar(run);
}

void RunAction::WriteMetadataSidecar(const G4Run* run) const {
  auto am = G4AnalysisManager::Instance();
  auto det = dynamic_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  if (!det) return;

  // Global production range cut (set by /run/setCut); all particles share it.
  G4double cut_mm = 0.;
  auto cuts = G4ProductionCutsTable::GetProductionCutsTable()
                  ->GetDefaultProductionCuts();
  if (cuts) cut_mm = cuts->GetProductionCut("gamma") / mm;

  const G4double thickness_mm = det->GetThickness() / mm;
  const G4double x0_mm = det->GetTargetX0() / mm;
  const G4int nevt = run->GetNumberOfEvent();

  // E5 provenance: record the EM-physics / MSC configuration (env-driven).
  auto envstr = [](const char* k, const char* dflt) {
    const char* v = std::getenv(k);
    return std::string((v && *v) ? v : dflt);
  };
  const std::string em_phys = envstr("MCS_EM_PHYSICS", "option4");
  const std::string msc_theta = envstr("MCS_MSC_THETALIMIT_DEG", "");
  const std::string msc_rf = envstr("MCS_MSC_MUHAD_RANGEFAC", "");

  // Reproducibility provenance: the actual master RNG seed used (set via the
  // macro's /random/setSeeds; default Geant4 seeding otherwise) and the code
  // version. git_commit is passed by the campaign wrapper (MCS_GIT_COMMIT).
  const long seed = G4Random::getTheSeed();
  const std::string git_commit = envstr("MCS_GIT_COMMIT", "");

  const std::string path = std::string(am->GetFileName()) + ".meta.json";
  std::ofstream f(path);
  f << "{\n"
    << "  \"em_physics\": \"" << em_phys << "\",\n"
    << "  \"msc_thetalimit_deg\": \"" << msc_theta << "\",\n"
    << "  \"msc_muhad_rangefac\": \"" << msc_rf << "\",\n"
    << "  \"material\": \"" << det->GetMaterialName() << "\",\n"
    << "  \"thickness_mm\": " << thickness_mm << ",\n"
    << "  \"density_g_cm3\": " << det->GetTargetDensity() << ",\n"
    << "  \"X0_mm\": " << x0_mm << ",\n"
    << "  \"t_over_X0\": " << (x0_mm > 0 ? thickness_mm / x0_mm : 0.) << ",\n"
    << "  \"maxStep_mm\": " << det->GetMaxStep() / mm << ",\n"
    << "  \"cut_mm\": " << cut_mm << ",\n"
    << "  \"geom\": \"" << det->GetGeomType() << "\",\n"
    << "  \"topology\": \"" << det->GetTopology() << "\",\n"
    << "  \"infill\": " << det->GetInfill() << ",\n"
    << "  \"cell_mm\": " << det->GetCellSize() / mm << ",\n"
    << "  \"tilt_deg\": " << det->GetTiltDeg() << ",\n"
    << "  \"depth_mm\": " << det->GetTargetDepth() / mm << ",\n"
    << "  \"ntile\": " << det->GetNTile() << ",\n"
    << "  \"reflective\": " << (det->GetReflective() ? 1 : 0) << ",\n"
    << "  \"seed\": " << seed << ",\n"
    << "  \"git_commit\": \"" << git_commit << "\",\n"
    << "  \"n_events\": " << nevt << "\n"
    << "}\n";
  f.close();
  G4cout << "[RunAction] wrote metadata sidecar: " << path << G4endl;
}
