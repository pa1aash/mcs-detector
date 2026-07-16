// PhysicsList implementation (S3) --- LOCKED configuration by default.
//
// S7 / E5 (MSC-model systematic): the EM constructor and the WentzelVI
// MSC<->single-Coulomb balance are now selectable via environment variables so
// the locked run and the alternative-MSC run share one binary. With NO env vars
// set the configuration is BIT-IDENTICAL to the S3 lock (option4, default
// WentzelVI parameters) --- the lock is the default; the knobs are opt-in.
//
//   MCS_EM_PHYSICS = option4 (default, LOCKED) | option3 (UrbanMsc-based)
//   MCS_MSC_THETALIMIT_DEG = <deg>   -> G4EmParameters::SetMscThetaLimit
//                                       (WentzelVI MSC / single-Coulomb split;
//                                        default in G4 is 3.1416 rad = pi)
//   MCS_MSC_MUHAD_RANGEFAC = <f>     -> SetMscMuHadRangeFactor (hadron MSC step;
//                                       default 0.2)
#include "PhysicsList.hh"

#include <cstdlib>
#include <string>

#include "G4EmParameters.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4SystemOfUnits.hh"

namespace {
std::string env_or(const char* k, const std::string& dflt) {
  const char* v = std::getenv(k);
  return (v && *v) ? std::string(v) : dflt;
}
}  // namespace

PhysicsList::PhysicsList(G4int verbose) : FTFP_BERT(verbose) {
  const std::string em = env_or("MCS_EM_PHYSICS", "option4");

  if (em == "option3") {
    // ALTERNATIVE MSC (E5 systematic): option3's UrbanMsc-based EM, instead of
    // the locked option4 WentzelVI + single-Coulomb split. Used ONLY for the
    // MSC model-systematic study.
    ReplacePhysics(new G4EmStandardPhysics_option3(verbose));
  } else {
    // LOCKED default: option4 = WentzelVI multiple scattering + single Coulomb
    // scattering (G4eCoulombScatteringModel) for protons/hadrons.
    ReplacePhysics(new G4EmStandardPhysics_option4(verbose));
  }

  // E5 WentzelVI-balance knobs (no-op when unset -> identical to the lock).
  auto* emp = G4EmParameters::Instance();
  const std::string th = env_or("MCS_MSC_THETALIMIT_DEG", "");
  if (!th.empty()) emp->SetMscThetaLimit(std::stod(th) * deg);
  const std::string rf = env_or("MCS_MSC_MUHAD_RANGEFAC", "");
  if (!rf.empty()) emp->SetMscMuHadRangeFactor(std::stod(rf));

  // Make the in-target G4UserLimits max step actually bite.
  RegisterPhysics(new G4StepLimiterPhysics());

  // Default production range cut (overridable from the macro via /run/setCut).
  SetDefaultCutValue(0.05 * mm);
}

void PhysicsList::SetCuts() {
  // Honour SetDefaultCutValue / any /run/setCut, then let FTFP_BERT apply.
  FTFP_BERT::SetCuts();
}
