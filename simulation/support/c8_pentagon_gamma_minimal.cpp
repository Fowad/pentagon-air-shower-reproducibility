/*
 * Paper-specific minimal CORSIKA 8 gamma-shower application for the
 * five-detector Pentagon response study.
 *
 * Scope is deliberately narrow: gamma primaries, SIBYLL 2.3d high-energy
 * hadronic interactions, UrQMD low-energy hadronic interactions, SIBYLL
 * decays, PROPOSAL electromagnetic transport, SOPHIA resonance-region
 * photonuclear interactions, particle cuts, EM thinning support, and an
 * absorbing observation plane with Parquet ground-particle output.
 *
 * It intentionally omits alternative HE models, Pythia8, TAUOLA, CONEX,
 * radio emission, StackInspector, interaction histograms, longitudinal and
 * production profiles, and neutrino-primary machinery.  Those are not needed
 * for the manuscript's baseline gamma-response question.
 */

#include <corsika/framework/core/Cascade.hpp>
#include <corsika/framework/core/EnergyMomentumOperations.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/geometry/Plane.hpp>
#include <corsika/framework/process/InteractionCounter.hpp>
#include <corsika/framework/process/ProcessSequence.hpp>
#include <corsika/framework/process/SwitchProcessSequence.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/utility/CorsikaFenv.hpp>

#include <corsika/media/CORSIKA7Atmospheres.hpp>
#include <corsika/media/Environment.hpp>
#include <corsika/media/LayeredSphericalAtmosphereBuilder.hpp>
#include <corsika/media/ShowerAxis.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/interfaces/IMagneticFieldModel.hpp>
#include <corsika/media/magnetic/GeomagneticModel.hpp>
#include <corsika/media/magnetic/UniformMagneticField.hpp>
#include <corsika/media/medium/MediumPropertyModel.hpp>
#include <corsika/media/refractivity/GladstoneDaleRefractiveIndex.hpp>

#include <corsika/modules/BetheBlochPDG.hpp>
#include <corsika/modules/ObservationPlane.hpp>
#include <corsika/modules/PROPOSAL.hpp>
#include <corsika/modules/ParticleCut.hpp>
#include <corsika/modules/Sibyll.hpp>
#include <corsika/modules/Sophia.hpp>
#include <corsika/modules/UrQMD.hpp>
#include <corsika/modules/sibyll/Decay.hpp>
#include <corsika/modules/thinning/EMThinning.hpp>
#include <corsika/modules/writers/EnergyLossWriter.hpp>
#include <corsika/modules/writers/PrimaryWriter.hpp>
#include <corsika/modules/writers/SubWriter.hpp>
#include <corsika/output/OutputManager.hpp>

#include <corsika/setup/SetupC7trackedParticles.hpp>
#include <corsika/setup/SetupStack.hpp>
#include <corsika/setup/SetupTrajectory.hpp>

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <random>
#include <sstream>
#include <string>
#include <tuple>

using namespace corsika;
using namespace std;

using EnvironmentInterface = media::IRefractiveIndexModel<
    media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>>;
using EnvType = media::Environment<EnvironmentInterface>;
using StackType = setup::Stack<EnvType>;
using TrackingType = setup::Tracking;
using Particle = StackType::particle_type;

template <typename T>
using MyExtraEnv = media::GladstoneDaleRefractiveIndex<
    media::MediumPropertyModel<media::UniformMagneticField<T>>>;

long registerRandomStreams(long seed) {
  // Only streams used by this minimal paper application are registered.
  RNGManager<>::getInstance().registerRandomStream("cascade");
  RNGManager<>::getInstance().registerRandomStream("sibyll");
  RNGManager<>::getInstance().registerRandomStream("sophia");
  RNGManager<>::getInstance().registerRandomStream("urqmd");
  RNGManager<>::getInstance().registerRandomStream("proposal");
  RNGManager<>::getInstance().registerRandomStream("thinning");
  if (seed == 0) {
    std::random_device rd;
    seed = rd();
    CORSIKA_LOG_INFO("random seed (auto) {}", seed);
  } else {
    CORSIKA_LOG_INFO("random seed {}", seed);
  }
  RNGManager<>::getInstance().setSeed(seed);
  return seed;
}

int main(int argc, char** argv) {
  CLI::App app{"Minimal gamma-shower CORSIKA 8 application for the Pentagon paper."};

  double energy_gev = 0.0;
  double zenith_deg = 0.0;
  double azimuth_deg = 0.0;
  double observation_level_m = 1200.0;
  double injection_height_m = 112.75e3;
  double geomagnetic_year = 2017.0;
  double site_latitude_deg = 35.704;
  double site_longitude_deg = 51.351;
  double emcut_gev = 0.5e-3;
  double hadcut_gev = 0.3;
  double mucut_gev = 0.3;
  double transition_gev = std::pow(10.0, 1.9);
  double emthin_fraction = 1.e-6;
  double max_weight_cli = 1.0;
  double max_deflection_angle = 0.2;
  int nevent = 1;
  long seed_cli = 0;
  bool multithin = false;
  bool compress_output = false;
  std::string filename;
  std::string verbosity = "info";

  app.add_option("-E,--energy", energy_gev, "Gamma primary energy in GeV")
      ->required()->check(CLI::PositiveNumber)->group("Primary");
  app.add_option("-z,--zenith", zenith_deg, "Primary zenith angle in degrees")
      ->default_val(0.0)->check(CLI::Range(0.0, 90.0))->group("Primary");
  app.add_option("-a,--azimuth", azimuth_deg, "Primary azimuth angle in degrees")
      ->default_val(0.0)->check(CLI::Range(0.0, 360.0))->group("Primary");

  app.add_option("--observation-level", observation_level_m,
                 "Observation altitude above Earth radius in m")
      ->default_val(1200.0)->check(CLI::Range(-1.e3, 1.e5))->group("Site");
  app.add_option("--injection-height", injection_height_m,
                 "Injection altitude above Earth radius in m")
      ->default_val(112.75e3)->check(CLI::Range(-1.e3, 1.e6))->group("Site");
  app.add_option("--geomagnetic-year", geomagnetic_year,
                 "Decimal year for IGRF-13 geomagnetic field")
      ->default_val(2017.0)->check(CLI::Range(1900.0, 2025.0))->group("Site");
  app.add_option("--site-latitude", site_latitude_deg,
                 "Geodetic latitude, north positive, degrees")
      ->default_val(35.704)->check(CLI::Range(-90.0, 90.0))->group("Site");
  app.add_option("--site-longitude", site_longitude_deg,
                 "Geodetic longitude, east positive, degrees")
      ->default_val(51.351)->check(CLI::Range(-180.0, 180.0))->group("Site");

  app.add_option("--emcut", emcut_gev,
                 "Minimum kinetic energy of photons/electrons/positrons in GeV")
      ->default_val(0.5e-3)->check(CLI::Range(0.000001, 1.e13))->group("Physics");
  app.add_option("--hadcut", hadcut_gev, "Minimum kinetic energy of hadrons in GeV")
      ->default_val(0.3)->check(CLI::Range(0.02, 1.e13))->group("Physics");
  app.add_option("--mucut", mucut_gev, "Minimum kinetic energy of muons in GeV")
      ->default_val(0.3)->check(CLI::Range(0.000001, 1.e13))->group("Physics");
  app.add_option("-T,--hadronModelTransitionEnergy", transition_gev,
                 "SIBYLL/UrQMD transition energy in GeV")
      ->default_val(std::pow(10.0, 1.9))->check(CLI::NonNegativeNumber)->group("Physics");
  app.add_option("--max-deflection-angle", max_deflection_angle,
                 "Maximum tracking deflection angle in radians")
      ->default_val(0.2)->check(CLI::Range(1.e-8, 1.0))->group("Physics");

  app.add_option("--emthin", emthin_fraction,
                 "Fraction of primary energy at which EM thinning starts")
      ->default_val(1.e-6)->check(CLI::Range(0.0, 1.0))->group("Thinning");
  app.add_option("--max-weight", max_weight_cli, "Maximum EM thinning weight")
      ->default_val(1.0)->check(CLI::NonNegativeNumber)->group("Thinning");
  app.add_flag("--multithin", multithin, "Keep thinned particles with zero weight")
      ->group("Thinning");

  app.add_option("-N,--nevent", nevent, "Number of showers")
      ->default_val(1)->check(CLI::PositiveNumber)->group("Output");
  app.add_option("-s,--seed", seed_cli, "Random-number seed")
      ->default_val(0)->check(CLI::NonNegativeNumber)->group("Output");
  app.add_option("-f,--filename", filename, "Output CORSIKA library directory")
      ->required()->check(CLI::NonexistentPath)->group("Output");
  app.add_flag("--compress", compress_output, "Compress output directory")
      ->group("Output");
  app.add_option("-v,--verbosity", verbosity, "warn, info, debug, or trace")
      ->default_val("info")->check(CLI::IsMember({"warn","info","debug","trace"}));

  CLI11_PARSE(app, argc, argv);

  if (verbosity == "warn") logging::set_level(logging::level::warn);
  else if (verbosity == "info") logging::set_level(logging::level::info);
  else if (verbosity == "debug") logging::set_level(logging::level::debug);
  else if (verbosity == "trace") {
#ifndef _C8_DEBUG_
    CORSIKA_LOG_ERROR("trace log level requires a Debug build");
    return EXIT_FAILURE;
#else
    logging::set_level(logging::level::trace);
#endif
  }

  auto const seed = registerRandomStreams(seed_cli);

  // Site/environment: fixed paper baseline = USStdBK atmosphere + IGRF-13.
  EnvType env;
  CoordinateSystemPtr const& rootCS = env.getCoordinateSystem();
  Point const center{rootCS, 0_m, 0_m, 0_m};
  Point const surface{rootCS, 0_m, 0_m, constants::EarthRadius::Mean};
  media::GeomagneticModel igrf(center, corsika_data("GeoMag/IGRF13.COF"));
  auto const magneticField = igrf.getField(geomagnetic_year, observation_level_m * 1_m,
                                           site_latitude_deg, site_longitude_deg);
  CORSIKA_LOG_INFO(
      "PAPER-MINIMAL configuration: gamma; SIBYLL-2.3d HE; UrQMD LE; SIBYLL decays; "
      "PROPOSAL EM; SOPHIA resonance photonuclear; USStdBK atmosphere");
  CORSIKA_LOG_INFO(
      "IGRF-13 field: year={} latitude={} longitude={} altitude={} m; "
      "(north,west,up)=({}, {}, {}) uT",
      geomagnetic_year, site_latitude_deg, site_longitude_deg, observation_level_m,
      magneticField.getX(rootCS) / 1_uT, magneticField.getY(rootCS) / 1_uT,
      magneticField.getZ(rootCS) / 1_uT);
  media::create_5layer_atmosphere<EnvironmentInterface, MyExtraEnv>(
      env, media::AtmosphereId::USStdBK, center, 1.000327, surface,
      media::Medium::AirDry1Atm, magneticField);

  Code const beamCode = convert_from_PDG(PDGCode(22));
  HEPEnergyType const primaryTotalEnergy = energy_gev * 1_GeV;
  auto const thetaRad = zenith_deg / 180.0 * M_PI;
  auto const phiRad = azimuth_deg / 180.0 * M_PI;
  auto const [nx, ny, nz] = std::make_tuple(std::sin(thetaRad) * std::cos(phiRad),
                                            std::sin(thetaRad) * std::sin(phiRad),
                                            -std::cos(thetaRad));
  auto const propDir = DirectionVector(rootCS, {nx, ny, nz});

  auto const observationHeight = observation_level_m * 1_m + constants::EarthRadius::Mean;
  auto const injectionHeight = injection_height_m * 1_m + constants::EarthRadius::Mean;
  // Keep sqrt unqualified so ADL selects the PhysicalUnits overload for a length^2 quantity.
  auto const t = -observationHeight * std::cos(thetaRad) +
                 sqrt(-static_pow<2>(std::sin(thetaRad) * observationHeight) +
                      static_pow<2>(injectionHeight));
  Point const showerCore{rootCS, 0_m, 0_m, observationHeight};
  Point const injectionPos =
      showerCore + DirectionVector{rootCS,
                                   {-std::sin(thetaRad) * std::cos(phiRad),
                                    -std::sin(thetaRad) * std::sin(phiRad),
                                    std::cos(thetaRad)}} * t;
  media::ShowerAxis const showerAxis{injectionPos, (showerCore - injectionPos) * 1.2, env};
  auto const dX = 10_g / square(1_cm);

  std::stringstream args;
  for (int i = 0; i < argc; ++i) args << argv[i] << " ";
  OutputManager output(filename, seed, args.str(), compress_output);

  // Energy-loss writer is retained because ParticleCut and continuous-loss processes
  // use it for energy accounting.  No longitudinal/production/interaction outputs.
  EnergyLossWriter dEdX{showerAxis, dX};
  output.add("energyloss", dEdX);

  auto const trackedParticles = corsika::setup::C7trackedParticles;
  auto const all_elements = corsika::media::get_all_elements_in_universe(env);
  corsika::sibyll::Interaction heIntModel(all_elements, trackedParticles);
  corsika::urqmd::UrQMD leIntModel{};

  // SIBYLL provides the decay process in this minimal baseline, avoiding Pythia8/TAUOLA.
  // In the locked CORSIKA SIBYLL Decay class, handleAllDecays_ defaults to true.
  // Do not call the setHandleAllDecay method: at this commit that redundant non-template
  // setter is not exported by the minimal SIBYLL link target on macOS.
  corsika::sibyll::Decay decaySibyll;

  // Photonuclear resonance region used by PROPOSAL.
  corsika::sophia::InteractionModel sophia;

  HEPEnergyType const emcut = emcut_gev * 1_GeV;
  HEPEnergyType const hadcut = hadcut_gev * 1_GeV;
  HEPEnergyType const mucut = mucut_gev * 1_GeV;
  HEPEnergyType const taucut = mucut; // no tau-specific study; same tracking floor as muons
  ParticleCut<SubWriter<decltype(dEdX)>> cut(emcut, emcut, hadcut, mucut, taucut,
                                             true /* cut neutrinos */, dEdX);

  auto const prod_threshold = std::min({emcut, hadcut, mucut, taucut});
  set_energy_production_threshold(Code::Electron, prod_threshold);
  set_energy_production_threshold(Code::Positron, prod_threshold);
  set_energy_production_threshold(Code::Photon, prod_threshold);
  set_energy_production_threshold(Code::MuMinus, prod_threshold);
  set_energy_production_threshold(Code::MuPlus, prod_threshold);
  set_energy_production_threshold(Code::TauMinus, prod_threshold);
  set_energy_production_threshold(Code::TauPlus, prod_threshold);

  HEPEnergyType const heThreshold = transition_gev * 1_GeV;
  corsika::proposal::Interaction emCascade(
      env, sophia, heIntModel.getHadronInteractionModel(), heThreshold);

  corsika::proposal::ContinuousProcess<SubWriter<decltype(dEdX)>> emContinuousProposal(
      env, dEdX);
  BetheBlochPDG<SubWriter<decltype(dEdX)>> emContinuousBethe{dEdX};
  struct EMHadronSwitch {
    bool operator()(const Particle& p) const { return is_hadron(p.getPID()); }
  };
  auto emContinuous = make_select(EMHadronSwitch(), emContinuousBethe,
                                  emContinuousProposal);

  struct EnergySwitch {
    HEPEnergyType cutE;
    bool operator()(const Particle& p) const { return p.getKineticEnergy() < cutE; }
  };
  // Interaction models are not themselves ProcessSequence processes.  Wrap them in the
  // same InteractionCounter process adapter used by the locked upstream c8_air_shower.
  // The counters are not registered as manuscript output; they are used here only as the
  // process wrapper required by the CORSIKA framework API.
  InteractionCounter leIntProcess{leIntModel};
  InteractionCounter heIntProcess{heIntModel};
  auto hadronSequence =
      make_select(EnergySwitch{heThreshold}, leIntProcess, heIntProcess);

  Plane const obsPlane(showerCore, DirectionVector(rootCS, {0., 0., 1.}));
  ObservationPlane<TrackingType, ParticleWriterParquet> observationLevel{
      obsPlane, DirectionVector(rootCS, {1., 0., 0.}), true, 1e-6 * 1_m, false};
  output.add("particles", observationLevel);
  PrimaryWriter<TrackingType, ParticleWriterParquet> primaryWriter(observationLevel);
  output.add("primary", primaryWriter);

  output.startOfLibrary();
  for (int i_shower = 1; i_shower <= nevent; ++i_shower) {
    CORSIKA_LOG_INFO("Minimal gamma shower {} / {}", i_shower, nevent);
    double const maxWeight = max_weight_cli > 0.0
                                 ? max_weight_cli
                                 : 0.5 * emthin_fraction * primaryTotalEnergy / 1_GeV;
    EMThinning thinning{emthin_fraction * primaryTotalEnergy, maxWeight, !multithin};

    auto sequence = make_sequence(hadronSequence, decaySibyll, emCascade,
                                  emContinuous, observationLevel, thinning, cut);
    TrackingType tracking(max_deflection_angle);
    StackType stack;
    Cascade EAS(env, tracking, sequence, output, stack);
    stack.clear();

    auto const primaryProperties =
        std::make_tuple(beamCode, primaryTotalEnergy, propDir.normalized(), injectionPos, 0_ns);
    stack.addParticle(primaryProperties);
    primaryWriter.recordPrimary(primaryProperties);

    CORSIKA_LOG_INFO("Primary: gamma, E={} GeV, zenith={} deg, azimuth={} deg",
                     primaryTotalEnergy / 1_GeV, zenith_deg, azimuth_deg);
    EAS.run();

    HEPEnergyType const Efinal = dEdX.getEnergyLost() + observationLevel.getEnergyGround();
    CORSIKA_LOG_INFO("Energy accounting (GeV): total={} dEdX={} ground={}",
                     Efinal / 1_GeV, dEdX.getEnergyLost() / 1_GeV,
                     observationLevel.getEnergyGround() / 1_GeV);
  }

  output.endOfLibrary();
  return EXIT_SUCCESS;
}
