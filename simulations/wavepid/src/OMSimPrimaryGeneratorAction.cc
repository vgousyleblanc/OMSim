/**
 * @file OMSimPrimaryGeneratorAction.cc
 * @brief Implementation of primary generator action for WavePID.
 *
 * Calculates muon starting position based on impact parameter to ensure
 * the Cherenkov light cone hits the DOM at the origin.
 *
 * Geometry:
 *   - DOM at origin (0, 0, 0)
 *   - Muon travels in +x direction
 *   - Muon track at perpendicular distance r (impact parameter) in z
 *   - Starting x position: -1.7*r (ensures Cherenkov cone reaches DOM)
 *
 * The factor 1.7 accounts for the Cherenkov angle in ice (~41 deg):
 *   tan(41°) ≈ 0.87, so light needs ~1.15*r to reach DOM
 *   Using 1.7 provides margin for the full cone to illuminate the DOM.
 */
#include "OMSimPrimaryGeneratorAction.hh"
#include <cmath>
#include <Randomize.hh>
#include <G4OpticalPhoton.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>
#include "OMSimCommandArgsTable.hh"
#include "OMSimLogger.hh"

#include <G4ParticleGun.hh>
#include <G4GeneralParticleSource.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4SystemOfUnits.hh>
#include <G4RandomTools.hh>

thread_local std::unique_ptr<G4ParticleGun> OMSimPrimaryGeneratorAction::m_particleGun;
thread_local std::unique_ptr<G4GeneralParticleSource> OMSimPrimaryGeneratorAction::m_particleSource;
G4Mutex OMSimPrimaryGeneratorAction::m_mutex;

// Unpolarized light: each primary optical photon gets a random polarization perpendicular to its own direction. A
// random 3D direction (as before) leaves a component along the momentum, and Geant4 then treats about two thirds of the
// light as p-polarized. A polarization set explicitly (e.g. /gps/polarization) is kept.
static void setUnpolarized(G4Event *p_event)
{
    for (G4int i = 0; i < p_event->GetNumberOfPrimaryVertex(); ++i)
        for (G4PrimaryParticle *particle = p_event->GetPrimaryVertex(i)->GetPrimary(); particle; particle = particle->GetNext())
        {
            if (particle->GetParticleDefinition() != G4OpticalPhoton::Definition() || particle->GetPolarization().mag2() > 0)
                continue;
            const G4ThreeVector direction = particle->GetMomentumDirection().unit();
            const G4ThreeVector e1 = direction.orthogonal().unit();
            const G4ThreeVector e2 = direction.cross(e1);
            const G4double phi = CLHEP::twopi * G4UniformRand();
            const G4ThreeVector polarization = std::cos(phi) * e1 + std::sin(phi) * e2;
            particle->SetPolarization(polarization.x(), polarization.y(), polarization.z());
        }
}

OMSimPrimaryGeneratorAction::OMSimPrimaryGeneratorAction()
    : m_useParticleGun(true)
{
    OMSimCommandArgsTable& args = OMSimCommandArgsTable::getInstance();

    // Check if a macro file is specified - if so, use GPS for flexibility
    std::string macroFile = args.get<std::string>("macro");
    if (!macroFile.empty()) {
        m_useParticleGun = false;
        log_info("Using GPS mode (macro file specified)");
    }

    if (m_useParticleGun) {
        setupParticleGun();
    } else {
        if (!m_particleSource) {
            m_particleSource = std::make_unique<G4GeneralParticleSource>();
        }
    }
}

void OMSimPrimaryGeneratorAction::setupParticleGun()
{
    OMSimCommandArgsTable& args = OMSimCommandArgsTable::getInstance();

    // Get parameters
    G4double impactParameter = args.get<G4double>("impact_parameter");
    G4double energy = args.get<G4double>("primary_energy");
    std::string particleName = args.get<std::string>("primary_particle");

    // Calculate position: muon at z=r, starting at x=-1.7*r, traveling +x
    // The factor 1.7 ensures the Cherenkov cone (41° in ice) hits the DOM
    G4double startX = -1.7 * impactParameter * m;
    m_position = G4ThreeVector(startX, 0.0, impactParameter * m);
    m_direction = G4ThreeVector(1.0, 0.0, 0.0);

    log_info("Impact parameter mode: d = {} m", impactParameter);
    log_info("Particle: {}, Energy: {} GeV", particleName, energy);
    log_info("Start position: ({}, {}, {}) m",
             m_position.x()/m, m_position.y()/m, m_position.z()/m);
    log_info("Direction: ({}, {}, {})",
             m_direction.x(), m_direction.y(), m_direction.z());

    // Create particle gun
    if (!m_particleGun) {
        m_particleGun = std::make_unique<G4ParticleGun>(1);
    }

    // Set particle type
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle(particleName);
    if (!particle) {
        log_error("Particle '{}' not found! Using mu-", particleName);
        particle = particleTable->FindParticle("mu-");
    }

    m_particleGun->SetParticleDefinition(particle);
    m_particleGun->SetParticleEnergy(energy * GeV);
    m_particleGun->SetParticlePosition(m_position);
    m_particleGun->SetParticleMomentumDirection(m_direction);
}

OMSimPrimaryGeneratorAction::~OMSimPrimaryGeneratorAction()
{
}

void OMSimPrimaryGeneratorAction::GeneratePrimaries(G4Event* p_event)
{
    std::lock_guard<G4Mutex> lock(m_mutex);

    if (m_useParticleGun && m_particleGun) {
        m_particleGun->GeneratePrimaryVertex(p_event);
    } else if (m_particleSource) {
        m_particleSource->GeneratePrimaryVertex(p_event);
        setUnpolarized(p_event);
    }
}
