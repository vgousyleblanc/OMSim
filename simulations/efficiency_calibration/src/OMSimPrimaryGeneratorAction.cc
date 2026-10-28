#include "OMSimPrimaryGeneratorAction.hh"
#include <cmath>
#include <Randomize.hh>
#include <G4OpticalPhoton.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>

#include <G4GeneralParticleSource.hh>
#include <G4ParticleTypes.hh>
#include <G4RandomTools.hh>

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
{
	if (!m_particleSource)
	{
	m_particleSource = std::make_unique<G4GeneralParticleSource>();
	}
}

OMSimPrimaryGeneratorAction::~OMSimPrimaryGeneratorAction()
{

}
 
void OMSimPrimaryGeneratorAction::GeneratePrimaries(G4Event *p_event)
{
	if (m_particleSource)
	{
		std::lock_guard<G4Mutex> lock(m_mutex);
		m_particleSource->GeneratePrimaryVertex(p_event);
		setUnpolarized(p_event);
	}
}
