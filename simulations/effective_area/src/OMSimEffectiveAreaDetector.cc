#include "OMSimEffectiveAreaDetector.hh"
#include "OMSimPDOM.hh"
#include "OMSimPOM.hh"
#include "OMSimLOM16.hh"
#include "OMSimLOM18.hh"
#include "OMSimDEGG.hh"
#include "OMSimMDOM.hh"
#include "OMSimCommandArgsTable.hh"
#include "OMSimHitManager.hh"
#include <G4Orb.hh>
#include "OMSimSensitiveDetector.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include <G4Cons.hh>
#include <G4SubtractionSolid.hh>
#include <G4Sphere.hh>
#include <G4IntersectionSolid.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4VisAttributes.hh>

/**
 * @brief Constructs the world volume (sphere).
 */
void OMSimEffectiveAreaDetector::constructWorld()
{
    m_worldSolid = new G4Orb("World", OMSimCommandArgsTable::getInstance().get<G4double>("world_radius") * m);
    m_worldLogical = new G4LogicalVolume(m_worldSolid, m_data->getMaterial("argWorld"), "World_log", 0, 0, 0);
    m_worldPhysical = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), m_worldLogical, "World_phys", 0, false, 0);
    G4VisAttributes *worldVis = new G4VisAttributes(G4Colour(0.45, 0.5, 0.35, 0.));
    m_worldLogical->SetVisAttributes(worldVis);
}

/**
 * @brief Constructs the selected detector from the command line argument and returns the physical world volume.
 * @return Pointer to the physical world volume
 */
void OMSimEffectiveAreaDetector::constructDetector()
{

    OMSimHitManager &hitManager = OMSimHitManager::getInstance();

    bool placeHarness = OMSimCommandArgsTable::getInstance().get<bool>("place_harness");

    OMSimOpticalModule *opticalModule = nullptr;

    switch (OMSimCommandArgsTable::getInstance().get<G4int>("detector_type"))
    {

    case 0:
    {
        log_critical("No custom detector implemented!");
        break;
    }
    case 1:
    {
        log_info("Constructing single PMT");
        const bool placeSinglePMTGelpad = OMSimCommandArgsTable::getInstance().get<bool>("single_pmt_gelpad");

        // Reuse the exact same gelpad geometry as the real P-OM (POM::m_gelpad_*, tunable in
        // OMSimPOM.hh) rather than a separate set of numbers, so this is genuinely the P-OM's
        // gelpad+PMT pair as a standalone entity. The cone's own envelope (small/large radius,
        // thickness) never changes - only the PMT moves. Pushing the PMT forward into the (fixed)
        // cone by m_gelpad_thickness-m_gelThicknessFrontPMT means its tip ends up exactly
        // m_gelThicknessFrontPMT short of the cone's far face, and the subtraction below carves
        // out whatever the PMT's dome actually occupies, so gel and PMT touch directly (no air)
        // everywhere the dome protrudes into the cone.
        const G4double kPMTTranslation = placeSinglePMTGelpad ? (POM::m_gelpad_thickness - POM::m_gelThicknessFrontPMT) : 0.0 * mm;

        OMSimPMTConstruction *managerPMT = new OMSimPMTConstruction();
        managerPMT->selectPMT("argPMT");
        managerPMT->includeHAcoating();
        managerPMT->construction();
        managerPMT->placeIt(G4ThreeVector(0, 0, kPMTTranslation), G4RotationMatrix(), m_worldLogical, "_0");
        hitManager.setNumberOfPMTs(1, 0);
        managerPMT->configureSensitiveVolume(this, "/PMT/0");

        if (placeSinglePMTGelpad)
        {
            // Glue a gelpad onto the PMT's tip - same shape as POM's polar gelpads, including
            // the outer face being cut by the same glass-sphere radius (POM::m_gelpad_sphere_radius),
            // not left flat. A flat exposed face has a fixed surface normal, so light arriving
            // off-axis hits it at a steep local angle and suffers heavy Fresnel reflection there -
            // the real P-OM avoids this because the gelpad's outer curve matches the glass dome
            // it sits against, so it presents close to normal incidence across a wide angular
            // range, same as the bare PMT's own dome does.
            G4VSolid *gelpadCone = new G4Cons("SinglePMT_Gelpad_cone",
                                               0, POM::m_gelpad_small_radius,
                                               0, POM::m_gelpad_large_radius,
                                               POM::m_gelpad_thickness / 2.0,
                                               0, 2 * CLHEP::pi);

            // Cone position is fixed, computed as if the PMT were still at the origin (matching
            // the untranslated tip distance) - it does not follow the PMT's translation.
            G4double pmtOffset = managerPMT->getDistancePMTCenterToTip();
            G4ThreeVector gelpadPosition(0, 0, pmtOffset + POM::m_gelpad_thickness / 2.0);
            G4ThreeVector pmtWorldPosition(0, 0, kPMTTranslation);

            // In POM, the gelpad's far face sits exactly on the glass sphere's surface
            // (rBaseGelpad + thickness/2 == m_glassInRad by construction) - the sphere is huge
            // (201.9mm) compared to the gelpad's own ~66mm width, so it only actually curves
            // anything if the gelpad sits right at its surface, not buried deep inside it.
            // Center the sphere on-axis so it passes through the cone's far face on-axis point.
            G4double farFaceZ = pmtOffset + POM::m_gelpad_thickness;
            G4ThreeVector sphereCenterWorld(0, 0, farFaceZ - POM::m_gelpad_sphere_radius);
            G4VSolid *gelpadSphereCut = new G4Sphere("SinglePMT_Gelpad_sphere_cut",
                                                      0, POM::m_gelpad_sphere_radius,
                                                      0, 2 * CLHEP::pi,
                                                      0, CLHEP::pi);
            G4Transform3D sphereTransform(G4RotationMatrix(), sphereCenterWorld - gelpadPosition);
            G4VSolid *gelpadShaped = new G4IntersectionSolid("SinglePMT_Gelpad_shaped", gelpadCone, gelpadSphereCut, sphereTransform);

            // Subtract the (translated) PMT solid so the gel wraps the dome instead of
            // overlapping it - PMT's own frame, seen from the cone's local frame, sits at
            // pmtWorldPosition - gelpadPosition.
            G4VSolid *gelpadFinal = new G4SubtractionSolid("SinglePMT_Gelpad", gelpadShaped,
                                                            managerPMT->getPMTSolid(), 0, pmtWorldPosition - gelpadPosition);

            G4LogicalVolume *gelpadLV = new G4LogicalVolume(gelpadFinal, m_data->getMaterial("RiAbs_Gel_Shin-Etsu"), "SinglePMT_GelpadLV");
            G4VPhysicalVolume *gelpadPhys = new G4PVPlacement(0, gelpadPosition, gelpadLV, "SinglePMT_GelpadPhys", m_worldLogical, false, 0,
                               OMSimCommandArgsTable::getInstance().get<bool>("check_overlaps"));

            if (OMSimCommandArgsTable::getInstance().get<bool>("single_pmt_reflector"))
            {
                // A thin conical shell wrapping just the gel's tapered side wall (not the front/
                // back faces, which must stay open for light to enter/exit) - tests whether a
                // real reflective coating there captures more light than relying purely on the
                // gel's own total-internal-reflection funneling at that same boundary.
                const G4double kShellThickness = 1.0 * mm;
                G4VSolid *reflectorSleeve = new G4Cons("SinglePMT_ReflectorSleeve",
                                                        POM::m_gelpad_small_radius, POM::m_gelpad_small_radius + kShellThickness,
                                                        POM::m_gelpad_large_radius, POM::m_gelpad_large_radius + kShellThickness,
                                                        POM::m_gelpad_thickness / 2.0,
                                                        0, 2 * CLHEP::pi);
                G4LogicalVolume *reflectorLV = new G4LogicalVolume(reflectorSleeve, m_data->getMaterial("NoOptic_Reflector"), "SinglePMT_ReflectorLV");
                reflectorLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.7, 0.7, 0.75)));
                G4VPhysicalVolume *reflectorPhys = new G4PVPlacement(0, gelpadPosition, reflectorLV, "SinglePMT_ReflectorPhys", m_worldLogical, false, 0,
                                                                      OMSimCommandArgsTable::getInstance().get<bool>("check_overlaps"));

                auto mirrorSurface = m_data->getOpticalSurface("Surf_PMTSideMirror");
                new G4LogicalBorderSurface("SinglePMT_Gel_to_Reflector", gelpadPhys, reflectorPhys, mirrorSurface);
                new G4LogicalBorderSurface("SinglePMT_Reflector_to_Gel", reflectorPhys, gelpadPhys, mirrorSurface);
            }
        }
        break;
    }
    case 2:
    {

        opticalModule = new mDOM(placeHarness);
        break;
    }
    case 3:
    {

        opticalModule = new DOM(placeHarness);
        break;
    }
    case 4:
    {

        opticalModule = new LOM16(placeHarness);
        break;
    }
    case 5:
    {

        opticalModule = new LOM18(placeHarness);
        break;
    }
    case 6:
    {
        opticalModule = new DEGG(placeHarness);
        break;
    }
    case 7:
    {
        opticalModule = new DOM(placeHarness, true);
        break;
    }
    case 8:
    {
        opticalModule = new POM(placeHarness);
        break;
    }
    }

    if (opticalModule)
    {
        opticalModule->placeIt(G4ThreeVector(0, 0, 0), G4RotationMatrix(), m_worldLogical, "");
        opticalModule->configureSensitiveVolume(this);
    }
}
