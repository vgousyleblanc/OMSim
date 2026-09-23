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
            // Glue a gelpad cone onto the PMT's tip - same shape as POM's polar gelpads. With no
            // glass sphere to clip against here, the far face is just left flat.
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

            // Subtract the (translated) PMT solid so the gel wraps the dome instead of
            // overlapping it - PMT's own frame, seen from the cone's local frame, sits at
            // pmtWorldPosition - gelpadPosition.
            G4VSolid *gelpadFinal = new G4SubtractionSolid("SinglePMT_Gelpad", gelpadCone,
                                                            managerPMT->getPMTSolid(), 0, pmtWorldPosition - gelpadPosition);

            G4LogicalVolume *gelpadLV = new G4LogicalVolume(gelpadFinal, m_data->getMaterial("RiAbs_Gel_Shin-Etsu"), "SinglePMT_GelpadLV");
            new G4PVPlacement(0, gelpadPosition, gelpadLV, "SinglePMT_GelpadPhys", m_worldLogical, false, 0,
                               OMSimCommandArgsTable::getInstance().get<bool>("check_overlaps"));
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
