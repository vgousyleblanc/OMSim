/** 
 *  @todo   - Clean up magic number variables and comment their meaning
 *          - Write documentation and parse current comments into Doxygen style
 */
#include "OMSimPOM_4inch.hh"
#include "OMSimTools.hh"
#include "OMSimLogger.hh"
#include "OMSimCommandArgsTable.hh"

#include <G4Cons.hh>
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include <G4Ellipsoid.hh>
#include <G4IntersectionSolid.hh>
#include <G4DisplacedSolid.hh>
#include <G4UnionSolid.hh>
#include <G4Polycone.hh>
#include <G4EllipticalCone.hh>
#include <G4TessellatedSolid.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4OpticalSurface.hh>
#include <G4MaterialPropertiesTable.hh>
#include <fstream>
#include <cstdlib>
#include <G4Torus.hh>

// Gelpad is the mother volume of its PMT (see createGelpadLogicalVolumes()); POM_GEL_SUBTRACT=1
// reverts to the old gel-minus-PMT subtraction.
static const bool m_gelMotherOfPMT = !std::getenv("POM_GEL_SUBTRACT");

// Gel thickness in front of every PMT tip. All 8 PMTs sit at the same radius: in the P-OM CAD
// (P-OM_internal_v13) every PMT tip is the same distance below the inner glass surface.
// m_EqPMTrOffset used to be 2.6mm, copied from the mDOM (whose equatorial PMTs really do sit further
// out); in the P-OM it pushed the outer ring's PMTs and gelpads 2.6mm closer to the glass than the
// polar ring's, leaving the outer PMTs with 2.6mm less front gel.
// Pad dimensions follow the 4-inch CAD ("gel pad 4 inch 35deg" of the P2DOR v1 hemisphere assembly, measured from
// its STL): back face at the PMT's widest point with radius 52.4 mm, wall 35deg from the PMT axis, 31.4 mm from the
// back face to the inner glass on the axis (rim meets the glass at ~66 mm radius), 0.85 mm gel in front of the PMT
// tip. In the CAD every pad/PMT axis points at the dome centre, both rings at the same radius (32.5 / 65 deg).
// They apply to both hemispheres of the full module and the single-hemisphere prototype.
// TEST: POM_GEL_FRONT_MM overrides the 0.85 mm gel thickness in front of the PMT.
const G4double POM_4inch::m_gelThicknessFrontPMT = (std::getenv("POM_GEL_FRONT_MM") ? std::atof(std::getenv("POM_GEL_FRONT_MM")) : 0.85) * mm;
// TEST: POM_PAD_SMALL_RADIUS_MM overrides the 52.4 mm small (back-face) radius.
const G4double POM_4inch::m_gelpad_small_radius = (std::getenv("POM_PAD_SMALL_RADIUS_MM") ? std::atof(std::getenv("POM_PAD_SMALL_RADIUS_MM")) : 52.4) * mm;
// TEST: POM_PAD_THICKNESS_MM overrides the 31.4 mm pad thickness (back face to the inner glass on the axis).
const G4double POM_4inch::m_gelpad_thickness = (std::getenv("POM_PAD_THICKNESS_MM") ? std::atof(std::getenv("POM_PAD_THICKNESS_MM")) : 31.4) * mm;
// TEST: POM_PAD_OPENING_DEG overrides the 35deg opening angle.
const G4double POM_4inch::m_polarPadOpeningAngle = (std::getenv("POM_PAD_OPENING_DEG") ? std::atof(std::getenv("POM_PAD_OPENING_DEG")) : 35.0) * deg;
const G4double POM_4inch::m_gelpad_large_radius = POM_4inch::m_gelpad_small_radius + std::tan(POM_4inch::m_polarPadOpeningAngle) * POM_4inch::m_gelpad_thickness;
const G4double POM_4inch::m_gelpad_sphere_radius = 201.9 * mm; // == m_glassInRad

POM_4inch::POM_4inch(G4bool p_placeHarness) : OMSimOpticalModule(new OMSimPMTConstruction()), m_placeHarness(p_placeHarness),
    // TEST: POM_GLASS_EXTEND_MM (single hemisphere only) inserts a straight glass section of that
    // height between the PVC band (top stays at z=40mm) and the dome, lifting dome, gelpads and
    // PMTs by the same amount while the frame keeps its absolute position.
    m_cylinderHeight(OMSimCommandArgsTable::getInstance().get<bool>("single_hemisphere")
                         ? 40.0 * mm + (std::getenv("POM_GLASS_EXTEND_MM") ? std::atof(std::getenv("POM_GLASS_EXTEND_MM")) : 0.0) * mm
                         : 80.0 * mm),
    m_singleHemisphere(OMSimCommandArgsTable::getInstance().get<bool>("single_hemisphere")),
    m_pomReflector(OMSimCommandArgsTable::getInstance().get<bool>("pom_reflector")),
    m_numberBuiltPMTs(m_singleHemisphere ? (m_numberPolarPMTs + m_numberEqPMTs) : m_totalNumberPMTs)
{
    log_info("Constructing POM");
    // The P-OM's PMTs are clear glass all round: no black HA coating on the back tube (unlike the mDOM's).
    // POM_HA_COATING=1 puts it back for comparison.
    if (std::getenv("POM_HA_COATING"))
        m_managerPMT->includeHAcoating();
    m_managerPMT->selectPMT(m_PMTModel);
    m_managerPMT->construction();
    m_PMToffset = m_managerPMT->getDistancePMTCenterToTip(); //
    m_maxPMTRadius = m_managerPMT->getMaxPMTRadius() + 2 * mm;//43
    std::cout << "PMT offset: " << m_PMToffset << std::endl;
    //if (p_placeHarness) m_harness = new POMHarness(this);
    construction();
    //if (p_placeHarness) integrateDetectorComponent(m_harness, G4ThreeVector(0, 0, 0), G4RotationMatrix(), "");
    log_trace("Finished constructing POM");
}

void POM_4inch::PrintVolumeTree(G4LogicalVolume* lv, G4int depth)
{
    if (!lv) return;

    for (G4int i = 0; i < depth; ++i)
        G4cout << "  ";

    G4cout << lv->GetName()
           << "  material=" << lv->GetMaterial()->GetName()
           << "  daughters=" << lv->GetNoDaughters()
           << G4endl;

    for (G4int i = 0; i < lv->GetNoDaughters(); ++i)
    {
        auto pv = lv->GetDaughter(i);
        for (G4int j = 0; j < depth + 1; ++j)
            G4cout << "  ";
        G4cout << "-> " << pv->GetName() << G4endl;
        PrintVolumeTree(pv->GetLogicalVolume(), depth + 1);
    }
}

void POM_4inch::construction()
{
    // The optical hierarchy is world (water) -> glass shell -> inner cavity -> gelpads/PMTs.
    //
    // This used to subtract the inner cavity's own solid from the outer glass solid to build an
    // explicit hollow shell (outerGlass MINUS innerAir), with innerAirLV then placed as a daughter
    // of that shell. That makes the shell and the cavity share an exactly-tangent boundary
    // everywhere (same shape, same radius, nothing in between) - a textbook G4 boolean degeneracy,
    // and not just a harmless --check_overlaps false positive here: photon tracking confirmed a
    // real navigation crack from it (rays crossing the glass near the dome's pole were exiting
    // into stray World material instead of reaching InnerAirPhys, some then re-entering the cavity
    // near the OTHER hemisphere's pole, i.e. crossing between hemispheres through this crack - not
    // through any real optical path). Padding the subtracted copy slightly larger (tried 0.05mm and
    // 2mm) didn't fix it either - it just turned the exact-tangent crack into a literal gap of the
    // same width, since a GEOMETRIC VOID now existed between the (unpadded) cavity daughter's own
    // boundary and the (padded) hole it was supposed to exactly fill.
    //
    // The actual fix: don't subtract the cavity from the glass at all. A solid glass sphere with
    // the air cavity placed as a plain daughter inside it is optically identical (daughters always
    // override their mother's material in the region they occupy, so the mother's own material -
    // glass - automatically fills everywhere else) and has no fragile double-boolean chain to get
    // wrong. The titanium ring still needs its own notch cut from the glass below (a real
    // geometric necessity, not a side effect of this pattern), so glassShellSolid is still built
    // with one subtraction, just one fewer than before.
    G4VSolid *innerAirSolid = pressureVessel(m_glassInRad, "InnerAirVoid", !m_singleHemisphere);
    G4VSolid *glassShellSolid = pressureVessel(m_glassOutRad, "GlassShellOuter", !m_singleHemisphere);

    // The cylinder band below the dome is never glass. Full module: a titanium band between the
    // two domes. Single-hemisphere prototype: the band is the (opaque) PVC frame the dome is
    // mounted on.
    G4VSolid *titaniumCylinder = nullptr;
    G4LogicalVolume *titaniumLV = nullptr;
    G4VSolid *pvcBand = nullptr;
    G4LogicalVolume *pvcBandLV = nullptr;
    {
        // The band is physically thicker than the glass, so its outer radius
        // (m_cylinder_outer_radius) extends past m_glassOutRad. That means it cannot be a
        // daughter of the glass shell (a daughter can never exceed its mother's own bounds) - it
        // sits alongside the shell as a sibling, both placed directly in the world, and the
        // glass shell has the same band carved out of it.
        // The cut (and the band itself) starts kCutPadding OUTSIDE the cavity radius rather than
        // exactly on it: exactly on it is a G4 boolean degeneracy (shared surface with
        // innerAirSolid), and starting inside it left the air cavity poking 0.01mm out of its
        // glass mother (an InnerAirPhys/GlassShellLV overlap). The 0.01mm glass skin left between
        // the cavity and the band is optically irrelevant.
        const G4double kCutPadding = 0.01 * mm;
        // With POM_GLASS_EXTEND_MM the band stops glassExtend below the dome (at z=40mm), leaving
        // that top part of the cylinder as glass wall; bandShift recentres the shortened band.
        const G4double glassExtend = (m_singleHemisphere && std::getenv("POM_GLASS_EXTEND_MM")) ? std::atof(std::getenv("POM_GLASS_EXTEND_MM")) * mm : 0.0;
        m_bandShift = G4ThreeVector(0, 0, -glassExtend / 2.0);
        // TEST (POM_PVC_OUTER_MM, single hemisphere): outer radius of the PVC flange/back plate below the equator
        // (default m_cylinder_outer_radius = 229 mm, i.e. only 13 mm beyond the glass).
        const G4double bandOuter = (m_singleHemisphere && std::getenv("POM_PVC_OUTER_MM")) ? std::atof(std::getenv("POM_PVC_OUTER_MM")) * mm : m_cylinder_outer_radius;
        G4Tubs *band = new G4Tubs("CylinderBand", m_glassInRad + kCutPadding, bandOuter, m_cylinderHeight - glassExtend / 2.0, 0, 2 * CLHEP::pi);
        G4Tubs *bandCutter = new G4Tubs("CylinderBandCutter", m_glassInRad + kCutPadding, m_cylinder_outer_radius,
                                        m_cylinderHeight - glassExtend / 2.0 + kCutPadding, 0, 2 * CLHEP::pi);
        glassShellSolid = new G4SubtractionSolid("GlassShell", glassShellSolid, bandCutter, 0, m_bandShift);

        if (m_singleHemisphere)
        {
            pvcBand = band;
            pvcBandLV = new G4LogicalVolume(pvcBand, m_data->getMaterial("Plastic"), "PVCFrameBandLV");
        }
        else
        {
            titaniumCylinder = band;
            titaniumLV = new G4LogicalVolume(titaniumCylinder, m_data->getMaterial("Titanium"), "TitaniumCylinderLV");
        }
    }

    setPMTAndGelpadPositions();

    auto glassMat = m_data->getMaterial("RiAbs_Glass_Vitrovex");
    // TEST: POM_GLASS_ABS_MM replaces the pressure-vessel glass absorption length by a constant
    // value (mm) at all wavelengths (Vitrovex data: ~380 mm at 460 nm).
    if (const char *absEnv = std::getenv("POM_GLASS_ABS_MM"))
    {
        const G4double absLength = std::atof(absEnv) * mm;
        const std::vector<G4double> energies = {1.0 * eV, 7.0 * eV};
        glassMat->GetMaterialPropertiesTable()->AddProperty("ABSLENGTH", energies, {absLength, absLength}, true);
        log_info("Glass absorption length set to {} mm", absLength / mm);
    }
    // TEST: POM_GEL_ABS_MM replaces the gelpad gel's absorption length by a constant value (mm) at all
    // wavelengths (P-ONE gel data: 420 mm at 400 nm, 549 mm at 460 nm).
    if (const char *gelAbsEnv = std::getenv("POM_GEL_ABS_MM"))
    {
        const G4double absLength = std::atof(gelAbsEnv) * mm;
        const std::vector<G4double> energies = {1.0 * eV, 7.0 * eV};
        m_data->getMaterial("argGel")->GetMaterialPropertiesTable()->AddProperty("ABSLENGTH", energies, {absLength, absLength}, true);
        log_info("Gel absorption length set to {} mm", absLength / mm);
    }
    // TEST: scattering inside the gelpad gel. POM_GEL_RAYLEIGH_MM = Rayleigh scattering length (mm);
    // POM_GEL_MIE_MM = Mie (Henyey-Greenstein) scattering length (mm), with POM_GEL_MIE_G = forward asymmetry g (default 0.9).
    if (const char *rayEnv = std::getenv("POM_GEL_RAYLEIGH_MM"))
    {
        const std::vector<G4double> energies = {1.0 * eV, 7.0 * eV};
        const G4double length = std::atof(rayEnv) * mm;
        m_data->getMaterial("argGel")->GetMaterialPropertiesTable()->AddProperty("RAYLEIGH", energies, {length, length}, true);
        log_info("Gel Rayleigh scattering length set to {} mm", length / mm);
    }
    if (const char *mieEnv = std::getenv("POM_GEL_MIE_MM"))
    {
        const std::vector<G4double> energies = {1.0 * eV, 7.0 * eV};
        const G4double length = std::atof(mieEnv) * mm;
        const G4double g = std::getenv("POM_GEL_MIE_G") ? std::atof(std::getenv("POM_GEL_MIE_G")) : 0.9;
        auto mpt = m_data->getMaterial("argGel")->GetMaterialPropertiesTable();
        mpt->AddProperty("MIEHG", energies, {length, length}, true);
        mpt->AddConstProperty("MIEHG_FORWARD", g, true);
        mpt->AddConstProperty("MIEHG_BACKWARD", g, true);
        mpt->AddConstProperty("MIEHG_FORWARD_RATIO", 1.0, true);
        log_info("Gel Mie scattering length set to {} mm (g = {})", length / mm, g);
    }
    auto airMat = m_data->getMaterial("Ri_Air");

    G4LogicalVolume *glassShellLV = new G4LogicalVolume(glassShellSolid, glassMat, "GlassShellLV");
    G4LogicalVolume *innerAirLV = new G4LogicalVolume(innerAirSolid, airMat, "InnerAirLV");

    // Place the cavity inside the glass shell. The titanium ring (full-module mode only) is NOT
    // placed here - it is wider than the glass shell (see comment above) and is instead appended
    // as its own top-level component below, so it gets placed as a sibling of the shell in the world.
    new G4PVPlacement(nullptr, G4ThreeVector(), innerAirLV, "InnerAirPhys", glassShellLV, false, 0, m_checkOverlaps);

    // TEST (POM_GLASS_SKIN_DEG = polar angle from the pole, deg): a 0.2 mm glass skin on the outside of
    // the dome from that angle down to the equator, so that a different outer finish can be given to the
    // lower part of the glass only (POM_GLASS_SKIN_SIGMA, applied in OMSimEffectiveAreaDetector).
    if (const char *skinEnv = std::getenv("POM_GLASS_SKIN_DEG"))
    {
        const G4double thetaStart = std::atof(skinEnv) * deg;
        const G4double thetaEnd = CLHEP::halfpi - std::asin(0.05 * mm / m_glassOutRad);
        G4Sphere *skin = new G4Sphere("GlassSkin", m_glassOutRad - 0.2 * mm, m_glassOutRad, 0, 2 * CLHEP::pi, thetaStart, thetaEnd - thetaStart);
        G4LogicalVolume *skinLV = new G4LogicalVolume(skin, glassMat, "GlassSkinLV");
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, m_cylinderHeight), skinLV, "GlassSkinPhys", glassShellLV, false, 0, m_checkOverlaps);
        log_info("Glass outer skin from {} deg to the equator", thetaStart / deg);
    }

    // TEST: POM_GEL_GLASS_GAP_MM leaves a thin air layer between every gelpad's front face and the
    // glass (gelpads clipped to a cavity shrunk by that amount) - models imperfect gel-glass contact.
    if (const char *gapEnv = std::getenv("POM_GEL_GLASS_GAP_MM"))
    {
        const G4double gap = std::atof(gapEnv) * mm;
        createGelpadLogicalVolumes(pressureVessel(m_glassInRad - gap, "GelClipShrunk", !m_singleHemisphere));
        log_info("Gelpads pulled back {} mm from the glass (air layer)", gap / mm);
    }
    else
        createGelpadLogicalVolumes(innerAirSolid);
    placePMTs(innerAirLV);
    placeGelpads(innerAirLV);
    if (m_pomReflector)
        placeGelpadReflectors(innerAirLV, innerAirSolid);
    InternalCADComponents(innerAirLV);

    if (!m_singleHemisphere)
    {
        // Photon tracking showed the glass sphere + gel-filled cavity acts as a thick ball lens:
        // parallel light entering one hemisphere's glass refracts to a focus a few hundred mm
        // away, and for marginal rays that focus sits close enough to the module that light can
        // cross the equator into the OTHER hemisphere's air cavity and reach its PMTs.
        //
        // The real module seals this with two solid plates, one per hemisphere, closing off each
        // glass half's open base completely - not just a ring around the electronics board (an
        // earlier, more conservative version of this only filled the open annulus around the
        // board, which left the board's own internal gaps - mounting holes, cutouts - exposed;
        // photon tracking confirmed light still leaked through those).
        //
        // An even earlier version tried to avoid colliding with the existing electronics board
        // and PMT frames by subtracting their real CAD-mesh shapes out of each plate first. That
        // back-fired: CAD_Frame_up/down's mesh, once correctly transformed into world coordinates
        // (a 90+90deg rotation - easy to get backwards, and the fix still wasn't enough), turns
        // out to have a genuinely large footprint that overlaps almost this entire disk, so the
        // subtraction voided out nearly the whole plate - confirmed with a same-seed, same-RNG
        // before/after photon-tracking comparison showing the "sealed" plate changed the outcome
        // of only 1 ray out of 300, i.e. it was blocking almost nothing.
        //
        // The electronics board and both frames are already NoOptic_Absorber (fully opaque), so a
        // plain full disk placed right through them is not a physics problem, only a geometry-QA
        // one: --check_overlaps will flag the overlap, but two coincident absorbers behave exactly
        // like one - no photon can tell the difference. Simplicity and an actual working seal win
        // here over a "clean" geometry that seals almost nothing.
        const G4double kPlateHalfHeight = 2.0 * mm;
        const G4Colour kPlateColour(0.9, 0.6, 0.1, 1.0); // distinct from the dark-gray frame/board, easy to spot in the viewer
        G4VisAttributes plateVis(kPlateColour);

        for (G4int side : {+1, -1})
        {
            // The two plates must touch exactly at z=0 with zero gap between them - an earlier
            // version left a 1mm gap "for realism" (sitting back-to-back, not literally touching)
            // and that alone let the large majority of cross-hemisphere rays leak straight through
            // the resulting 2mm-wide open-air slice spanning the full radius at the equator.
            const G4ThreeVector plateCenter(0, 0, side * kPlateHalfHeight);
            G4VSolid *plateSolid = new G4Tubs("EquatorPlate_raw", 0, m_glassInRad, kPlateHalfHeight, 0, 2 * CLHEP::pi);
            G4String name = "EquatorPlate_" + std::string(side > 0 ? "up" : "down");
            G4LogicalVolume *plateLV = new G4LogicalVolume(plateSolid, m_data->getMaterial("NoOptic_Absorber"), name + "LV");
            plateLV->SetVisAttributes(plateVis);
            new G4PVPlacement(nullptr, plateCenter, plateLV, name + "Phys", innerAirLV, false, 0, false);
        }
    }

    if (m_singleHemisphere)
    {
        // Omitting the lower hemisphere leaves the lower half of the
        // cavity (z = -m_cylinderHeight..0, radius 0..m_glassInRad) as air, which doesn't match
        // reality: that volume is occupied by the mounting flange, glued directly onto the
        // glass tube's inner surface and its open bottom rim - no air gap. Fill it with a solid
        // plastic block so it touches the glass directly on both the flat rim (bottom) and the
        // cylindrical wall (sides); only the annulus between m_glassInRad and m_glassOutRad
        // remains real glass tube wall, unaffected. The block stops short of the electronics
        // board (still placed at the equator regardless of hemisphere count), rather than
        // reaching all the way to z=0, since that board genuinely extends below the equator.
        G4double flangeTopZ = 0.0;
        if (auto electronics = G4PhysicalVolumeStore::GetInstance()->GetVolume("CAD_Electronics_physical"))
        {
            G4ThreeVector boundMin, boundMax;
            electronics->GetLogicalVolume()->GetSolid()->BoundingLimits(boundMin, boundMax);
            const G4double kSafetyGap = 2.0 * mm;
            flangeTopZ = std::min(flangeTopZ, boundMin.z() - kSafetyGap);
        }
        const G4double kFlangeHalfHeight = (flangeTopZ + m_cylinderHeight) / 2.0;
        const G4double kFlangeCenterZ = -m_cylinderHeight + kFlangeHalfHeight;
        G4Tubs *capSolid = new G4Tubs("HemisphereCap", 0, m_glassInRad, kFlangeHalfHeight, 0, 2 * CLHEP::pi);
        G4LogicalVolume *capLV = new G4LogicalVolume(capSolid, m_data->getMaterial("Plastic"), "HemisphereCapLV");
        capLV->SetVisAttributes(m_pom_frame);
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, kFlangeCenterZ), capLV, "HemisphereCapPhys", innerAirLV, false, 0, m_checkOverlaps);

        // TEST (POM_INTERIOR_REFL = 0-1): the PVC inside the module (bottom plate + flange wall) as a
        // diffuse (Lambertian) reflector instead of black.
        if (const char *reflEnv = std::getenv("POM_INTERIOR_REFL"))
        {
            auto pvcSurface = new G4OpticalSurface("PVCInteriorSurface", unified, ground, dielectric_metal, 0.0);
            auto mpt = new G4MaterialPropertiesTable();
            const std::vector<G4double> energies = {1.0 * eV, 7.0 * eV};
            const G4double refl = std::atof(reflEnv);
            mpt->AddProperty("REFLECTIVITY", energies, {refl, refl});
            mpt->AddProperty("SPECULARLOBECONSTANT", energies, {0.0, 0.0});
            mpt->AddProperty("SPECULARSPIKECONSTANT", energies, {0.0, 0.0});
            mpt->AddProperty("BACKSCATTERCONSTANT", energies, {0.0, 0.0});
            pvcSurface->SetMaterialPropertiesTable(mpt);
            new G4LogicalSkinSurface("PVCCapSkin", capLV, pvcSurface);
            if (pvcBandLV) new G4LogicalSkinSurface("PVCBandSkin", pvcBandLV, pvcSurface);
            log_info("Interior PVC (bottom plate + flange): diffuse reflector, reflectivity {}", refl);
        }
    }

    // Keep this component registered so the module can be placed in the world later.
    appendComponent(glassShellSolid, glassShellLV, G4ThreeVector(0, 0, 0), G4RotationMatrix(), "PressureVessel_" + std::to_string(m_index));
    if (titaniumLV)
    {
        titaniumLV->SetVisAttributes(m_pom_flange);
        appendComponent(titaniumCylinder, titaniumLV, G4ThreeVector(0, 0, 0), G4RotationMatrix(), "TitaniumFlange_" + std::to_string(m_index));
    }
    // Titanium flange collar of each glass hemisphere, as in the CAD (FLA-TI-17-14, which the real module matches): the
    // glass sits in the flange, whose collar is a cylinder r = 217.2-230 mm (1.3 mm water gap to the glass) rising 6 mm
    // above the glass equator from the band top. The full module has one per hemisphere, both fitting into the
    // titanium ring in the middle (TitaniumFlange); the single-hemisphere prototype has one, on the PVC band.
    // POM_TI_COLLAR_MM overrides the height (0 = no collar); POM_TI_COLLAR_REFL (0-1) makes it a polished reflector.
    // The TEST spherical ring (POM_TI_RING_MM, below) replaces the collar unless a collar height is given as well.
    const G4double collarHeight = (std::getenv("POM_TI_COLLAR_MM") ? std::atof(std::getenv("POM_TI_COLLAR_MM"))
                                   : (std::getenv("POM_TI_RING_MM") ? 0.0 : 6.0)) * mm;
    if (collarHeight > 0)
    {
        G4Tubs *collar = new G4Tubs("TitaniumCollar", 217.2 * mm, 230.0 * mm, collarHeight / 2.0, 0, 2 * CLHEP::pi);
        G4LogicalVolume *collarLV = new G4LogicalVolume(collar, m_data->getMaterial("NoOptic_Absorber"), "TitaniumCollarLV");
        collarLV->SetVisAttributes(m_pom_flange);
        if (const char *reflEnv = std::getenv("POM_TI_COLLAR_REFL"))
        {
            auto tiSurface = new G4OpticalSurface("TitaniumCollarSurface", unified, polished, dielectric_metal);
            auto mpt = new G4MaterialPropertiesTable();
            const std::vector<G4double> energies = {1.0 * eV, 7.0 * eV};
            const G4double refl = std::atof(reflEnv);
            mpt->AddProperty("REFLECTIVITY", energies, {refl, refl});
            tiSurface->SetMaterialPropertiesTable(mpt);
            new G4LogicalSkinSurface("TitaniumCollarSkin", collarLV, tiSurface);
        }
        const G4double collarZ = m_cylinderHeight + collarHeight / 2.0 + 0.01 * mm;
        appendComponent(collar, collarLV, G4ThreeVector(0, 0, collarZ), G4RotationMatrix(), "TitaniumCollar_" + std::to_string(m_index));
        if (!m_singleHemisphere)
            appendComponent(collar, collarLV, G4ThreeVector(0, 0, -collarZ), G4RotationMatrix(), "TitaniumCollarLower_" + std::to_string(m_index));
        log_info("Titanium collar (CAD): {} mm above the glass equator, r = 217.2-230 mm, {}", collarHeight / mm,
                 m_singleHemisphere ? "one" : "one per hemisphere");
    }
    G4VSolid *tiRingSolid = nullptr;   // TEST titanium ring (POM_TI_RING_MM), kept for the cross-section diagnostic
    G4ThreeVector tiRingPos;
    if (pvcBandLV)
    {
        pvcBandLV->SetVisAttributes(m_pom_frame);
        appendComponent(pvcBand, pvcBandLV, m_bandShift, G4RotationMatrix(), "PVCFrameBand_" + std::to_string(m_index));

        // TEST (POM_TI_RING_MM = height in mm above the dome's equator): the titanium ring glued
        // around the bottom of the dome and screwed into the PVC frame - an opaque shell that hugs
        // the outside of the glass (same centre as the dome, just outside its outer radius) from the
        // top of the PVC band (z = m_cylinderHeight) up that far. POM_TI_RING_THICK_MM: radial
        // thickness (default 6 mm).
        if (const char *ringEnv = std::getenv("POM_TI_RING_MM"))
        {
            const G4double ringHeight = std::atof(ringEnv) * mm;
            const char *ringThickEnv = std::getenv("POM_TI_RING_THICK_MM");
            const G4double ringInner = m_glassOutRad + 0.05 * mm;
            const G4double ringOuter = ringInner + (ringThickEnv ? std::atof(ringThickEnv) : 6.0) * mm;
            const G4double thetaTop = std::acos(ringHeight / ringInner);                  // polar angle of the ring's upper edge
            const G4double thetaBottom = CLHEP::halfpi - std::asin(0.01 * mm / ringInner); // 0.01 mm above the band top
            G4Sphere *ring = new G4Sphere("TitaniumRing", ringInner, ringOuter, 0, 2 * CLHEP::pi, thetaTop, thetaBottom - thetaTop);
            G4LogicalVolume *ringLV = new G4LogicalVolume(ring, m_data->getMaterial("NoOptic_Absorber"), "TitaniumRingLV");
            ringLV->SetVisAttributes(m_pom_flange);
            // POM_TI_RING_REFL (0-1): polished reflective titanium surface instead of black.
            if (const char *reflEnv = std::getenv("POM_TI_RING_REFL"))
            {
                auto tiSurface = new G4OpticalSurface("TitaniumRingSurface", unified, polished, dielectric_metal);
                auto mpt = new G4MaterialPropertiesTable();
                const std::vector<G4double> energies = {1.0 * eV, 7.0 * eV};
                const G4double refl = std::atof(reflEnv);
                mpt->AddProperty("REFLECTIVITY", energies, {refl, refl});
                tiSurface->SetMaterialPropertiesTable(mpt);
                new G4LogicalSkinSurface("TitaniumRingSkin", ringLV, tiSurface);
                log_info("Titanium ring surface: polished metal, reflectivity {}", refl);
            }
            tiRingSolid = ring;
            tiRingPos = G4ThreeVector(0, 0, m_cylinderHeight);
            appendComponent(ring, ringLV, tiRingPos, G4RotationMatrix(), "TitaniumRing_" + std::to_string(m_index));
            log_info("Titanium ring hugging the dome: {} mm up the glass, r = {:.1f}-{:.1f} mm", ringHeight / mm, ringInner / mm, ringOuter / mm);
        }
    }

    glassShellLV->SetVisAttributes(m_glassVis);
    innerAirLV->SetVisAttributes(m_airVis);

    // Optional geometry diagnostic, only when the POM_XSECTION environment variable is set: a
    // 1mm-grid Y=0 cross-section of the built geometry via Inside() checks, written to
    // /tmp/pom_crosssection.csv, plus the gel-to-glass gap along each PMT axis.
    // Optional diagnostic, only when POM_FRAME_MAP is set: for every direction (theta, phi) from
    // the dome centre, march radially inward from the inner glass surface and record the first
    // cavity daughter hit (gelpad, frame, ...) and its depth below the glass, or "open" if nothing
    // is hit within POM_FRAME_MAP_DEPTH mm (default 60). Written to /tmp/pom_frame_map.csv; used to locate the frame's openings.
    if (std::getenv("POM_FRAME_MAP"))
    {
        std::ofstream mapFile("/tmp/pom_frame_map.csv");
        mapFile << "theta,phi,label,depth\n";
        const G4double mapDepth = std::getenv("POM_FRAME_MAP_DEPTH") ? std::atof(std::getenv("POM_FRAME_MAP_DEPTH")) : 60.0;
        const G4ThreeVector center(0, 0, m_cylinderHeight);
        for (G4double th = 0; th <= 90; th += 1)
        {
            for (G4double ph = 0; ph < 360; ph += 1)
            {
                const G4ThreeVector dir(sin(th * deg) * cos(ph * deg), sin(th * deg) * sin(ph * deg), cos(th * deg));
                G4String label = "open";
                G4double depth = -1;
                for (G4double d = 0.05; d <= mapDepth && label == "open"; d += 0.5)
                {
                    const G4ThreeVector p = center + (m_glassInRad - d * mm) * dir;
                    if (innerAirSolid->Inside(p) == kOutside) continue;
                    for (G4int i = 0; i < innerAirLV->GetNoDaughters(); ++i)
                    {
                        auto pv = innerAirLV->GetDaughter(i);
                        G4RotationMatrix rotInv;
                        if (pv->GetObjectRotation()) rotInv = pv->GetObjectRotation()->inverse();
                        if (pv->GetLogicalVolume()->GetSolid()->Inside(rotInv * (p - pv->GetObjectTranslation())) != kOutside)
                        {
                            label = pv->GetName();
                            depth = d;
                            break;
                        }
                    }
                }
                mapFile << th << "," << ph << "," << label << "," << depth << "\n";
            }
        }
        std::cout << "DIAG frame map written to /tmp/pom_frame_map.csv" << std::endl;
    }

    if (std::getenv("POM_XSECTION"))
    {
        // POM_XSECTION_PHI (deg, default 0) rotates the cut plane about the module axis; x is then
        // the in-plane coordinate along azimuth phi. POM_XSECTION_STEP (mm, default 1) sets the grid.
        const G4double cutPhi = (std::getenv("POM_XSECTION_PHI") ? std::atof(std::getenv("POM_XSECTION_PHI")) : 0.0) * deg;
        const G4double gridStep = std::getenv("POM_XSECTION_STEP") ? std::atof(std::getenv("POM_XSECTION_STEP")) : 1.0;
        std::ofstream diagFile("/tmp/pom_crosssection.csv");
        diagFile << "x,z,label\n";
        G4double zMin = m_singleHemisphere ? -m_cylinderHeight - 10 : -260;
        G4double zMax = m_cylinderHeight + m_glassOutRad + 10;
        for (G4double x = -240; x <= 240; x += gridStep)
        {
            for (G4double z = zMin; z <= zMax; z += gridStep)
            {
                G4ThreeVector testPoint(x * cos(cutPhi) * mm, x * sin(cutPhi) * mm, z * mm);
                G4String label = "outside_cavity";
                if (innerAirSolid->Inside(testPoint) != kOutside)
                {
                    label = "air";
                    for (G4int i = 0; i < innerAirLV->GetNoDaughters(); ++i)
                    {
                        auto pv = innerAirLV->GetDaughter(i);
                        G4RotationMatrix rotInv;
                        if (pv->GetObjectRotation()) rotInv = pv->GetObjectRotation()->inverse();
                        G4ThreeVector local = rotInv * (testPoint - pv->GetObjectTranslation());
                        if (pv->GetLogicalVolume()->GetSolid()->Inside(local) != kOutside)
                        {
                            label = pv->GetName();
                            break;
                        }
                    }
                }
                else if (glassShellSolid->Inside(testPoint) != kOutside)
                    label = "glass";
                else if (titaniumCylinder && titaniumCylinder->Inside(testPoint) != kOutside)
                    label = "titanium";
                else if (pvcBand && pvcBand->Inside(testPoint - m_bandShift) != kOutside)
                    label = "pvc";
                else if (tiRingSolid && tiRingSolid->Inside(testPoint - tiRingPos) != kOutside)
                    label = "titanium_ring";
                diagFile << x << "," << z << "," << label << "\n";
            }
        }
        std::cout << "DIAG cross-section written to /tmp/pom_crosssection.csv" << std::endl;

        // Fine radial scan from the dome center along each PMT axis (and +-8deg off-axis):
        // where does the gel end and where does the glass begin?
        for (G4int k = 0; k < m_numberBuiltPMTs; ++k)
        {
            G4double cz = (m_singleHemisphere || k < m_totalNumberPMTs / 2) ? m_cylinderHeight : -m_cylinderHeight;
            G4ThreeVector center(0, 0, cz);
            G4ThreeVector axis = (m_positionsGelpad[k] - center).unit();
            for (G4double offDeg : {0.0, -8.0, 8.0})
            {
                G4ThreeVector dir = axis;
                dir.rotate(offDeg * deg, axis.orthogonal());
                G4double lastGel = -1, firstGlass = -1;
                for (G4double r = 150; r <= 216; r += 0.01)
                {
                    G4ThreeVector p = center + r * mm * dir;
                    if (m_gelPadLogical[k]->GetSolid()->Inside(p) != kOutside) lastGel = r;
                    if (firstGlass < 0 && innerAirSolid->Inside(p) == kOutside) firstGlass = r;
                }
                std::cout << "DIAG gap PMT" << k << " off=" << offDeg << "deg: gel ends r=" << lastGel
                          << " glass starts r=" << firstGlass << " gap=" << firstGlass - lastGel << " mm" << std::endl;
            }
        }
    }

    auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World_log");
    if (worldLV)
    {
        G4cout << "=== POM volume tree validation ===" << G4endl;
        PrintVolumeTree(worldLV, 0);
    }

    if (auto gel = G4PhysicalVolumeStore::GetInstance()->GetVolume("GelPad_0"))
        G4cout << "GelPad_0 mother: " << gel->GetMotherLogical()->GetName() << G4endl;
    if (auto pmt = G4PhysicalVolumeStore::GetInstance()->GetVolume("PMT_0"))
        G4cout << "PMT_0 mother:    " << pmt->GetMotherLogical()->GetName() << G4endl;
    if (auto glass = G4PhysicalVolumeStore::GetInstance()->GetVolume("GlassShellPhys"))
        G4cout << "GlassShellPhys mother: " << glass->GetMotherLogical()->GetName() << G4endl;

    // Slightly red, more opaque than the shared default (m_gelpadVis) - just for POM, so the
    // gelpad is easy to pick out visually (e.g. from underneath the reflector sleeve) without
    // changing the shared color other modules (LOM16 etc.) still use. Heap-allocated (not a
    // local/stack G4VisAttributes) since the logical volume only stores a pointer to it.
    G4VisAttributes *gelpadVisPOM = new G4VisAttributes(G4Colour(0.8, 0.15, 0.15, 0.5));
    for (int i = 0; i <= m_numberBuiltPMTs - 1; i++)
    {
        m_gelPadLogical[i]->SetVisAttributes(gelpadVisPOM);
    }
}

// ---------------- Component functions --------------------------------------------------------------------------------

G4UnionSolid *POM_4inch::pressureVessel(const G4double pOutRad, G4String pSuffix, G4bool p_bothHemispheres)
{
    G4Tubs *cylinderSolid = new G4Tubs("Cylinder solid" + pSuffix, 0, pOutRad, m_cylinderHeight, 0, 2 * CLHEP::pi);

    // Each hemisphere's flat cut face is placed to land exactly on the cylinder's own flat end
    // face (world z=+-m_cylinderHeight) - a textbook G4 boolean-solid degeneracy (the same one
    // the TitaniumCutter code elsewhere in this file already works around, see its own comment):
    // two solids of different types sharing an exactly-coincident surface confuses the CSG
    // inside/outside classification right at that seam, opening a real navigation gap there -
    // confirmed by photon tracking showing tracks pass straight from the glass shell into stray
    // "World" material around this exact z value, then re-emerge somewhere else in the cavity.
    // Padding the z-cut by kSeamPadding makes each hemisphere extend slightly INTO the cylinder's
    // own volume instead of exactly touching it - a deliberate small overlap, not a knife-edge.
    const G4double kSeamPadding = 0.01 * mm;
    G4Ellipsoid *topHalfSphere = new G4Ellipsoid("SphereTop solid" + pSuffix, pOutRad, pOutRad, pOutRad, -kSeamPadding, pOutRad);

    // place hemispheres at the cylinder ends so they start where the cylinder ends
    G4UnionSolid *topUnion = new G4UnionSolid("temp" + pSuffix, cylinderSolid, topHalfSphere, 0, G4ThreeVector(0, 0, m_cylinderHeight));

    // In single-hemisphere mode the lower hemisphere is omitted entirely: the cylinder's own
    // flat end face at -m_cylinderHeight becomes the vessel's boundary there, matching a
    // prototype whose missing hemisphere is sealed by a flat cap rather than a glass dome.
    if (!p_bothHemispheres)
        return topUnion;

    G4Ellipsoid *bottomHalfSphere = new G4Ellipsoid("SphereBottom solid" + pSuffix, pOutRad, pOutRad, pOutRad, -pOutRad, kSeamPadding);
    G4UnionSolid *unionSolid = new G4UnionSolid("OM body" + pSuffix, topUnion, bottomHalfSphere, 0, G4ThreeVector(0, 0, -m_cylinderHeight));
    return unionSolid;
}

void POM_4inch::appendEquatorBand()
{
    G4Box *cuttingBox = new G4Box("Cutter", m_glassOutRad + 10 * mm, m_glassOutRad + 10 * mm, m_equatorialBandWidth / 2.);
    G4UnionSolid *outerSolid = pressureVessel(m_glassOutRad + 0.01 * mm + m_equatorialBandThickness, "BandOuter");

    G4IntersectionSolid *intersectionSolid = new G4IntersectionSolid("BandThicknessBody", cuttingBox, outerSolid, 0, G4ThreeVector(0, 0, 0));
    G4SubtractionSolid *equatorBandSolid = new G4SubtractionSolid("Equatorband_solid", intersectionSolid,
                                                                   m_components.at("PressureVessel_" + std::to_string(m_index)).VSolid, 0, G4ThreeVector(0, 0, 0));

    G4LogicalVolume *equatorBandLogical = new G4LogicalVolume(equatorBandSolid, m_data->getMaterial("NoOptic_Absorber"), "Equatorband_log");
    equatorBandLogical->SetVisAttributes(m_absorberVis);
    appendComponent(equatorBandSolid, equatorBandLogical, G4ThreeVector(0, 0, 0), G4RotationMatrix(), "TeraTape");
}

/*
G4SubtractionSolid *POM_4inch::substractHarnessPCA(G4VSolid *p_solid)
{
    Component plug = m_harness->getComponent("CAD_PCA");
    G4Transform3D plugTransform = G4Transform3D(plug.Rotation, plug.Position);
    G4SubtractionSolid *solidSubstracted = new G4SubtractionSolid(p_solid->GetName() + "_PCASubstracted", p_solid, plug.VSolid, plugTransform);
    return solidSubstracted;
}
// ---------------- Module specific funtions below --------------------------------------------------------------------------------
*/
void POM_4inch::InternalCADComponents(G4LogicalVolume *p_innerVolume)
{
    G4RotationMatrix lRotationInternal;
    G4RotationMatrix lRotation_frame_up;
    G4RotationMatrix lRotation_flange_up;
    G4RotationMatrix lRotation_frame_down;
    //lRotation_frame_up.rotateY(90*deg);
    //lRotation_flange_up.rotateX(90*deg);
    lRotation_frame_up.rotateX(90*deg);
    lRotation_frame_down.rotateX(-90*deg);
    //lRotation_frame_down.rotateZ(-90*deg);
    //lRotation_frame_down.rotateY(-45*deg);
    G4ThreeVector lOriginInternal(0 * mm, 0 * mm, 0 * mm);
    G4ThreeVector lOriginInternal_up, lOriginInternal_down;
    // Position of each frame's CENTRE (centre of its bounding box) in module coordinates (x, y, z = module
    // axis). The 4-inch frame mesh has its corner at the origin (it spans 0-378 x 0-211 x 0-378 mm), so
    // frameTranslation() below shifts it to rotate about, and be placed by, its own centre.
    lOriginInternal_up = G4ThreeVector(0 * mm, 0 * mm, 170 * mm);
    lOriginInternal_down = G4ThreeVector(0 * mm, 0 * mm, -170 * mm);
    // Placement translation that puts the solid's bounding-box centre, after rotation, at p_centre.
    auto frameTranslation = [](G4VSolid *p_solid, const G4RotationMatrix &p_rotation, const G4ThreeVector &p_centre)
    {
        G4ThreeVector lMin, lMax;
        p_solid->BoundingLimits(lMin, lMax);
        const G4ThreeVector lTranslation = p_centre - p_rotation * (0.5 * (lMin + lMax));
        log_info("{}: centre at ({}, {}, {}) mm, size ({:.1f}, {:.1f}, {:.1f}) mm before rotation", p_solid->GetName(),
                 p_centre.x() / mm, p_centre.y() / mm, p_centre.z() / mm, (lMax - lMin).x() / mm, (lMax - lMin).y() / mm, (lMax - lMin).z() / mm);
        return lTranslation;
    };
    //Support structure
    //Tools::AppendCADComponent(this, 1.0, lOriginInternal, lRotationInternal, "POM/SupportStructure_250213.obj", "CAD_SupportStructure", m_data->getMaterial("NoOptic_Stahl"), m_steelVis, m_data->getOpticalSurface("Surf_StainlessSteelGround"));
    //{
    //auto comp = getComponent("CAD_SupportStructure");
    //new G4PVPlacement(G4Transform3D(lRotationInternal, G4ThreeVector()),
    //    comp.VLogical, "CAD_SupportStructure_physical", p_innerVolume, false, 0, m_checkOverlaps);
    //deleteComponent("CAD_SupportStructure");
    //}
    //log_info("Adding glass hemisphere and frame!");
    //Tools::AppendCADComponent(this, 1.0, lOriginInternal, lRotationInternal, "POM/GlasHemisphere.obj", "CAD_Glass", m_data->getMaterial("RiAbs_Glass_Vitrovex"), m_boardVis);
    log_info("Adding frame!");
    // TEST: POM_FRAME_SCALE shrinks/grows the upper frame about the dome centre (0,0,m_cylinderHeight),
    // keeping every frame feature at the same (theta, phi) seen from the centre. AppendCADComponent
    // scales about the mesh origin (vertex*s + offset*s, then rotated, mesh y -> module z), so the
    // offset is shifted by (1-s)/s * dome-centre height to make the dome centre the fixed point.
    //const char *frameScaleEnv = std::getenv("POM_FRAME_SCALE");
    //const G4double frameScale = frameScaleEnv ? std::atof(frameScaleEnv) : 1.0;
    //G4ThreeVector frameOffsetScaled = lOriginInternal_up;
    //if (m_singleHemisphere && frameScale != 1.0)
    //{
    //    frameOffsetScaled += G4ThreeVector(0, (1.0 - frameScale) / frameScale * m_cylinderHeight, 0);
    //    log_info("Frame scaled by {} about the dome centre", frameScale);
    //}
    Tools::AppendCADComponent(this, 1.0, G4ThreeVector(), lRotation_frame_up, "POM/p-om frame_4_inch_v1_4.obj", "CAD_Frame_up",m_data->getMaterial("NoOptic_Absorber"), m_pom_frame);
    //Tools::AppendCADComponent(this, 10.0, lOriginInternal, lRotation_flange_up, "POM/p-om_flange_glass.obj", "CAD_Flange",m_data->getMaterial("NoOptic_Absorber"), m_pom_flange);
    //Tools::AppendCADComponent(this, 10.0, lOriginInternal, lRotation_flange_up, "POM/p-om_flange_frame.obj", "CAD_glass_flange",m_data->getMaterial("Titanium"),m_pom_flange);
    {
    auto comp = getComponent("CAD_Frame_up");
    // At the CAD position the frame reaches 201.0 mm from the dome centre, just inside the glass
    // (201.9 mm). Clip it to the dome's inner sphere anyway (same center/radius the gelpad's own
    // sphere-cut uses) so it can never exceed InnerAirLV's curved boundary, whatever the placement.
    //G4ThreeVector sphereCenterInFrameLocal = lRotation_frame_up.inverse() * G4ThreeVector(0, 0, m_cylinderHeight);
    //G4VSolid *frameClipSphere = new G4Sphere("CAD_Frame_up_clip_sphere", 0, m_glassInRad, 0, 2 * CLHEP::pi, 0, CLHEP::pi);
    //G4VSolid *frameClippedSolid = new G4IntersectionSolid("CAD_Frame_up_clipped", comp.VSolid, frameClipSphere,
     //   G4Transform3D(G4RotationMatrix(), sphereCenterInFrameLocal));
    // With POM_GLASS_EXTEND_MM the dome (and the gelpads on it) is lifted relative to the frame, so
    // the frame can cut into the gelpads; notch the frame around every gelpad instead of overlapping.
    // Gelpads are placed unrotated at the cavity origin, so in frame-local coordinates they sit at
    // the inverse of the frame's own placement rotation.
    //if (std::getenv("POM_GLASS_EXTEND_MM"))
    //{
    //    for (auto gelLV : m_gelPadLogical)
    //        frameClippedSolid = new G4SubtractionSolid("CAD_Frame_up_notched", frameClippedSolid, gelLV->GetSolid(),
    //                                                   G4Transform3D(lRotation_frame_up.inverse(), G4ThreeVector()));
    //}
    G4LogicalVolume *frameClippedLV = new G4LogicalVolume(comp.VSolid, comp.VLogical->GetMaterial(), "CAD_Frame_up");
    frameClippedLV->SetVisAttributes(comp.VLogical->GetVisAttributes());
    // TEST (POM_FRAME_REFL = 0-1): PMT frame as a diffuse (Lambertian) reflector instead of black.
    if (const char *frameReflEnv = std::getenv("POM_FRAME_REFL"))
    {
        auto frameSurface = new G4OpticalSurface("PMTFrameSurface", unified, ground, dielectric_metal, 0.0);
        auto mpt = new G4MaterialPropertiesTable();
        const std::vector<G4double> energies = {1.0 * eV, 7.0 * eV};
        const G4double refl = std::atof(frameReflEnv);
        mpt->AddProperty("REFLECTIVITY", energies, {refl, refl});
        mpt->AddProperty("SPECULARLOBECONSTANT", energies, {0.0, 0.0});
        mpt->AddProperty("SPECULARSPIKECONSTANT", energies, {0.0, 0.0});
        mpt->AddProperty("BACKSCATTERCONSTANT", energies, {0.0, 0.0});
        frameSurface->SetMaterialPropertiesTable(mpt);
        new G4LogicalSkinSurface("PMTFrameSkin", frameClippedLV, frameSurface);
        log_info("PMT frame: diffuse reflector, reflectivity {}", refl);
    }
    new G4PVPlacement(G4Transform3D(lRotation_frame_up, frameTranslation(comp.VSolid, lRotation_frame_up, lOriginInternal_up)),
        frameClippedLV, "CAD_Frame_up_physical", p_innerVolume, false, 0, m_checkOverlaps);
    deleteComponent("CAD_Frame_up");
    }
    // The lower frame holds the lower-hemisphere PMTs, which don't exist in single-hemisphere
    // mode (replaced by the flat cap) - skip it there rather than placing an absorber with
    // nothing to hold.
    if (!m_singleHemisphere)
    {
        log_info("Adding flange!");
        Tools::AppendCADComponent(this, 1.0, G4ThreeVector(), lRotation_frame_down, "POM/p-om frame_4_inch_v1_4.obj", "CAD_Frame_down",m_data->getMaterial("NoOptic_Absorber"), m_pom_frame);
        auto comp = getComponent("CAD_Frame_down");
        // Same clip as CAD_Frame_up above, mirrored to the lower hemisphere's dome center.
        //G4ThreeVector sphereCenterInFrameLocalDown = lRotation_frame_down.inverse() * G4ThreeVector(0, 0, -m_cylinderHeight);
        //G4VSolid *frameClipSphereDown = new G4Sphere("CAD_Frame_down_clip_sphere", 0, m_glassInRad, 0, 2 * CLHEP::pi, 0, CLHEP::pi);
        //G4VSolid *frameClippedSolidDown = new G4IntersectionSolid("CAD_Frame_down_clipped", comp.VSolid, frameClipSphereDown,
        //    G4Transform3D(G4RotationMatrix(), sphereCenterInFrameLocalDown));
        //G4LogicalVolume *frameClippedLVDown = new G4LogicalVolume(frameClippedSolidDown, comp.VLogical->GetMaterial(), "CAD_Frame_down_clipped_Logical");
        //frameClippedLVDown->SetVisAttributes(comp.VLogical->GetVisAttributes());
        new G4PVPlacement(G4Transform3D(lRotation_frame_down, frameTranslation(comp.VSolid, lRotation_frame_down, lOriginInternal_down)),
            comp.VLogical, "CAD_Frame_down_physical", p_innerVolume, false, 0, m_checkOverlaps);
        deleteComponent("CAD_Frame_down");
      
    }
    /*Electronics
    if (!m_singleHemisphere)
    {
        log_info("Simplified LOM electronics are defined as absorber!");
        Tools::AppendCADComponent(this, 1.0, lOriginInternal, lRotationInternal, "POM/Electronics_250213.obj", "CAD_Electronics", m_data->getMaterial("NoOptic_Absorber"), m_boardVis);
        auto comp = getComponent("CAD_Electronics");
        new G4PVPlacement(G4Transform3D(lRotationInternal, G4ThreeVector()),
            comp.VLogical, "CAD_Electronics_physical", p_innerVolume, false, 0, m_checkOverlaps);
        deleteComponent("CAD_Electronics");
    } */
}


void POM_4inch::setPMTAndGelpadPositions()
{
    G4double rPMT;
    G4double rGelpad;
    G4double zOffsetPMT;
    G4RotationMatrix rot;
    G4RotationMatrix rot2;
    G4double thetaPMT;
    G4double phiPMT;

    G4double rBasePMT = m_glassInRad - m_gelThicknessFrontPMT - m_PMToffset;
    // Anchor the gelpad's outer face to the glass directly (independent of PMT position);
    // the PMT's own placement (via m_gelThicknessFrontPMT) then sets the front gel clearance.
    G4double rBaseGelpad = m_glassInRad - m_gelpad_thickness / 2.0;

    std::vector<G4double> thetaPMTlist = {m_thetaPolar, m_thetaEquatorial, 180. * deg - m_thetaEquatorial, 180. * deg - m_thetaPolar};
    // PMT/gelpad rings are placed on a sphere of radius rBasePMT/rBaseGelpad centered at
    // zOffsetPMT. This center must coincide with the dome's own optical sphere center
    // (m_cylinderHeight, same point pressureVessel() unions the hemisphere cap onto and
    // createGelpadLogicalVolumes() uses to cut each gelpad's outer face to the glass curvature)
    // - otherwise every PMT's true angular position, as seen from the dome's actual center, is
    // not where its (theta,phi) label says it is. Using m_cylinderHeight-m_frame_offset here
    // instead (an old mismatch, m_frame_offset=12mm) bent that true angle away from the nominal
    // ring angle - worse for the equatorial ring (65deg + its own extra radial offset) than the
    // polar ring (32.5deg): +3.3deg vs +2.0deg, worked out from these same constants. Confirmed
    // against a fine (1deg) scan showing a smooth several-degree saddle (not a sharp shadow edge)
    // right around the equatorial ring's nominal angle - consistent with a true-angle/optical
    // mismatch, not a discrete mechanical shadow.
    std::vector<G4double> zOffsetPMTlist = {m_cylinderHeight, m_cylinderHeight, -m_cylinderHeight, -m_cylinderHeight};
    std::vector<int> countPMTlist = {m_numberPolarPMTs, m_numberEqPMTs, m_numberEqPMTs, m_numberPolarPMTs};
    // Position-table value: polar rings (j=0,3) phase 0 -> phi=45/135/225/315; equatorial rings
    // (j=1,2) phase 0.5 -> phi=0/90/180/270. A clocking-offset sweep against real angular-scan
    // data found that a uniform +45deg rotation of the WHOLE module (not a per-ring phase swap,
    // which moves the two rings in opposite directions and breaks the frame's alignment with
    // them) best matches the real per-channel peak heights, consistent with the user's
    // confirmation that the inner/polar ring has members sitting exactly perpendicular to the
    // light source. That's now applied as a phi offset on the SCAN/SOURCE side (e.g. --phi 45
    // instead of --phi 0) rather than baked into the geometry here - simpler, and the module's
    // own construction (PMTs, frame, gelpads) stays exactly as built, with no risk of the frame
    // and PMTs falling out of alignment with each other.
    std::vector<G4double> phaseShiftList = {0, 0.5, 0.5, 0.};

    // Rings are ordered upper-polar, upper-equatorial, lower-equatorial, lower-polar; in
    // single-hemisphere mode only the two upper rings exist (the lower hemisphere is capped).
    const int numberOfRings = m_singleHemisphere ? 2 : 4;
    for (int j = 0; j < numberOfRings; j++)
    {
        thetaPMT = thetaPMTlist[j];
        zOffsetPMT = zOffsetPMTlist[j];
        rPMT = rBasePMT;
        rGelpad = rBaseGelpad;
        if (j >= 1 && j <= 2)
        {
            // TEST: POM_EQ_R_OFFSET_MM overrides the outer ring's radial offset (m_EqPMTrOffset = 0, as in the
            // CAD; 2.6 reproduces the old mDOM-derived geometry).
            const G4double eqOffset = std::getenv("POM_EQ_R_OFFSET_MM") ? std::atof(std::getenv("POM_EQ_R_OFFSET_MM")) * mm : m_EqPMTrOffset;
            rPMT += eqOffset;
            rGelpad += eqOffset;
        }

        for (int i = 0; i < countPMTlist[j]; i++)
        {
            rot = G4RotationMatrix();
            rot2 = G4RotationMatrix();
            // lPMTphi = (i + 0.5 * (j % 2)) * 360. * deg / lPMTCountList[j];
            phiPMT = (phaseShiftList[j]+i)*90.0*deg+45*deg;//(i + phaseShiftList[j]) * 360. * deg / countPMTlist[j];
            std::cout<<"Angle"<<phiPMT / deg <<" "<<thetaPMT / deg <<std::endl;
            G4double lPMTrho = rPMT * sin(thetaPMT);
            G4ThreeVector pmtCenter(lPMTrho * cos(phiPMT), lPMTrho * sin(phiPMT), rPMT * cos(thetaPMT) + zOffsetPMT);
            m_positionsPMT.push_back(pmtCenter);

            G4double lGelpadRho = rGelpad * sin(thetaPMT);
            m_positionsGelpad.push_back(G4ThreeVector(lGelpadRho * cos(phiPMT),
                                                       lGelpadRho * sin(phiPMT),
                                                       rGelpad * cos(thetaPMT) + zOffsetPMT));
            m_thetaPMT.push_back(thetaPMT);
            m_phiPMT.push_back(phiPMT);
            //m_zOffsetPMT.push_back(zOffsetPMT);
        }
    }
}
    

// ToDo:
// sin(90 +- ...) -> cos(...)
/*
void POM_4inch::setPMTAndGelpadPositions()
{
  const G4double totalLength = m_data->getValueWithUnit(m_PMTModel, "jOuterShape.jTotalLenght");
    const G4double cetreModuleToBottomPMTPolar = 70.9 * mm+m_cylinderHeight+20;//Gelpad thickness change
    const G4double centreModuleToBottomPMTEquatorial = 25.4 * mm+m_cylinderHeight+20;

    const G4double zCentreBottomPMT = totalLength - m_PMToffset; // Distance from bottom base of PMT to center of the Geant4 solid
    G4double thetaPMT, phiPMT, xPMT, yPMT, zPMT, xGelPad, yGelPad, zGelPad, rhoPMT, rhoGelpad;
    // calculate PMT and pads positions in the usual 4 for loop way
    for (int i = 0; i <= m_totalNumberPMTs - 1; i++)
    {

        // upper polar
        if (i >= 0 && i <= m_numberPolarPMTs - 1)
        {
            thetaPMT = m_thetaPolar;
            phiPMT = m_polarEquatorialPMTphiPhase + i * 90.0 * deg;

            rhoPMT = (147.7406  - 4.9) * mm; // For Position of PMT to Vessel wall (147.7406). 4.9mm is distance from PMT photocathode to vessel inner surface
            rhoGelpad = rhoPMT + 2 * m_gelPadDZ;

            zPMT = cetreModuleToBottomPMTPolar + zCentreBottomPMT * sin(90 * deg - thetaPMT);
            zGelPad = zPMT + 2 * m_gelPadDZ * sin(90 * deg - thetaPMT);
        }

        // upper equatorial
        if (i >= m_numberPolarPMTs && i <= m_numberPolarPMTs + m_numberEqPMTs - 1)
        {
            thetaPMT = m_thetaEquatorial;
            phiPMT = i * 90.0 * deg;

            rhoPMT = (120.1640 - 5) * mm; // For Position of PMT to Vessel wall (147.7406). 4.9mm is distance from PMT photocathode to vessel inner surface
            rhoGelpad = rhoPMT + 2 * m_gelPadDZ;

            zPMT = centreModuleToBottomPMTEquatorial + zCentreBottomPMT * sin(90 * deg - thetaPMT);
        }

        // lower equatorial
        if (i >= m_numberPolarPMTs + m_numberEqPMTs && i <= m_numberPolarPMTs + m_numberEqPMTs + m_numberEqPMTs - 1)
        {
            thetaPMT = 180 * deg - m_thetaEquatorial; // 118*deg; // 118.5*deg;
            phiPMT = i * 90.0 * deg;

            rhoPMT = (120.1640 - 5) * mm; // For Position of PMT to Vessel wall (147.7406). 4.9mm is distance from PMT photocathode to vessel inner surface
            rhoGelpad = rhoPMT + 2 * m_gelPadDZ;

            zPMT = (-centreModuleToBottomPMTEquatorial) - zCentreBottomPMT * sin(90 * deg - (180 * deg - thetaPMT)); // rodo to cos ...
        }

        // lower polar
        if (i >= m_totalNumberPMTs - m_numberPolarPMTs && i <= m_totalNumberPMTs - 1)
        {
            thetaPMT = 180 * deg - m_thetaPolar; // 144*deg; //152*deg;
            phiPMT = m_polarEquatorialPMTphiPhase + ((i)*90.0) * deg;

            rhoPMT = (147.7406- 4.9) * mm; // For Position of PMT to Vessel wall (147.7406). 4.9mm is distance from PMT photocathode to vessel inner surface
            rhoGelpad = rhoPMT + 2 * m_gelPadDZ;

            zPMT = (-cetreModuleToBottomPMTPolar) - zCentreBottomPMT * sin(90 * deg - (180 * deg - thetaPMT));
            zGelPad = zPMT - 2 * m_gelPadDZ * sin(90 * deg - (180 * deg - thetaPMT));
        }

        // PMTs
        xPMT = (rhoPMT * mm) * sin(thetaPMT) * cos(phiPMT);
        yPMT = (rhoPMT * mm) * sin(thetaPMT) * sin(phiPMT);

        // Gelpads
        
        
        xGelPad = rhoGelpad * sin(thetaPMT) * cos(phiPMT);
        yGelPad = rhoGelpad * sin(thetaPMT) * sin(phiPMT);

        // save positions in arrays
        m_positionsPMT.push_back(G4ThreeVector(xPMT, yPMT, zPMT));
        m_positionsGelpad.push_back(G4ThreeVector(xGelPad, yGelPad, zGelPad));

        // save angles in arrays
        m_thetaPMT.push_back(thetaPMT);
        m_phiPMT.push_back(phiPMT);
    }


     

    // Adjust PMT centers so the PMT tip (center +/- m_PMToffset along radial direction)
    // meets the gelpad surface. For each PMT, move its center to: gelpad_pos - unit*(m_PMToffset)
    /*for (size_t k = 0; k < m_positionsPMT.size(); ++k)
    {
        G4ThreeVector vec = m_positionsGelpad[k] - m_positionsPMT[k];
        if (vec.mag() > 0) {
            G4ThreeVector unit = vec.unit();
            m_positionsPMT[k] = m_positionsGelpad[k] - unit * m_PMToffset;
        }
    }
       


// Todo
// 4*m_gelPadDZ -> 2 2*m_gelPadDZ -> 1 ...  not needed if m_gelPadDZ is long enough
// tra and transformers declaration uniformely.
// rename some stuff for clarity THIS is WHERE gelpads are made 
*/
/*
void POM_4inch::createGelpadLogicalVolumes(G4VSolid *p_gelSolid)
{
    // getting the PMT solid
    G4VSolid *solidPMT = m_managerPMT->getPMTSolid();

    // Definition of helper volumes
    G4Cons *gelPadBasicSolid;
    G4EllipticalCone *gelPadBasicSolidEquatorial;
    G4IntersectionSolid *cutCone;
    G4LogicalVolume *gelPadLogical;
    G4SubtractionSolid *cutConeFinal;

    // Definition of semiaxes in (elliptical section) cone for titled gel pads
    G4double dx = std::cos(m_equatorialTiltAngle) / (2 * (1 + std::sin(m_equatorialTiltAngle) * std::tan(m_equatorialPadOpeningAngle))) * m_maxPMTRadius * 2; // semiaxis y at -ztop
    G4double ztop = 2 * m_gelPadDZ;
    G4double dy = m_maxPMTRadius;   //PMT dependant
    G4double Dy = dy + 2 * ztop * std::tan(m_equatorialPadOpeningAngle); // semiaxis x at +ztop
    G4double Dx = dx * Dy / dy;                                 //     Dx/dx=Dy/dy always ;  //semiaxis y at -ztop
    G4double xsemiaxis = (Dx - dx) / (2 * ztop);                // Best way it can be defined
    G4double ysemiaxis = (Dy - dy) / (2 * ztop);                // Best way it can be defined
    G4double zmax = (Dx + dx) / (2 * xsemiaxis);                // Best way it can be defined

    // For the placement of tilted  equatorial pads
    G4double dz = -dx * std::sin(m_equatorialTiltAngle);
    G4double dY3 = m_maxPMTRadius - (ztop * std::sin(m_equatorialTiltAngle) + dx * std::cos(m_equatorialTiltAngle));

    // create logical volume for each gelpad
    for (int k = 0; k <= m_totalNumberPMTs - 1; k++)
    {
        G4Transform3D *tra;
        G4Transform3D *tra2;

        m_converter.str("");
        m_converter2.str("");
        m_converter << "GelPad_" << k << "_solid";
        m_converter2 << "Gelpad_final" << k << "_logical";

        // polar gel pads: cone + overflow tube intersected with inner glass sphere to create spherical top
        if (k <= m_numberPolarPMTs - 1 or k >= m_totalNumberPMTs - m_numberPolarPMTs)
        {
            // parameters
            G4double gelpad_opening_angle = m_polarPadOpeningAngle; // use polar opening angle,50
            G4double gelpad_small_radius = m_maxPMTRadius; // base radius near PMT, 41
            G4double gelpad_thickness = m_gelPadDZ; // full thickness in Z,24
            G4double gelpad_large_radius = gelpad_small_radius + std::tan(gelpad_opening_angle) * gelpad_thickness;

            // overflow tube parameters
            G4double gelpad_overflow_max_radius = 250.0 * mm;
            G4double gelpad_overflow_height = 50.0 * mm;
            // place overflow so it slightly extends above cone
            G4double gelpad_overflow_offset = gelpad_thickness / 2.0 - gelpad_overflow_height / 2.0;

            // create cone
            G4VSolid *gelpad_cone = new G4Cons(m_converter.str() + "_cone",
                                               0 * mm,
                                               gelpad_small_radius,
                                               0 * mm,
                                               gelpad_large_radius,
                                               gelpad_thickness / 2.0,
                                               0,
                                               2 * CLHEP::pi);

            // overflow tube
            G4Tubs *gelpad_overflow_tubs = new G4Tubs(m_converter.str() + "_overflow",
                                                      0,
                                                      gelpad_overflow_max_radius,
                                                      gelpad_overflow_height / 2.0,
                                                      0,
                                                      2 * CLHEP::pi);

            // union cone + overflow
            G4VSolid *gelpad_cone_and_tubs = new G4UnionSolid(m_converter.str() + "_cone_and_tubs",
                                                              gelpad_cone,
                                                              gelpad_overflow_tubs,
                                                              nullptr,
                                                              G4ThreeVector(0, 0, gelpad_overflow_offset));

            // sphere to cut top to match inner glass curvature
            // compute the sphere placement in gelpad-local coordinates so that
            // after rotating and translating the gelpad into place the sphere
            // center sits at the global origin (center of the inner glass).
            G4double sphere_radius = m_glassInRad; // inner glass radius
            G4RotationMatrix rotation_for_calc;
            rotation_for_calc.rotateY(m_thetaPMT[k]);
            rotation_for_calc.rotateZ(m_phiPMT[k]);
            // compute the sphere center in gelpad-local coordinates so that
            // R * local_center + placement = 0  => local_center = R^{-1} * (-placement)
            G4ThreeVector local_center = rotation_for_calc.inverse() * ( - m_positionsGelpad[k] );

            G4VSolid *gelpad_sphere_cut = new G4Sphere(m_converter.str() + "_sphere_cut",
                                                       0,
                                                       sphere_radius,
                                                       0,
                                                       2 * CLHEP::pi,
                                                       0,
                                                       CLHEP::pi);

            // intersect cone+overflow with sphere to create spherical top
            G4VSolid *gelpad_solid = nullptr;
            {
                log_info("Creating gelpad internal intersection for '{}'", m_converter.str() + "_gelpad");
                if (dynamic_cast<G4TessellatedSolid *>(gelpad_cone_and_tubs)) {
                    log_warning("operand 'gelpad_cone_and_tubs' is a tessellated solid");
                }
                if (dynamic_cast<G4TessellatedSolid *>(gelpad_sphere_cut)) {
                    log_warning("operand 'gelpad_sphere_cut' is a tessellated solid");
                }
                gelpad_solid = new G4IntersectionSolid(m_converter.str() + "_gelpad",
                                                      gelpad_cone_and_tubs,
                                                      gelpad_sphere_cut,
                                                      new G4RotationMatrix(),
                                                      G4ThreeVector(0,0,0));
            }
            
            // rotation and position of gelpad
            G4RotationMatrix *rotation = new G4RotationMatrix();
            rotation->rotateY(m_thetaPMT[k]);
            rotation->rotateZ(m_phiPMT[k]);

            tra = new G4Transform3D(*rotation, G4ThreeVector(m_positionsGelpad[k]));
            G4Transform3D transformers = G4Transform3D(*rotation, G4ThreeVector(m_positionsPMT[k]));

            // subtract PMT solid to create final gelpad volume
            log_info("Creating intersection between p_gelSolid ('{}') and gelpad_solid ('{}') for '{}'", p_gelSolid->GetName(), gelpad_solid->GetName(), m_converter.str());
            if (dynamic_cast<G4TessellatedSolid *>(p_gelSolid)) log_warning("p_gelSolid is tessellated");
            if (dynamic_cast<G4TessellatedSolid *>(gelpad_solid)) log_warning("gelpad_solid is tessellated");
            cutCone = new G4IntersectionSolid(m_converter.str(), p_gelSolid, gelpad_solid, *tra);
            log_info("Subtracting PMT solid ('{}') from intersection ('{}')", solidPMT->GetName(), cutCone->GetName());
            if (dynamic_cast<G4TessellatedSolid *>(solidPMT)) log_warning("solidPMT is tessellated");
            cutConeFinal = new G4SubtractionSolid(m_converter.str(), cutCone, solidPMT, transformers);
            
            gelPadLogical = new G4LogicalVolume(cutConeFinal, m_data->getMaterial("RiAbs_Gel_Shin-Etsu"), m_converter2.str());
        };

        // upper equatorial
        if (k >= m_numberPolarPMTs && k <= m_numberPolarPMTs + m_numberEqPMTs - 1)
        {
            // the tilted pad is a cone with elliptical section. Moreover, there is an union with a disc called "Tube". This union does not affect the geometry, but it is used to shift the center of the geant solid so it can be placed in the same way as the PMT. In addition to that, Cut_tube is used to do a better substraction of the residual part that lies below the photocatode.
            gelPadBasicSolidEquatorial = new G4EllipticalCone("cone", xsemiaxis, ysemiaxis, zmax, ztop);

            // rotation and position of tilted gelpad
            G4RotationMatrix *rotation = new G4RotationMatrix();
            rotation->rotateY((180 * deg + m_equatorialTiltAngle));
            tra2 = new G4Transform3D(*rotation, G4ThreeVector(-dY3, 0, ztop * std::cos(m_equatorialTiltAngle) + dz));

            // place tilted cone
            G4Tubs *Tube = new G4Tubs("tube", 0, m_maxPMTRadius, 0.1 * mm, 0., 2 * CLHEP::pi);
            G4UnionSolid *Tilted_cone = new G4UnionSolid("union", Tube, gelPadBasicSolidEquatorial, *tra2);

            // For subtraction of gelpad reaching below the photocathode
            rotation = new G4RotationMatrix();
            tra = new G4Transform3D(*rotation, G4ThreeVector(0, 0, -30 * mm)); // 30 -> thickness/width of subtraction box ... just has to be large eough to cut overlapCLHEPng gelpads under photocathode

            // rotation and position of PMT
            G4RotationMatrix *rot = new G4RotationMatrix();
            rot->rotateY(m_thetaPMT[k]);
            rot->rotateZ(m_phiPMT[k]);
            G4Transform3D transformers = G4Transform3D(*rot, G4ThreeVector(m_positionsPMT[k]));

            // creating volumes ... basic cone, tilted cone, subtract PMT, logical volume of gelpad
            G4Tubs *Cut_tube = new G4Tubs("tub", 0, 2 * m_maxPMTRadius, 30 * mm, 0, 2 * CLHEP::pi);
            G4SubtractionSolid *Tilted_final = new G4SubtractionSolid("cut", Tilted_cone, Cut_tube, *tra);
            log_info("Creating intersection between p_gelSolid ('{}') and Tilted_final ('{}') for '{}'", p_gelSolid->GetName(), Tilted_final->GetName(), m_converter.str());
            if (dynamic_cast<G4TessellatedSolid *>(p_gelSolid)) log_warning("p_gelSolid is tessellated");
            if (dynamic_cast<G4TessellatedSolid *>(Tilted_final)) log_warning("Tilted_final is tessellated");
            cutCone = new G4IntersectionSolid(m_converter.str(), p_gelSolid, Tilted_final, transformers);
            log_info("Subtracting PMT solid ('{}') from intersection ('{}')", solidPMT->GetName(), cutCone->GetName());
            if (dynamic_cast<G4TessellatedSolid *>(solidPMT)) log_warning("solidPMT is tessellated");
            cutConeFinal = new G4SubtractionSolid(m_converter.str(), cutCone, solidPMT, transformers);
            gelPadLogical = new G4LogicalVolume(cutConeFinal, m_data->getMaterial("RiAbs_Gel_Shin-Etsu"), m_converter2.str());
        };

        // lower equatorial
        if (k >= m_numberPolarPMTs + m_numberEqPMTs && k <= m_totalNumberPMTs - m_numberPolarPMTs - 1)
        {
            gelPadBasicSolidEquatorial = new G4EllipticalCone("cone", xsemiaxis, ysemiaxis, zmax, ztop);

            // rotation and position of tilted gelpad
            G4RotationMatrix *rotation = new G4RotationMatrix();
            rotation->rotateY((180 * deg - m_equatorialTiltAngle));
            tra2 = new G4Transform3D(*rotation, G4ThreeVector(dY3, 0, ztop * std::cos(m_equatorialTiltAngle) + dz));

            // place tilted cone
            G4Tubs *Tube = new G4Tubs("tube", 0, m_maxPMTRadius, 0.1 * mm, 0., 2 * CLHEP::pi);
            G4UnionSolid *Tilted_cone = new G4UnionSolid("union", Tube, gelPadBasicSolidEquatorial, *tra2);

            // For subtraction of gelpad reaching below the photocathode
            rotation = new G4RotationMatrix();
            tra = new G4Transform3D(*rotation, G4ThreeVector(0, 0, -30 * mm));

            // rotation and position of PMT
            G4RotationMatrix *rot = new G4RotationMatrix();
            rot->rotateY(m_thetaPMT[k]);
            rot->rotateZ(m_phiPMT[k]);
            G4Transform3D transformers = G4Transform3D(*rot, G4ThreeVector(m_positionsPMT[k]));

            // creating volumes ... basic cone, tilted cone, subtract PMT, logical volume of gelpad
            G4Tubs *Cut_tube = new G4Tubs("tub", 0, 2 * m_maxPMTRadius, 30 * mm, 0, 2 * CLHEP::pi);
            G4SubtractionSolid *Tilted_final = new G4SubtractionSolid("cut", Tilted_cone, Cut_tube, *tra);
            log_info("Creating intersection between p_gelSolid ('{}') and Tilted_final ('{}') for '{}'", p_gelSolid->GetName(), Tilted_final->GetName(), m_converter.str());
            if (dynamic_cast<G4TessellatedSolid *>(p_gelSolid)) log_warning("p_gelSolid is tessellated");
            if (dynamic_cast<G4TessellatedSolid *>(Tilted_final)) log_warning("Tilted_final is tessellated");
            cutCone = new G4IntersectionSolid(m_converter.str(), p_gelSolid, Tilted_final, transformers);
            log_info("Subtracting PMT solid ('{}') from intersection ('{}')", solidPMT->GetName(), cutCone->GetName());
            if (dynamic_cast<G4TessellatedSolid *>(solidPMT)) log_warning("solidPMT is tessellated");
            cutConeFinal = new G4SubtractionSolid(m_converter.str(), cutCone, solidPMT, transformers);
            gelPadLogical = new G4LogicalVolume(cutConeFinal, m_data->getMaterial("RiAbs_Gel_Shin-Etsu"), m_converter2.str());
        };

        // save logicalvolume of gelpads in array
        m_gelPadLogical.push_back(gelPadLogical); //
    }
}
*/
    void POM_4inch::createGelpadLogicalVolumes(G4VSolid *p_gelSolid)
{   
    
    // getting the PMT solid
    G4VSolid *solidPMT = m_managerPMT->getPMTSolid();
    G4VSolid *cutCone;
    G4LogicalVolume *gelPadLogical;

    // create logical volume for each gelpad
    for (int k = 0; k <= m_numberBuiltPMTs - 1; k++)
    {
        G4Transform3D *tra;

        m_converter.str("");
        m_converter2.str("");
        m_converter << "GelPad_" << k << "_solid";
        m_converter2 << "Gelpad_final_" << k << "_logical";

        // TEST: POM_EQ_PAD_SMALL_RADIUS_MM / POM_EQ_PAD_OPENING_DEG give the outer (equatorial)
        // ring its own pad size; by default every pad uses the polar values.
        G4double padSmallRadius = m_gelpad_small_radius;
        G4double padLargeRadius = m_gelpad_large_radius;
        const bool equatorialPad = std::fabs(std::cos(m_thetaPMT[k])) < std::cos(50 * deg);
        if (equatorialPad && (std::getenv("POM_EQ_PAD_SMALL_RADIUS_MM") || std::getenv("POM_EQ_PAD_OPENING_DEG")))
        {
            const char *rEnv = std::getenv("POM_EQ_PAD_SMALL_RADIUS_MM");
            const char *aEnv = std::getenv("POM_EQ_PAD_OPENING_DEG");
            padSmallRadius = rEnv ? std::atof(rEnv) * mm : m_gelpad_small_radius;
            const G4double opening = aEnv ? std::atof(aEnv) * deg : m_polarPadOpeningAngle;
            padLargeRadius = padSmallRadius + std::tan(opening) * m_gelpad_thickness;
        }
        // The cone runs from the pad back (local z = -T/2, padSmallRadius) to kConeOvershoot past the glass, where
        // the dome sphere below cuts it: its flat top would otherwise touch that sphere exactly at the pad axis
        // (rGelpad = m_glassInRad - T/2), a tangency that broke the drawn polyhedron (see below).
        const G4double kConeOvershoot = 1.0 * mm;
        const G4double halfT = m_gelpad_thickness / 2.0;
        const G4double wallSlope = (padLargeRadius - padSmallRadius) / m_gelpad_thickness;
        const G4double coneZ[2] = {-halfT, halfT + kConeOvershoot};
        const G4double coneRin[2] = {0, 0};
        const G4double coneRout[2] = {padSmallRadius, padLargeRadius + wallSlope * kConeOvershoot};
        G4VSolid *gelpadCone = new G4Polycone(m_converter.str() + "_cone", 0, 2 * CLHEP::pi, 2, coneZ, coneRin, coneRout);
        // TEST (POM_PAD_FLARE_H_MM, POM_PAD_FLARE_DR_MM): the pad wall flares out near the glass -
        // over the last FLARE_H mm before the face the wall opens up so the face radius is FLARE_DR
        // larger than the straight cone's (two-segment approximation of the curved lip).
        if (const char *flareHEnv = std::getenv("POM_PAD_FLARE_H_MM"))
        {
            const G4double hf = std::atof(flareHEnv) * mm;
            const char *flareDrEnv = std::getenv("POM_PAD_FLARE_DR_MM");
            const G4double dr = (flareDrEnv ? std::atof(flareDrEnv) : 5.0) * mm;
            const G4double rKnee = padSmallRadius + (padLargeRadius - padSmallRadius) * (m_gelpad_thickness - hf) / m_gelpad_thickness;
            const G4double flareSlope = (padLargeRadius + dr - rKnee) / hf;
            const G4double zPlanes[3] = {-halfT, halfT - hf, halfT + kConeOvershoot};
            const G4double rInner[3] = {0, 0, 0};
            const G4double rOuter[3] = {padSmallRadius, rKnee, padLargeRadius + dr + flareSlope * kConeOvershoot};
            gelpadCone = new G4Polycone(m_converter.str() + "_flared", 0, 2 * CLHEP::pi, 3, zPlanes, rInner, rOuter);
            if (k == 0)
                log_info("Gelpad flare: last {} mm open up by {} mm at the face (straight-cone face radius {:.1f} mm)", hf / mm, dr / mm, padLargeRadius / mm);
        }

        // No overflow ring: tested directly against real angular-scan data (optimal-assignment
        // match, mean correlation across all 7 live channels) and confirmed to make the match
        // WORSE (0.86 with vs 0.93 without, at the same phi=45 clocking) - so it's left out,
        // not just removed for rendering reasons. A 250mm flat disk union here (tried and
        // reverted) was the earlier attempt, in case this needs revisiting with a different
        // size/shape.
                G4RotationMatrix *rotation = new G4RotationMatrix();
                rotation->rotateY(m_thetaPMT[k]);
                rotation->rotateZ(m_phiPMT[k]);

        // The pad face is cut by the dome's inner glass sphere itself, so the gel touches the glass directly (an
        // inward-padded cut once left a 1 mm air layer: glass->air TIR blocked all light beyond ~49deg). With
        // POM_GEL_GLASS_GAP_MM the cut sphere shrinks by the gap, leaving that air layer on purpose.
        //
        // The pad used to be built with a sphere padded 1 mm OUTSIDE the glass and then clipped to the whole cavity
        // (cylinder + two half-ellipsoids). Navigation was right, but the drawn polyhedron of that boolean chain is
        // fragile: depending on the viewer's lineSegmentsPerCircle (vis.mac uses 50) single pads - or at 72 all of
        // them - came out as the whole cavity, i.e. gel drawn over the entire module. Cutting the overshooting cone
        // with the dome sphere alone (no tangency, no cavity union) gives the same solid and draws reliably. The pads
        // lie entirely inside their dome (outer-ring rim at ~80deg from the dome pole), so the cavity clip did nothing
        // else.
        const G4double gelGlassGap = std::getenv("POM_GEL_GLASS_GAP_MM") ? std::atof(std::getenv("POM_GEL_GLASS_GAP_MM")) * mm : 0.0;
        const G4double gelCutRadius = m_gelpad_sphere_radius - gelGlassGap;
        G4VSolid *gelpadSphereCut = new G4Orb(m_converter.str() + "_sphere_cut", gelCutRadius);
        // In single-hemisphere mode every built PMT (k < m_numberBuiltPMTs) belongs to the
        // upper ring, since the lower hemisphere doesn't exist - always use the upper sphere.
        G4double hemisphereCenterZ = (m_singleHemisphere || k < m_totalNumberPMTs / 2)
                         ? m_cylinderHeight
                         : -m_cylinderHeight;
        G4ThreeVector sphereCenterWorld(0, 0, hemisphereCenterZ);
        G4ThreeVector sphereCenterLocal = rotation->inverse() * (sphereCenterWorld - m_positionsGelpad[k]);
        G4Transform3D sphereTransform(G4RotationMatrix(), sphereCenterLocal);
        G4VSolid *gelpadSolid = new G4IntersectionSolid(m_converter.str() + "_profile",
                                gelpadCone,
                                gelpadSphereCut,
                                sphereTransform);

        // TEST (POM_OVERFLOW_MM = width beyond the pad face, POM_OVERFLOW_THICK_MM = thickness,
        // default 3 mm): interface gel squeezed out around the pad onto the glass - a thin gel shell
        // on the inner glass surface, centred on the pad axis, joined to the pad (ending at the same
        // cut sphere as the pad face, so it sits exactly on the glass).
        if (const char *overflowEnv = std::getenv("POM_OVERFLOW_MM"))
        {
            const char *ovThickEnv = std::getenv("POM_OVERFLOW_THICK_MM");
            const G4double ovThick = (ovThickEnv ? std::atof(ovThickEnv) : 3.0) * mm;
            const G4double faceHalfAngle = std::atan(padLargeRadius / m_glassInRad);
            const G4double thetaMax = faceHalfAngle + std::atof(overflowEnv) * mm / m_glassInRad;
            G4VSolid *overflow = new G4Sphere(m_converter.str() + "_overflow", gelCutRadius - ovThick, gelCutRadius,
                                              0, 2 * CLHEP::pi, 0, thetaMax);
            gelpadSolid = new G4UnionSolid(m_converter.str() + "_profile_overflow", gelpadSolid, overflow, sphereTransform);
            if (k == 0)
                log_info("Gelpad overflow: {} mm beyond the pad face, {} mm thick (cap half-angle {:.1f} deg)", std::atof(overflowEnv), ovThick / mm, thetaMax / deg);
        }

        // TEST: gel at the PMT edge. The R15458 window is widest (r = 39.87 mm) 18.4 mm behind the PMT tip.
        //   POM_GEL_WRAP_MM = h: the gel continues h mm down the bulb past that edge, 1.5 mm thick.
        //   POM_EDGE_VOID_MM = d: an air void of radius d (torus) in the gel right at the window edge.
        {
            const G4double halfT = m_gelpad_thickness / 2.0;
            const G4double zEdge = halfT - m_gelThicknessFrontPMT - 18.4 * mm;   // local z of the window edge
            const G4double rEdge = 39.87 * mm;
            if (const char *wrapEnv = std::getenv("POM_GEL_WRAP_MM"))
            {
                const G4double h = std::atof(wrapEnv) * mm;
                G4VSolid *wrap = new G4Tubs(m_converter.str() + "_wrap", 0, rEdge + 1.5 * mm, h / 2.0, 0, 2 * CLHEP::pi);
                gelpadSolid = new G4UnionSolid(m_converter.str() + "_profile_wrap", gelpadSolid, wrap, nullptr, G4ThreeVector(0, 0, zEdge - h / 2.0));
                if (k == 0) log_info("Gel wraps {} mm down the PMT bulb past the window edge", h / mm);
            }
            if (const char *voidEnv = std::getenv("POM_EDGE_VOID_MM"))
            {
                const G4double d = std::atof(voidEnv) * mm;
                G4VSolid *edgeVoid = new G4Torus(m_converter.str() + "_edgevoid", 0, d, rEdge, 0, 2 * CLHEP::pi);
                gelpadSolid = new G4SubtractionSolid(m_converter.str() + "_profile_void", gelpadSolid, edgeVoid, nullptr, G4ThreeVector(0, 0, zEdge));
                if (k == 0) log_info("Air void of radius {} mm in the gel at the PMT window edge", d / mm);
            }
        }

        // TEST (POM_RIM_GAP_MM = width in mm, along the glass, of a ring at the pad rim): the gel has lifted
        // off the glass there - a thin air layer (POM_RIM_GAP_THICK_MM, default 0.1 mm) between glass and gel,
        // while the centre of the pad face stays in optical contact.
        if (const char *rimEnv = std::getenv("POM_RIM_GAP_MM"))
        {
            const char *rimThickEnv = std::getenv("POM_RIM_GAP_THICK_MM");
            const G4double rimThick = (rimThickEnv ? std::atof(rimThickEnv) : 0.1) * mm;
            // pad rim = where the cone wall meets the glass sphere (local z = pad axis, dome centre at
            // z = -(m_glassInRad - thickness/2), cone radius padSmallRadius at z = -thickness/2)
            const G4double halfT = m_gelpad_thickness / 2.0, zDome = -(m_glassInRad - halfT);
            const G4double slope = (padLargeRadius - padSmallRadius) / m_gelpad_thickness;
            G4double zLo = -halfT, zHi = halfT + 10 * mm;
            for (int it = 0; it < 60; it++)
            {
                const G4double z = 0.5 * (zLo + zHi), r = padSmallRadius + slope * (z + halfT);
                (r * r + (z - zDome) * (z - zDome) > m_glassInRad * m_glassInRad ? zHi : zLo) = z;
            }
            const G4double rimR = padSmallRadius + slope * (zLo + halfT);
            const G4double faceHalfAngle = std::atan2(rimR, zLo - zDome);
            const G4double thetaContact = std::max(0.0, faceHalfAngle - std::atof(rimEnv) * mm / m_glassInRad);
            G4VSolid *rimShell = new G4Sphere(m_converter.str() + "_rim_gap", m_glassInRad - rimThick, m_glassInRad + 1 * mm,
                                              0, 2 * CLHEP::pi, thetaContact, CLHEP::halfpi - thetaContact);
            gelpadSolid = new G4SubtractionSolid(m_converter.str() + "_profile_rimgap", gelpadSolid, rimShell, sphereTransform);
            if (k == 0)
                log_info("Gelpad rim not in contact with the glass: {} mm wide ring (contact within {:.1f} deg of the pad axis, face edge at {:.1f} deg), {} mm air",
                         std::atof(rimEnv), thetaContact / deg, faceHalfAngle / deg, rimThick / mm);
        }

        tra = new G4Transform3D(*rotation, G4ThreeVector(m_positionsGelpad[k]));
        G4Transform3D transformers = G4Transform3D(*rotation, G4ThreeVector(m_positionsPMT[k]));
        G4Transform3D gelpadTransform = G4Transform3D(*rotation, G4ThreeVector(m_positionsGelpad[k]));

        // The local profile, placed at the pad position in module coordinates (the gelpad is placed at the cavity origin).
        cutCone = new G4DisplacedSolid(m_converter.str(), gelpadSolid, gelpadTransform);
        // Gel as MOTHER of the PMT (default), not gel minus PMT: with the subtraction, photons
        // crossing the gelpad's side wall at grazing angles (~74-80deg w.r.t. the PMT axis) were
        // navigated straight from the gel into the PMT's internal vacuum volumes, skipping the
        // PMT glass entirely - no overlap involved, a G4 navigation failure on the boolean
        // gel-minus-PMT surface shared with the PMT. Reproduced and confirmed in the single-PMT
        // bench rig, where it suppressed the measured gain just past the pad edge. The gel
        // envelope now also covers the PMT (+ HA coating) and the PMT is placed inside it
        // (placePMTs()), so the gel/PMT boundary is a plain mother/daughter boundary.
        // POM_GEL_SUBTRACT=1 restores the old subtraction for comparison.
        G4VSolid *gelpadFinal;
        if (m_gelMotherOfPMT)
        {
            G4VSolid *envelope = new G4UnionSolid(m_converter.str() + "_plus_pmt", cutCone, solidPMT, transformers);
            if (!std::getenv("POM_HA_COATING"))
                gelpadFinal = envelope;
            else
            {
                auto coating = m_managerPMT->getComponent("HACoating");
                gelpadFinal = new G4UnionSolid(m_converter.str() + "_envelope", envelope, coating.VSolid,
                                               transformers * G4Transform3D(coating.Rotation, coating.Position));
            }
        }
        else
            gelpadFinal = new G4SubtractionSolid(m_converter.str() + "_minus_pmt", cutCone, solidPMT, transformers);
        gelPadLogical = new G4LogicalVolume(gelpadFinal, m_data->getMaterial("argGel"), m_converter2.str());

        // save logical volume of gelpads in array
        m_gelPadLogical.push_back(gelPadLogical);
    }
        
}
    

void POM_4inch::placePMTs(G4LogicalVolume *p_innerVolume)
{
    for (int k = 0; k <= m_numberBuiltPMTs - 1; k++)
    {
        m_converter.str("");
        m_converter << "_" << k;

        G4RotationMatrix *rot = new G4RotationMatrix();
        rot->rotateY(m_thetaPMT[k]);
        rot->rotateZ(m_phiPMT[k]);
        G4Transform3D transformers = G4Transform3D(*rot, G4ThreeVector(m_positionsPMT[k]));

        // With the gel as mother, the PMT sits inside its own gelpad (same coordinates, since
        // the gelpad is placed at the origin of p_innerVolume) - see createGelpadLogicalVolumes().
        G4LogicalVolume *mother = m_gelMotherOfPMT ? m_gelPadLogical[k] : p_innerVolume;
        m_managerPMT->placeIt(transformers, mother, m_converter.str());
    }
}

void POM_4inch::placeGelpads(G4LogicalVolume *p_innerVolume)
{
    for (int k = 0; k <= m_numberBuiltPMTs - 1; k++)
    {
        m_converter.str("");
        m_converter << "GelPad_" << k;

        G4RotationMatrix *rot = new G4RotationMatrix();
        rot->rotateY(m_thetaPMT[k]);
        rot->rotateZ(m_phiPMT[k]);
        G4Transform3D transformers = G4Transform3D(*rot, G4ThreeVector(m_positionsPMT[k]));

        G4VPhysicalVolume *gelPhys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), m_gelPadLogical[k], m_converter.str(), p_innerVolume, false, 0, m_checkOverlaps);

        // TEST (POM_LAYER_MAT = existing material name, with POM_GEL_GLASS_GAP_MM = layer thickness): the gap
        // between pad face and glass is filled by a coupling layer of that material over the pad face
        // (e.g. RiAbs_Glass_Vitrovex for n ~ 1.48, water for n = 1.33).
        if (std::getenv("POM_LAYER_MAT") && std::getenv("POM_GEL_GLASS_GAP_MM"))
        {
            G4Material *layerMat = m_data->getMaterial(std::getenv("POM_LAYER_MAT"));
            if (k == 0) log_info("Coupling layer between pads and glass: {}, {} mm thick", std::getenv("POM_LAYER_MAT"), std::atof(std::getenv("POM_GEL_GLASS_GAP_MM")));
            const G4double gap = std::atof(std::getenv("POM_GEL_GLASS_GAP_MM")) * mm;
            // pad rim angle seen from the dome centre (same construction as the rim-gap option)
            const G4double halfT = m_gelpad_thickness / 2.0, zDome = -(m_glassInRad - halfT);
            const G4double slope = (m_gelpad_large_radius - m_gelpad_small_radius) / m_gelpad_thickness;
            G4double zLo = -halfT, zHi = halfT + 10 * mm;
            for (int it = 0; it < 60; it++)
            {
                const G4double z = 0.5 * (zLo + zHi), r = m_gelpad_small_radius + slope * (z + halfT);
                (r * r + (z - zDome) * (z - zDome) > m_glassInRad * m_glassInRad ? zHi : zLo) = z;
            }
            const G4double faceHalfAngle = std::atan2(m_gelpad_small_radius + slope * (zLo + halfT), zLo - zDome);
            G4Sphere *layer = new G4Sphere(m_converter.str() + "_layer", m_glassInRad - gap, m_glassInRad, 0, 2 * CLHEP::pi, 0, faceHalfAngle);
            // exactly from the gel face to the glass: any leftover micro-gap would act as an air layer (TIR)
            G4LogicalVolume *layerLV = new G4LogicalVolume(layer, layerMat, m_converter.str() + "_layerLV");
            new G4PVPlacement(G4Transform3D(*rot, G4ThreeVector(0, 0, m_cylinderHeight)), layerLV, m_converter.str() + "_layer", p_innerVolume, false, 0, m_checkOverlaps);
        }

        // TEST (only when POM_WALL_SIGMA is set): rough gel<->air side walls. POM_WALL_SIGMA is the
        // microfacet roughness sigma_alpha in radians; POM_WALL_DIFFUSE (0-1, default 0) is the
        // fraction of reflections that are Lambertian instead of specular about the tilted facet.
        if (const char *sigmaEnv = std::getenv("POM_WALL_SIGMA"))
        {
            static G4OpticalSurface *wallSurface = nullptr;
            if (!wallSurface)
            {
                const char *diffuseEnv = std::getenv("POM_WALL_DIFFUSE");
                const G4double diffuse = diffuseEnv ? std::atof(diffuseEnv) : 0.0;
                wallSurface = new G4OpticalSurface("GelpadWallRough", unified, ground, dielectric_dielectric, std::atof(sigmaEnv));
                auto mpt = new G4MaterialPropertiesTable();
                const std::vector<G4double> energies = {1.0 * eV, 7.0 * eV};
                mpt->AddProperty("SPECULARLOBECONSTANT", energies, {1.0 - diffuse, 1.0 - diffuse});
                mpt->AddProperty("SPECULARSPIKECONSTANT", energies, {0.0, 0.0});
                mpt->AddProperty("BACKSCATTERCONSTANT", energies, {0.0, 0.0});
                wallSurface->SetMaterialPropertiesTable(mpt);
                log_info("Gelpad walls: sigma_alpha={} rad, diffuse fraction={}", std::atof(sigmaEnv), diffuse);
            }
            if (auto airPhys = G4PhysicalVolumeStore::GetInstance()->GetVolume("InnerAirPhys"))
            {
                new G4LogicalBorderSurface(m_converter.str() + "_rough_out", gelPhys, airPhys, wallSurface);
                new G4LogicalBorderSurface(m_converter.str() + "_rough_in", airPhys, gelPhys, wallSurface);
            }
        }
        // TEST (POM_WALL_LOSS = absorbed fraction per hit, 0-1): a lossy gel<->air pad surface (polished, but a
        // fraction of the photons reaching it is absorbed - e.g. residue or damage on the cast wall). Affects only
        // light that hits the wall (or the back annulus), i.e. mostly light reaching the PMT via the wall.
        else if (const char *lossEnv = std::getenv("POM_WALL_LOSS"))
        {
            static G4OpticalSurface *lossSurface = nullptr;
            if (!lossSurface)
            {
                lossSurface = new G4OpticalSurface("GelpadWallLossy", unified, polished, dielectric_dielectric);
                auto mpt = new G4MaterialPropertiesTable();
                const std::vector<G4double> energies = {1.0 * eV, 7.0 * eV};
                const G4double keep = 1.0 - std::atof(lossEnv);
                mpt->AddProperty("REFLECTIVITY", energies, {keep, keep});
                lossSurface->SetMaterialPropertiesTable(mpt);
                log_info("Gelpad walls: {} of the photons reaching them are absorbed", std::atof(lossEnv));
            }
            if (auto airPhys = G4PhysicalVolumeStore::GetInstance()->GetVolume("InnerAirPhys"))
            {
                new G4LogicalBorderSurface(m_converter.str() + "_loss_out", gelPhys, airPhys, lossSurface);
                new G4LogicalBorderSurface(m_converter.str() + "_loss_in", airPhys, gelPhys, lossSurface);
            }
        }
    }
}

void POM_4inch::placeGelpadReflectors(G4LogicalVolume *p_innerVolume, G4VSolid *p_gelSolid)
{
    // Thin conical shell wrapping each gelpad's tapered side wall only (not the front/back
    // faces, which must stay open for light to enter/exit) - array-level test of whether a real
    // reflective coating there (Surf_PMTSideMirror) further enhances collection, and whether
    // packing one around every PMT introduces inter-PMT shadowing/crosstalk that a single
    // isolated PMT test cannot show.
    //
    // The raw sleeve (a plain radial offset of the gelpad's own cone) can poke past the cavity
    // boundary: the gelpad's outer surface already touches that boundary with essentially no
    // spare room (its far face sits exactly on the glass sphere by construction), so adding
    // extra radius without clipping overlaps the cavity's mother volume. Intersect it against
    // the same cavity solid (p_gelSolid) the gelpad itself is clipped against, exactly mirroring
    // createGelpadLogicalVolumes()'s own clipping step.
    auto mirrorSurface = m_data->getOpticalSurface("Surf_PMTSideMirror");
    auto reflectorMat = m_data->getMaterial("NoOptic_Reflector");

    for (int k = 0; k <= m_numberBuiltPMTs - 1; k++)
    {
        m_converter.str("");
        m_converter << "GelpadReflector_" << k;

        G4VSolid *sleeve = new G4Cons(m_converter.str() + "_cone",
                                      m_gelpad_small_radius, m_gelpad_small_radius + m_reflectorConeSheetThickness,
                                      m_gelpad_large_radius, m_gelpad_large_radius + m_reflectorConeSheetThickness,
                                      m_gelpad_thickness / 2.0, 0, 2 * CLHEP::pi);

        G4RotationMatrix *rot = new G4RotationMatrix();
        rot->rotateY(m_thetaPMT[k]);
        rot->rotateZ(m_phiPMT[k]);
        G4Transform3D transform(*rot, G4ThreeVector(m_positionsGelpad[k]));

        G4VSolid *clippedSleeve = new G4IntersectionSolid(m_converter.str(), p_gelSolid, sleeve, transform);
        // The PMT bulb can be wider than the pad's back radius (36mm measured pad vs ~38-39mm bulb),
        // so the sleeve would cut into the glass - remove the PMT volume from it.
        clippedSleeve = new G4SubtractionSolid(m_converter.str() + "_minus_pmt", clippedSleeve, m_managerPMT->getPMTSolid(),
                                               G4Transform3D(*rot, G4ThreeVector(m_positionsPMT[k])));

        G4LogicalVolume *sleeveLV = new G4LogicalVolume(clippedSleeve, reflectorMat, m_converter.str() + "_log");
        sleeveLV->SetVisAttributes(m_pom_flange);

        G4VPhysicalVolume *sleevePhys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), sleeveLV, m_converter.str() + "_phys", p_innerVolume, false, 0, m_checkOverlaps);

        m_converter2.str("");
        m_converter2 << "GelPad_" << k;
        if (auto gelpadPhys = G4PhysicalVolumeStore::GetInstance()->GetVolume(m_converter2.str()))
        {
            new G4LogicalBorderSurface(m_converter.str() + "_gel_to_reflector", gelpadPhys, sleevePhys, mirrorSurface);
            new G4LogicalBorderSurface(m_converter.str() + "_reflector_to_gel", sleevePhys, gelpadPhys, mirrorSurface);
        }
    }
}
