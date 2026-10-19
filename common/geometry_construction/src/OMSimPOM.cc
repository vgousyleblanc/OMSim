/** 
 *  @todo   - Clean up magic number variables and comment their meaning
 *          - Write documentation and parse current comments into Doxygen style
 */
#include "OMSimPOM.hh"
#include "OMSimTools.hh"
#include "OMSimLogger.hh"
#include "OMSimCommandArgsTable.hh"

#include <G4Cons.hh>
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include <G4Ellipsoid.hh>
#include <G4IntersectionSolid.hh>
#include <G4UnionSolid.hh>
#include <G4Polycone.hh>
#include <G4EllipticalCone.hh>
#include <G4TessellatedSolid.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4LogicalBorderSurface.hh>
#include <fstream>
#include <cstdlib>

// Gelpad is the mother volume of its PMT (see createGelpadLogicalVolumes()); POM_GEL_SUBTRACT=1
// reverts to the old gel-minus-PMT subtraction.
static const bool m_gelMotherOfPMT = !std::getenv("POM_GEL_SUBTRACT");

// Gel thickness in front of every PMT tip. All 8 PMTs sit at the same radius: in the P-OM CAD
// (P-OM_internal_v13) every PMT tip is the same distance below the inner glass surface.
// m_EqPMTrOffset used to be 2.6mm, copied from the mDOM (whose equatorial PMTs really do sit further
// out); in the P-OM it pushed the outer ring's PMTs and gelpads 2.6mm closer to the glass than the
// polar ring's, leaving the outer PMTs with 2.6mm less front gel.
const G4double POM::m_gelThicknessFrontPMT = 6 * mm;
const G4double POM::m_gelpad_small_radius = 40.0 * mm; //40
const G4double POM::m_gelpad_thickness = 22.0 * mm; //22
const G4double POM::m_polarPadOpeningAngle = 50.0 * deg;
const G4double POM::m_gelpad_large_radius = POM::m_gelpad_small_radius + std::tan(POM::m_polarPadOpeningAngle) * POM::m_gelpad_thickness;
const G4double POM::m_gelpad_sphere_radius = 201.9 * mm; // == m_glassInRad

POM::POM(G4bool p_placeHarness) : OMSimOpticalModule(new OMSimPMTConstruction()), m_placeHarness(p_placeHarness),
    m_cylinderHeight(OMSimCommandArgsTable::getInstance().get<bool>("single_hemisphere") ? 40.0*mm : 80.0*mm),
    m_singleHemisphere(OMSimCommandArgsTable::getInstance().get<bool>("single_hemisphere")),
    m_pomReflector(OMSimCommandArgsTable::getInstance().get<bool>("pom_reflector")),
    m_numberBuiltPMTs(m_singleHemisphere ? (m_numberPolarPMTs + m_numberEqPMTs) : m_totalNumberPMTs)
{
    log_info("Constructing POM");
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

void POM::PrintVolumeTree(G4LogicalVolume* lv, G4int depth)
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

void POM::construction()
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
        G4Tubs *band = new G4Tubs("CylinderBand", m_glassInRad + kCutPadding, m_cylinder_outer_radius, m_cylinderHeight, 0, 2 * CLHEP::pi);
        G4Tubs *bandCutter = new G4Tubs("CylinderBandCutter", m_glassInRad + kCutPadding, m_cylinder_outer_radius,
                                        m_cylinderHeight + kCutPadding, 0, 2 * CLHEP::pi);
        glassShellSolid = new G4SubtractionSolid("GlassShell", glassShellSolid, bandCutter, 0, G4ThreeVector());

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
    auto airMat = m_data->getMaterial("Ri_Air");

    G4LogicalVolume *glassShellLV = new G4LogicalVolume(glassShellSolid, glassMat, "GlassShellLV");
    G4LogicalVolume *innerAirLV = new G4LogicalVolume(innerAirSolid, airMat, "InnerAirLV");

    // Place the cavity inside the glass shell. The titanium ring (full-module mode only) is NOT
    // placed here - it is wider than the glass shell (see comment above) and is instead appended
    // as its own top-level component below, so it gets placed as a sibling of the shell in the world.
    new G4PVPlacement(nullptr, G4ThreeVector(), innerAirLV, "InnerAirPhys", glassShellLV, false, 0, m_checkOverlaps);

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
    }

    // Keep this component registered so the module can be placed in the world later.
    appendComponent(glassShellSolid, glassShellLV, G4ThreeVector(0, 0, 0), G4RotationMatrix(), "PressureVessel_" + std::to_string(m_index));
    if (titaniumLV)
    {
        titaniumLV->SetVisAttributes(m_pom_flange);
        appendComponent(titaniumCylinder, titaniumLV, G4ThreeVector(0, 0, 0), G4RotationMatrix(), "TitaniumFlange_" + std::to_string(m_index));
    }
    if (pvcBandLV)
    {
        pvcBandLV->SetVisAttributes(m_pom_frame);
        appendComponent(pvcBand, pvcBandLV, G4ThreeVector(0, 0, 0), G4RotationMatrix(), "PVCFrameBand_" + std::to_string(m_index));
    }

    glassShellLV->SetVisAttributes(m_glassVis);
    innerAirLV->SetVisAttributes(m_airVis);

    // Optional geometry diagnostic, only when the POM_XSECTION environment variable is set: a
    // 1mm-grid Y=0 cross-section of the built geometry via Inside() checks, written to
    // /tmp/pom_crosssection.csv, plus the gel-to-glass gap along each PMT axis.
    if (std::getenv("POM_XSECTION"))
    {
        std::ofstream diagFile("/tmp/pom_crosssection.csv");
        diagFile << "x,z,label\n";
        G4double zMin = m_singleHemisphere ? -m_cylinderHeight - 10 : -260;
        G4double zMax = 260;
        for (G4double x = -240; x <= 240; x += 1)
        {
            for (G4double z = zMin; z <= zMax; z += 1)
            {
                G4ThreeVector testPoint(x * mm, 0, z * mm);
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
                else if (pvcBand && pvcBand->Inside(testPoint) != kOutside)
                    label = "pvc";
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

G4UnionSolid *POM::pressureVessel(const G4double pOutRad, G4String pSuffix, G4bool p_bothHemispheres)
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

void POM::appendEquatorBand()
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
G4SubtractionSolid *POM::substractHarnessPCA(G4VSolid *p_solid)
{
    Component plug = m_harness->getComponent("CAD_PCA");
    G4Transform3D plugTransform = G4Transform3D(plug.Rotation, plug.Position);
    G4SubtractionSolid *solidSubstracted = new G4SubtractionSolid(p_solid->GetName() + "_PCASubstracted", p_solid, plug.VSolid, plugTransform);
    return solidSubstracted;
}
// ---------------- Module specific funtions below --------------------------------------------------------------------------------
*/
void POM::InternalCADComponents(G4LogicalVolume *p_innerVolume)
{
    G4RotationMatrix lRotationInternal;
    G4RotationMatrix lRotation_frame_up;
    G4RotationMatrix lRotation_flange_up;
    G4RotationMatrix lRotation_frame_down;
    lRotation_frame_up.rotateZ(90*deg);
    lRotation_flange_up.rotateY(90*deg);
    lRotation_frame_up.rotateY(90*deg);
    lRotation_frame_down.rotateZ(90*deg);
    lRotation_frame_down.rotateY(-90*deg);
    G4ThreeVector lOriginInternal(0 * mm, 0 * mm, 0 * mm);
    G4ThreeVector lOriginInternal_up, lOriginInternal_down;
    if (m_singleHemisphere) {
        lOriginInternal_up = G4ThreeVector(0 * mm, 45 * mm, 0*mm); // m_cylinderHeight-m_frame_offset);
        lOriginInternal_down = G4ThreeVector(0 * mm, 45 * mm, 0*mm); // -(m_cylinderHeight-m_frame_offset));
    }
    else {
        lOriginInternal_up = G4ThreeVector(0 * mm,m_cylinderHeight+45 * mm,0 *mm);
        lOriginInternal_down = G4ThreeVector(0 * mm, -m_cylinderHeight-45 * mm,0*mm);
    }
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
    Tools::AppendCADComponent(this, 1.0, lOriginInternal_up, lRotation_frame_up, "POM/PMTFrame.obj", "CAD_Frame_up",m_data->getMaterial("NoOptic_Absorber"), m_pom_frame);
    //Tools::AppendCADComponent(this, 10.0, lOriginInternal, lRotation_flange_up, "POM/p-om_flange_glass.obj", "CAD_Flange",m_data->getMaterial("NoOptic_Absorber"), m_pom_flange);
    //Tools::AppendCADComponent(this, 10.0, lOriginInternal, lRotation_flange_up, "POM/p-om_flange_frame.obj", "CAD_glass_flange",m_data->getMaterial("Titanium"),m_pom_flange);
    {
    auto comp = getComponent("CAD_Frame_up");
    // The raw CAD mesh pokes ~3mm past the cavity's curved (glass) boundary at its outer edge -
    // a real daughter-exceeds-mother overlap (confirmed via --check_overlaps), not just a
    // cosmetic coincident-surface case like the equator plates. Clip it to the dome's true inner
    // sphere (same center/radius the gelpad's own sphere-cut uses) so it can never exceed
    // InnerAirLV's own curved boundary, regardless of mesh/placement details.
    G4ThreeVector sphereCenterInFrameLocal = lRotation_frame_up.inverse() * G4ThreeVector(0, 0, m_cylinderHeight);
    G4VSolid *frameClipSphere = new G4Sphere("CAD_Frame_up_clip_sphere", 0, m_glassInRad, 0, 2 * CLHEP::pi, 0, CLHEP::pi);
    G4VSolid *frameClippedSolid = new G4IntersectionSolid("CAD_Frame_up_clipped", comp.VSolid, frameClipSphere,
        G4Transform3D(G4RotationMatrix(), sphereCenterInFrameLocal));
    G4LogicalVolume *frameClippedLV = new G4LogicalVolume(frameClippedSolid, comp.VLogical->GetMaterial(), "CAD_Frame_up_clipped_Logical");
    frameClippedLV->SetVisAttributes(comp.VLogical->GetVisAttributes());
    new G4PVPlacement(G4Transform3D(lRotation_frame_up, G4ThreeVector()),
        frameClippedLV, "CAD_Frame_up_physical", p_innerVolume, false, 0, m_checkOverlaps);
    deleteComponent("CAD_Frame_up");
    }
    // The lower frame holds the lower-hemisphere PMTs, which don't exist in single-hemisphere
    // mode (replaced by the flat cap) - skip it there rather than placing an absorber with
    // nothing to hold.
    if (!m_singleHemisphere)
    {
        log_info("Adding flange!");
        Tools::AppendCADComponent(this, 1.0, lOriginInternal_down, lRotation_frame_down, "POM/PMTFrame.obj", "CAD_Frame_down",m_data->getMaterial("Plastic"), m_pom_frame);
        auto comp = getComponent("CAD_Frame_down");
        // Same clip as CAD_Frame_up above, mirrored to the lower hemisphere's dome center.
        G4ThreeVector sphereCenterInFrameLocalDown = lRotation_frame_down.inverse() * G4ThreeVector(0, 0, -m_cylinderHeight);
        G4VSolid *frameClipSphereDown = new G4Sphere("CAD_Frame_down_clip_sphere", 0, m_glassInRad, 0, 2 * CLHEP::pi, 0, CLHEP::pi);
        G4VSolid *frameClippedSolidDown = new G4IntersectionSolid("CAD_Frame_down_clipped", comp.VSolid, frameClipSphereDown,
            G4Transform3D(G4RotationMatrix(), sphereCenterInFrameLocalDown));
        G4LogicalVolume *frameClippedLVDown = new G4LogicalVolume(frameClippedSolidDown, comp.VLogical->GetMaterial(), "CAD_Frame_down_clipped_Logical");
        frameClippedLVDown->SetVisAttributes(comp.VLogical->GetVisAttributes());
        new G4PVPlacement(G4Transform3D(lRotation_frame_down, G4ThreeVector()),
            frameClippedLVDown, "CAD_Frame_down_physical", p_innerVolume, false, 0, m_checkOverlaps);
        deleteComponent("CAD_Frame_down");
    }
    //Electronics
    if (!m_singleHemisphere)
    {
        log_info("Simplified LOM electronics are defined as absorber!");
        Tools::AppendCADComponent(this, 1.0, lOriginInternal, lRotationInternal, "POM/Electronics_250213.obj", "CAD_Electronics", m_data->getMaterial("NoOptic_Absorber"), m_boardVis);
        auto comp = getComponent("CAD_Electronics");
        new G4PVPlacement(G4Transform3D(lRotationInternal, G4ThreeVector()),
            comp.VLogical, "CAD_Electronics_physical", p_innerVolume, false, 0, m_checkOverlaps);
        deleteComponent("CAD_Electronics");
    }
}


void POM::setPMTAndGelpadPositions()
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
            rPMT += m_EqPMTrOffset;
            rGelpad += m_EqPMTrOffset;
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
void POM::setPMTAndGelpadPositions()
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
void POM::createGelpadLogicalVolumes(G4VSolid *p_gelSolid)
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
    void POM::createGelpadLogicalVolumes(G4VSolid *p_gelSolid)
{   
    
    // getting the PMT solid
    G4VSolid *solidPMT = m_managerPMT->getPMTSolid();
    G4IntersectionSolid *cutCone;
    G4LogicalVolume *gelPadLogical;

    // create logical volume for each gelpad
    for (int k = 0; k <= m_numberBuiltPMTs - 1; k++)
    {
        G4Transform3D *tra;

        m_converter.str("");
        m_converter2.str("");
        m_converter << "GelPad_" << k << "_solid";
        m_converter2 << "Gelpad_final_" << k << "_logical";

        G4VSolid *gelpadCone = new G4Cons(m_converter.str() + "_cone",
                          0,
                  m_gelpad_small_radius,
                          0,
                  m_gelpad_large_radius,
                  m_gelpad_thickness / 2.0,
                          0,
                          2 * CLHEP::pi);

        // No overflow ring: tested directly against real angular-scan data (optimal-assignment
        // match, mean correlation across all 7 live channels) and confirmed to make the match
        // WORSE (0.86 with vs 0.93 without, at the same phi=45 clocking) - so it's left out,
        // not just removed for rendering reasons. A 250mm flat disk union here (tried and
        // reverted) was the earlier attempt, in case this needs revisiting with a different
        // size/shape.
                G4RotationMatrix *rotation = new G4RotationMatrix();
                rotation->rotateY(m_thetaPMT[k]);
                rotation->rotateZ(m_phiPMT[k]);

        // The gelpad cone's outer face lands EXACTLY at m_glassInRad by construction
        // (rGelpad = m_glassInRad - gelpad_thickness/2, cone half-height = gelpad_thickness/2),
        // so a cutting sphere of radius exactly m_gelpad_sphere_radius (== m_glassInRad) is
        // exactly tangent to it all the way around the rim - the same class of G4 boolean
        // degeneracy as the glass shell/cavity crack and the titanium cutter seam fixed
        // elsewhere in this file. Confirmed by bisection (not just theory) to be the actual
        // cause of a real OGL rendering bug: every gelpad rendered as one large, ill-defined
        // blob covering most of the dome's interior instead of its true small cone shape,
        // specifically once zOffsetPMTlist (setPMTAndGelpadPositions()) was corrected to place
        // the PMT/gelpad ring exactly on the dome's own optical center - before that fix, an
        // unrelated 12mm mismatch (m_frame_offset) happened to keep this boundary just off
        // being exactly tangent, masking the issue. A first attempt at padding by only 0.01mm
        // did NOT fix the rendering (confirmed by testing), even though it's enough to resolve
        // the pure CSG point-classification degeneracy used for navigation - G4's OGL
        // surface/polyhedron tessellation evidently uses a coarser internal tolerance for
        // "is this a real cut or a degenerate sliver." 1mm (still geometrically negligible
        // against the ~130mm gelpad) reliably clears it.
        //
        // The padding goes OUTWARD (sphere 1mm beyond the glass), not inward: an inward cut left
        // a uniform 1mm air layer between every gelpad and the glass, and glass->air total
        // internal reflection then blocked all light arriving more than ~49deg from the local
        // glass normal - a hard cutoff the real directly-glued gel doesn't have. With the sphere
        // outside the glass it no longer trims anything; the gelpad's outer face is set by the
        // clip against the cavity (p_gelSolid) below, i.e. exactly the glass's inner surface.
        const G4double kGelpadCutPadding = 1.0 * mm;
        G4VSolid *gelpadSphereCut = new G4Sphere(m_converter.str() + "_sphere_cut",
                             0,
                             m_gelpad_sphere_radius + kGelpadCutPadding,
                             0,
                             2 * CLHEP::pi,
                             0,
                             CLHEP::pi);
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

        tra = new G4Transform3D(*rotation, G4ThreeVector(m_positionsGelpad[k]));
        G4Transform3D transformers = G4Transform3D(*rotation, G4ThreeVector(m_positionsPMT[k]));
        G4Transform3D gelpadTransform = G4Transform3D(*rotation, G4ThreeVector(m_positionsGelpad[k]));

        // Apply the local profile once at the PMT position, then clip it to the inner cavity.
        cutCone = new G4IntersectionSolid(m_converter.str(), p_gelSolid, gelpadSolid, gelpadTransform);
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
            auto coating = m_managerPMT->getComponent("HACoating");
            G4VSolid *envelope = new G4UnionSolid(m_converter.str() + "_plus_pmt", cutCone, solidPMT, transformers);
            gelpadFinal = new G4UnionSolid(m_converter.str() + "_envelope", envelope, coating.VSolid,
                                           transformers * G4Transform3D(coating.Rotation, coating.Position));
        }
        else
            gelpadFinal = new G4SubtractionSolid(m_converter.str() + "_minus_pmt", cutCone, solidPMT, transformers);
        gelPadLogical = new G4LogicalVolume(gelpadFinal, m_data->getMaterial("argGel"), m_converter2.str());

        // save logical volume of gelpads in array
        m_gelPadLogical.push_back(gelPadLogical);
    }
        
}
    

void POM::placePMTs(G4LogicalVolume *p_innerVolume)
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

void POM::placeGelpads(G4LogicalVolume *p_innerVolume)
{
    for (int k = 0; k <= m_numberBuiltPMTs - 1; k++)
    {
        m_converter.str("");
        m_converter << "GelPad_" << k;

        G4RotationMatrix *rot = new G4RotationMatrix();
        rot->rotateY(m_thetaPMT[k]);
        rot->rotateZ(m_phiPMT[k]);
        G4Transform3D transformers = G4Transform3D(*rot, G4ThreeVector(m_positionsPMT[k]));

        new G4PVPlacement(0, G4ThreeVector(0, 0, 0), m_gelPadLogical[k], m_converter.str(), p_innerVolume, false, 0, m_checkOverlaps);
    }
}

void POM::placeGelpadReflectors(G4LogicalVolume *p_innerVolume, G4VSolid *p_gelSolid)
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
