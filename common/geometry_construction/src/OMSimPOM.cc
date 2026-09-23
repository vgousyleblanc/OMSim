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
#include <G4Polycone.hh>
#include <G4EllipticalCone.hh>
#include <G4TessellatedSolid.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4PhysicalVolumeStore.hh>

const G4double POM::m_gelThicknessFrontPMT = 2.0 * mm;
const G4double POM::m_gelpad_small_radius = 40.0 * mm;
const G4double POM::m_gelpad_thickness = 24.0 * mm;
const G4double POM::m_polarPadOpeningAngle = 50.0 * deg;
const G4double POM::m_gelpad_large_radius = POM::m_gelpad_small_radius + std::tan(POM::m_polarPadOpeningAngle) * POM::m_gelpad_thickness;

POM::POM(G4bool p_placeHarness) : OMSimOpticalModule(new OMSimPMTConstruction()), m_placeHarness(p_placeHarness),
    m_singleHemisphere(OMSimCommandArgsTable::getInstance().get<bool>("single_hemisphere")),
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
    // Build a real hollow shell so the optical hierarchy is:
    // world (water) -> glass shell -> inner cavity -> gelpads/PMTs.
    // This mirrors the intended POM geometry, with the shell split into two half-spheres
    // and a titanium ring rather than a single air mother enclosing glass pieces.
    G4VSolid *innerAirSolid = pressureVessel(m_glassInRad, "InnerAirVoid", !m_singleHemisphere);
    G4VSolid *outerGlassSolid = pressureVessel(m_glassOutRad, "GlassShellOuter", !m_singleHemisphere);
    G4VSolid *glassShellSolid = new G4SubtractionSolid("GlassShellMinusCavity", outerGlassSolid, innerAirSolid, 0, G4ThreeVector());

    // The single-hemisphere prototype has no titanium ring at all: the plastic cap is glued
    // directly onto the open end of the glass, so the titanium cylinder/cutter are only built
    // in full-module mode.
    G4Tubs *titaniumCylinder = nullptr;
    G4LogicalVolume *titaniumLV = nullptr;
    if (!m_singleHemisphere)
    {
        // The titanium ring is physically thicker than the glass and the glass hemispheres are
        // glued to it, so its real outer radius (m_cylinder_outer_radius) genuinely extends past
        // m_glassOutRad. That means it cannot be a daughter of the glass shell (a daughter can
        // never exceed its mother's own bounds) - it must sit alongside the shell as a sibling,
        // both placed directly in the world. The glass shell still needs this same cylindrical
        // band carved out of it, since no glass exists there at all (only titanium).
        titaniumCylinder = new G4Tubs("TitaniumCylinder", m_glassInRad, m_cylinder_outer_radius, m_cylinderHeight, 0, 2 * CLHEP::pi);
        // Cut the glass shell with a slightly padded copy of the cavity+titanium region, not the
        // exact solids themselves: their inner boundary exactly coincides with innerAirSolid's own
        // surface (both at m_glassInRad), which is a classic G4 boolean-solid degeneracy - it
        // confuses inside/outside classification right at that shared surface.
        const G4double kCutPadding = 0.01 * mm;
        G4Tubs *titaniumCutter = new G4Tubs("TitaniumCutter", m_glassInRad - kCutPadding, m_cylinder_outer_radius,
                                             m_cylinderHeight + kCutPadding, 0, 2 * CLHEP::pi);
        glassShellSolid = new G4SubtractionSolid("GlassShell", glassShellSolid, titaniumCutter, 0, G4ThreeVector());

        auto tiMat = m_data->getMaterial("Titanium");
        titaniumLV = new G4LogicalVolume(titaniumCylinder, tiMat, "TitaniumCylinderLV");
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
    InternalCADComponents(innerAirLV);

    if (m_singleHemisphere)
    {
        // Omitting the lower hemisphere and the titanium ring leaves the lower half of the
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
    if (!m_singleHemisphere)
        appendComponent(titaniumCylinder, titaniumLV, G4ThreeVector(0, 0, 0), G4RotationMatrix(), "TitaniumFlange_" + std::to_string(m_index));

    glassShellLV->SetVisAttributes(m_glassVis);
    innerAirLV->SetVisAttributes(m_airVis);
    if (!m_singleHemisphere)
        titaniumLV->SetVisAttributes(m_pom_flange);

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

    for (int i = 0; i <= m_numberBuiltPMTs - 1; i++)
    {
        m_gelPadLogical[i]->SetVisAttributes(m_gelpadVis);
    }
}

// ---------------- Component functions --------------------------------------------------------------------------------

G4UnionSolid *POM::pressureVessel(const G4double pOutRad, G4String pSuffix, G4bool p_bothHemispheres)
{
    G4Tubs *cylinderSolid = new G4Tubs("Cylinder solid" + pSuffix, 0, pOutRad, m_cylinderHeight, 0, 2 * CLHEP::pi);

    G4Ellipsoid *topHalfSphere = new G4Ellipsoid("SphereTop solid" + pSuffix, pOutRad, pOutRad, pOutRad, 0, pOutRad);

    // place hemispheres at the cylinder ends so they start where the cylinder ends
    G4UnionSolid *topUnion = new G4UnionSolid("temp" + pSuffix, cylinderSolid, topHalfSphere, 0, G4ThreeVector(0, 0, m_cylinderHeight));

    // In single-hemisphere mode the lower hemisphere is omitted entirely: the cylinder's own
    // flat end face at -m_cylinderHeight becomes the vessel's boundary there, matching a
    // prototype whose missing hemisphere is sealed by a flat cap rather than a glass dome.
    if (!p_bothHemispheres)
        return topUnion;

    G4Ellipsoid *bottomHalfSphere = new G4Ellipsoid("SphereBottom solid" + pSuffix, pOutRad, pOutRad, pOutRad, -pOutRad, 0);
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
    G4ThreeVector lOriginInternal_up(0 * mm,m_cylinderHeight-m_frame_offset, 0 * mm);
    G4ThreeVector lOriginInternal_down(0 * mm,m_cylinderHeight-m_frame_offset, 0 * mm);
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
    new G4PVPlacement(G4Transform3D(lRotation_frame_up, G4ThreeVector()),
        comp.VLogical, "CAD_Frame_up_physical", p_innerVolume, false, 0, m_checkOverlaps);
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
        new G4PVPlacement(G4Transform3D(lRotation_frame_down, G4ThreeVector()),
            comp.VLogical, "CAD_Frame_down_physical", p_innerVolume, false, 0, m_checkOverlaps);
        deleteComponent("CAD_Frame_down");
    }
    //Electronics
    
    log_info("Simplified LOM electronics are defined as absorber!");
    Tools::AppendCADComponent(this, 1.0, lOriginInternal, lRotationInternal, "POM/Electronics_250213.obj", "CAD_Electronics", m_data->getMaterial("NoOptic_Absorber"), m_boardVis);
    {
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
    std::vector<G4double> zOffsetPMTlist = {m_cylinderHeight-m_frame_offset, m_cylinderHeight-m_frame_offset, -m_cylinderHeight+m_frame_offset, -m_cylinderHeight+m_frame_offset};
    std::vector<int> countPMTlist = {m_numberPolarPMTs, m_numberEqPMTs, m_numberEqPMTs, m_numberPolarPMTs};
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
    G4SubtractionSolid *cutConeFinal;

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

                G4RotationMatrix *rotation = new G4RotationMatrix();
                rotation->rotateY(m_thetaPMT[k]);
                rotation->rotateZ(m_phiPMT[k]);

        G4VSolid *gelpadSphereCut = new G4Sphere(m_converter.str() + "_sphere_cut",
                             0,
                             m_gelpad_sphere_radius,
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
        cutConeFinal = new G4SubtractionSolid(m_converter.str() + "_minus_pmt", cutCone, solidPMT, transformers);
        gelPadLogical = new G4LogicalVolume(cutConeFinal, m_data->getMaterial("RiAbs_Gel_Shin-Etsu"), m_converter2.str());

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

        m_managerPMT->placeIt(transformers, p_innerVolume, m_converter.str());
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
