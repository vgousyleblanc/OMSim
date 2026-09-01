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


POM::POM(G4bool p_placeHarness) : OMSimOpticalModule(new OMSimPMTConstruction()), m_placeHarness(p_placeHarness)
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

void POM::construction()
{
    //G4VSolid *glassSolid = pressureVessel(m_glassOutRad, "Glass");
    G4VSolid *airSolid = pressureVessel(m_glassInRad, "air"); // Fill entire vessel with gel as logical volume (not placed) for intersectionsolids with gelpads, fill with air 

    // Set positions and rotations of PMTs and gelpads
    setPMTAndGelpadPositions();
    //Materials 
    auto glassMat = m_data->getMaterial("RiAbs_Glass_Vitrovex");
    auto tiMat    = m_data->getMaterial("Titanium");
    auto airMat   = m_data->getMaterial("Ri_Air");
    G4Ellipsoid *GlassTop =
        new G4Ellipsoid(
            "GlassTop",
            m_glassOutRad,
            m_glassOutRad,
            m_glassOutRad,
            0,
            m_glassOutRad
        );
        G4LogicalVolume *GlassTopLV =
        new G4LogicalVolume(
            GlassTop,
            glassMat,
            "GlassTopLV"
        );
        G4Ellipsoid *GlassBottom =
        new G4Ellipsoid(
            "GlassBottom",
            m_glassOutRad,
            m_glassOutRad,
            m_glassOutRad,
            -m_glassOutRad,
            0
        );
        G4LogicalVolume *GlassBottomLV =new G4LogicalVolume(GlassBottom,glassMat,
            "GlassBottomLV"
        );
    G4Tubs *titaniumCylinder =
        new G4Tubs(
            "TitaniumCylinder",
            m_glassInRad,
            m_glassOutRad+5*mm,
            m_cylinderHeight,
            0,
            2 * CLHEP::pi
        );
    //G4LogicalVolume* glassHalfTopLV =
    //    new G4LogicalVolume(topHalfSphere, glassMat, "GlassTopLV");

    //G4LogicalVolume* glassBottomLV =new G4LogicalVolume(bottomHalfSphere, glassMat, "GlassBottomLV");

    G4LogicalVolume* titaniumLV = new G4LogicalVolume(titaniumCylinder, tiMat, "TitaniumCylinderLV");

    // Main volumes
    //glassSolid = substractToVolume(glassSolid, G4ThreeVector(0, 0, 0), G4RotationMatrix(), "Glass");
    //airSolid = substractToVolume(airSolid, G4ThreeVector(0, 0, 0), G4RotationMatrix(), "Gel");
    
    // Logicals
    //G4LogicalVolume *lglassLogical = new G4LogicalVolume(glassSolid, m_data->getMaterial("RiAbs_Glass_Vitrovex"), " Glass_log"); // Vessel
    G4LogicalVolume *p_innerVolume = new G4LogicalVolume(airSolid, airMat, "InnerVolume"); // Inner volume of vessel (mothervolume of all internal components)
    
    G4VSolid *vesselEnvelope =
    pressureVessel(m_glassOutRad+5.0*mm, "VesselEnvelope");

    G4LogicalVolume *vesselLV =
        new G4LogicalVolume(
            vesselEnvelope,
            airMat,
            "VesselLV"
        );

    // The envelope is only a geometrical mother.
    // Do not visualize it.
    //vesselLV->SetVisAttributes(
    //    G4VisAttributes::GetInvisible()
    //);
    
   //subtract PCA
   /*
   if (m_placeHarness)
    {
        p_innerVolume = new G4LogicalVolume(substractHarnessPCA(airSolid),
                                         m_data->getMaterial("Ri_Air"),
                                         "InnerVolume");
        lglassLogical = new G4LogicalVolume(substractHarnessPCA(glassSolid),
                                           m_data->getMaterial("RiAbs_Glass_Vitrovex"),
                                           "Glass_log");
    }
    */
    vesselLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    new G4PVPlacement(nullptr,G4ThreeVector(),p_innerVolume,"AirPhys",vesselLV,false,0,m_checkOverlaps);
    new G4PVPlacement(nullptr,G4ThreeVector(),titaniumLV,"TitaniumPhys",vesselLV,false,0,m_checkOverlaps);
    new G4PVPlacement(nullptr,G4ThreeVector(0, 0, m_cylinderHeight),GlassTopLV,"GlassTopPhys",vesselLV,false,0,m_checkOverlaps);
    new G4PVPlacement(nullptr,G4ThreeVector(0, 0, -m_cylinderHeight),GlassBottomLV,"GlassBottomPhys",vesselLV,false,0,m_checkOverlaps);
    
    createGelpadLogicalVolumes(airSolid);                                                                                       // logicalvolumes of all gelpads saved globally to be placed below

    // Placements
    //new G4PVPlacement(0, G4ThreeVector(0, 0, 0), p_innerVolume, "Gel_physical", lglassLogical, false, 0, m_checkOverlaps);

    placePMTs(p_innerVolume);
    placeGelpads(p_innerVolume);

    InternalCADComponents(p_innerVolume);
    appendComponent(vesselEnvelope,vesselLV,G4ThreeVector(0, 0, 0),G4RotationMatrix(),"PressureVessel_" + std::to_string(m_index));

    //appendComponent(glassSolid, lglassLogical, G4ThreeVector(0, 0, 0), G4RotationMatrix(), "PressureVessel_" + std::to_string(m_index));
    
    //if (m_placeHarness) appendEquatorBand();

    // ---------------- visualisation attributes --------------------------------------------------------------------------------
    //lglassLogical->SetVisAttributes(m_glassVis);
    GlassTopLV->SetVisAttributes(m_glassVis);
    GlassBottomLV->SetVisAttributes(m_glassVis);
    titaniumLV->SetVisAttributes(m_pom_flange);
    p_innerVolume->SetVisAttributes(m_airVis); 
    //vesselLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    //GlassTopLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    //GlassBottomLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    //p_innerVolume->SetVisAttributes(G4VisAttributes::GetInvisible());
    for (int i = 0; i <= m_totalNumberPMTs - 1; i++)
    {
        m_gelPadLogical[i]->SetVisAttributes(m_gelpadVis); 
    }
}

// ---------------- Component functions --------------------------------------------------------------------------------

G4UnionSolid *POM::pressureVessel(const G4double pOutRad, G4String pSuffix)
{
    G4Tubs *cylinderSolid = new G4Tubs("Cylinder solid" + pSuffix, 0, pOutRad, m_cylinderHeight, 0, 2 * CLHEP::pi);

    G4Ellipsoid *topHalfSphere = new G4Ellipsoid("SphereTop solid" + pSuffix, pOutRad, pOutRad, pOutRad, 0, pOutRad);
    G4Ellipsoid *bottomHalfSphere = new G4Ellipsoid("SphereBottom solid" + pSuffix, pOutRad, pOutRad, pOutRad, -pOutRad, 0);

    // place hemispheres at the cylinder ends so they start where the cylinder ends
    G4UnionSolid *topUnion = new G4UnionSolid("temp" + pSuffix, cylinderSolid, topHalfSphere, 0, G4ThreeVector(0, 0, m_cylinderHeight));
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
    log_info("Adding flange!");
    Tools::AppendCADComponent(this, 1.0, lOriginInternal_down, lRotation_frame_down, "POM/PMTFrame.obj", "CAD_Frame_down",m_data->getMaterial("Plastic"), m_pom_frame);
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
G4double rPMT;       // radius for PMT positioning
G4double rpad; // radius for RefCone positioning
G4double zOffsetPMT;
G4double rReflector;
G4RotationMatrix rot;
G4RotationMatrix rot2;
G4double thetaPMT;
G4double phiPMT;

    G4double rBasePMT = m_glassInRad - m_PMToffset-2+5;// Overlaps ?
    G4double reflectorpad = m_glassInRad -2;

    std::vector<G4double> thetaPMTlist = {m_thetaPolar, m_thetaEquatorial, 180. * deg - m_thetaEquatorial, 180. * deg - m_thetaPolar};
    std::vector<G4double> zOffsetPMTlist = {m_cylinderHeight-m_frame_offset, m_cylinderHeight-m_frame_offset, -m_cylinderHeight+m_frame_offset, -m_cylinderHeight+m_frame_offset};
    //std::vector<G4double> zOffsetPMTlist = {m_cylinderHeight, m_cylinderHeight - m_EqPMTzOffset, -m_cylinderHeight + m_EqPMTzOffset, -m_cylinderHeight};
    std::vector<int> countPMTlist = {m_numberPolarPMTs, m_numberEqPMTs, m_numberEqPMTs, m_numberPolarPMTs};
    std::vector<G4double> phaseShiftList = {0, 0.5, 0.5, 0.};

    for (int j = 0; j < 4; j++)
    {
        thetaPMT = thetaPMTlist[j];
        zOffsetPMT = zOffsetPMTlist[j];
        rPMT = rBasePMT;
        rReflector = reflectorpad;
        //if (j >= 1 && j <= 2)
       // {
       //     rPMT += m_EqPMTrOffset;
       //     rReflector += m_EqPMTrOffset;
       // }

        for (int i = 0; i < countPMTlist[j]; i++)
        {
            rot = G4RotationMatrix();
            rot2 = G4RotationMatrix();
            // lPMTphi = (i + 0.5 * (j % 2)) * 360. * deg / lPMTCountList[j];
            phiPMT = (phaseShiftList[j]+i)*90.0*deg+45*deg;//(i + phaseShiftList[j]) * 360. * deg / countPMTlist[j];
            std::cout<<"Angle"<<phiPMT / deg <<" "<<thetaPMT / deg <<std::endl;
            G4double lPMTrho = rPMT * sin(thetaPMT);
            m_positionsPMT.push_back(G4ThreeVector(lPMTrho * cos(phiPMT), lPMTrho * sin(phiPMT), rPMT * cos(thetaPMT) + zOffsetPMT));

            G4double lRefConeRho = rReflector * sin(thetaPMT);
            m_positionsGelpad.push_back(G4ThreeVector(lRefConeRho * cos(phiPMT), lRefConeRho * sin(phiPMT), rReflector * cos(thetaPMT) + zOffsetPMT));
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
    for (int k = 0; k <= m_totalNumberPMTs - 1; k++)
    {
        G4Transform3D *tra;
        G4Transform3D *tra2;
    G4VSolid* GelPadBasic= new G4Cons("GelPadBasic",
                                        0,               //inside radius at -pDz
                                        m_gelpad_small_radius,   //outside radius at -pDz
                                        0  * mm,               //inside radius at +pDz
                                        m_gelpad_large_radius,   //outside radius at +pDz
                                        m_gelThickness/2 ,  //half length in Z
                                        0,                     //starting angle of the segment in radians
                                        2*CLHEP::pi);   

    // Tube to model overflow of interface gel between gelpad and glass
    G4Tubs* gelpad_overflow_tubs = new G4Tubs("gelpad_overflow_tub",
                                             0,
                                             m_gelpad_overflow_max_radius,
                                             m_gelpad_overflow_height / 2,
                                             0,
                                             360 * degree);

    // gelpad cone and tubs
    G4VSolid* gelpad_cone_and_tubs = new G4UnionSolid("gelpad_cone_and_tubs",
                                                      GelPadBasic,
                                                      gelpad_overflow_tubs,
                                                      nullptr,
                                                      G4ThreeVector(0,0,0));

    // gelpad sphere models inside of glass dome (to cut edges of cone)
    // slightly larger than real counterpart
    G4VSolid* gelpad_sphere_cut = new G4Sphere("gelpad_sphere_cut",
                                               0,               //r min
                                               m_glassInRad,   //r max
                                               0,               //start phi
                                               2*CLHEP::pi,            // end phi
                                               0,               //start theta
                                               CLHEP::pi);             //end theta

    // Intersection between gelpad cone and sphere
    // creates spherical top surface on gelpad
    //G4double z_translation = - m_glassInRad+ (m_gelThickness / 2);
    
    //G4double z_translation=m_zOffsetPMT[k];
    G4VSolid* gelpad_solid = new G4IntersectionSolid("gelpad",
                                                gelpad_cone_and_tubs,                //solid 1
                                                gelpad_sphere_cut,                   //solid 2
                                                new G4RotationMatrix(0,0,0),         //rotation (identity)
                                                G4ThreeVector(0, 0, 0)); //translation
                                   
        m_converter.str("");
        m_converter2.str("");
        m_converter << "GelPad_" << k << "_solid";
        m_converter2 << "Gelpad_final" << k << "_logical";
        // polar gel pads
        // rotation and position of gelpad
        G4RotationMatrix *rotation = new G4RotationMatrix();
        rotation->rotateY(m_thetaPMT[k]);
        rotation->rotateZ(m_phiPMT[k]);
        tra = new G4Transform3D(*rotation, G4ThreeVector(m_positionsGelpad[k]));
        G4Transform3D transformers = G4Transform3D(*rotation, G4ThreeVector(m_positionsPMT[k]));
        // creating volumes ... basic cone, subtract PMT, logical volume of gelpad
        cutCone = new G4IntersectionSolid(m_converter.str(), p_gelSolid, GelPadBasic, *tra);
        cutConeFinal = new G4SubtractionSolid(m_converter.str(), cutCone, solidPMT, transformers);
        gelPadLogical = new G4LogicalVolume(cutConeFinal, m_data->getMaterial("RiAbs_Gel_Shin-Etsu"), m_converter2.str());
        // save logicalvolume of gelpads in array
        m_gelPadLogical.push_back(gelPadLogical); //
    }
        
}
    

void POM::placePMTs(G4LogicalVolume *p_innerVolume)
{
    for (int k = 0; k <= m_totalNumberPMTs - 1; k++)
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
    for (int k = 0; k <= m_totalNumberPMTs - 1; k++)
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
