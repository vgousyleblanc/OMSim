/** @file OMSimPOM.hh
 *  @brief Construction of POM.
 *  @ingroup common
 */

#pragma once

#include "OMSimPMTConstruction.hh"
#include "OMSimOpticalModule.hh"

class POMHarness;

class POM : public OMSimOpticalModule
{
public:
    POM(G4bool p_placeHarness = true);
    ~POM(){};
    void construction();
    void PrintVolumeTree(G4LogicalVolume* lv, G4int depth);
    double getPressureVesselWeight() {return (5.38+5.35)*kg;};
    int getNumberOfPMTs() { return m_numberBuiltPMTs;};

    G4String getName()
    {
        std::stringstream ss;
        ss << "POM/" << m_index;
        return ss.str();
    }

    // Gelpad parameters, exposed as static constants (tunable here, in one place) so other
    // constructions (e.g. the single_pmt_gelpad test rig) can reuse the exact same gelpad
    // geometry as the real P-OM, rather than duplicating these numbers elsewhere.
    static const G4double m_gelThicknessFrontPMT;
    static const G4double m_gelpad_small_radius;
    static const G4double m_gelpad_thickness;
    static const G4double m_polarPadOpeningAngle;
    static const G4double m_gelpad_large_radius;

private:
    //selection variables
    POMHarness *m_harness;
    G4bool m_placeHarness = false;
    G4bool m_harnessUnion = false; //it should be true for the first module that you build, and then false
    G4SubtractionSolid *substractHarnessPCA(G4VSolid *pSolid);

    G4UnionSolid* pressureVessel(const G4double pOutRad, G4String pSuffix, G4bool p_bothHemispheres = true);

    //POM specific functions
    void InternalCADComponents(G4LogicalVolume* lInnerVolumeLogical);
    void appendEquatorBand();
    
    //for gelpad and PMT creation
    void placePMTsAndGelpads(G4VSolid* lGelSolid, G4LogicalVolume* lGelLogical);
    void setPMTAndGelpadPositions();
    void createGelpadLogicalVolumes(G4VSolid* lGelSolid);
    void placePMTs(G4LogicalVolume* lInnerVolumeLogical);
    void placeGelpads(G4LogicalVolume* lInnerVolumeLogical);

    //vectors for positions and rotations
    std::vector<G4ThreeVector> m_positionsPMT;
    std::vector<G4ThreeVector> m_positionsGelpad;
    std::vector<G4ThreeVector> m_zOffsetPMT;
    std::vector<G4double> m_thetaPMT;
    std::vector<G4double> m_phiPMT;

    //helper variables
    std::stringstream m_converter;
    std::stringstream m_converter2;

    //logical of gelpads
    std::vector<G4LogicalVolume*> m_gelPadLogical;

    G4String m_PMTModel =  "pmt_Hamamatsu_R15458_CT"; // Change PMT model,pmt_Hamamatsu_4inch" pmt_Hamamatsu_R15458_CT


    const G4double m_xInternalCAD = 68.248*mm;
    const G4double m_yInternalCAD = 0*mm;
    const G4double m_zInternalCAD = -124.218*mm;

    const G4double m_gelPadDZ = 12.0*mm; //old 30, gelpad thickness 24
  

    const G4double m_cylinder_outer_radius=229*mm;
    const G4double m_cylinderAngle = 2.5*deg;
    const G4double m_EqPMTrOffset = 2.6*mm;//2.6
    const G4double m_EqPMTzOffset = 15*mm;//62.5

    const G4double m_cylinderHeight = 80.0*mm;
    const G4double m_glassOutRad = 215.9*mm; // Outer vessel radius including 12 mm glass thickness
    const G4double m_glassThick = 14.0*mm;
    const G4double m_glassInRad = 201.9*mm;   // Inner air cavity radius, 202
    const G4double m_gelpad_sphere_radius = m_glassInRad;
    const G4int m_numberPolarPMTs = 4;
    const G4int m_numberEqPMTs = 4;
    
    const G4double m_frame_offset=12*mm;

    const G4double m_reflectorConeSheetThickness = 0.5*mm;
    const G4double m_reflectorConeToHolder = 1.55*mm;
    const G4double m_thetaPolar = 32.5*deg; //36,32.5
    const G4double m_thetaEquatorial = 65.0*deg; //62,65
  
    const G4double m_polarEquatorialPMTphiPhase = 0*deg;//45,0
    const G4double m_equatorialTiltAngle = 15.0*deg;
    
    


    const G4double m_equatorialPadOpeningAngle = 50.0*deg;
    
    const G4int m_totalNumberPMTs = (m_numberPolarPMTs + m_numberEqPMTs) * 2;

    // If true, build only the upper hemisphere (upper polar + upper equatorial PMTs) with a
    // flat plastic cap sealing the equator, instead of the full module - for comparison
    // against a half-module prototype measurement. Read from the "single_hemisphere" arg.
    const G4bool m_singleHemisphere;
    // Number of PMTs actually built (m_totalNumberPMTs, or half of it in single-hemisphere mode).
    const G4int m_numberBuiltPMTs;

    G4double m_PMToffset;
    G4double m_maxPMTRadius;

    const G4double m_equatorialBandWidth = 35 * mm; //Total width (both halves)
    const G4double m_equatorialBandThickness = 0.5 * mm; //Thickness since its a 3D object
};
