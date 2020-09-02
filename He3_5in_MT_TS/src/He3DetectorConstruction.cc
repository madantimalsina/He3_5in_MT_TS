//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// 
/// \file He3DetectorConstruction.cc
/// \brief Implementation of the He3DetectorConstruction class

#include "He3DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* He3DetectorConstruction::fMagFieldMessenger = nullptr; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

He3DetectorConstruction::He3DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fAbsorberPV(nullptr),
   fGapPV(nullptr),
   fAbsorberPV1(nullptr),
   fAbsorberPV2(nullptr),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

He3DetectorConstruction::~He3DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* He3DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void He3DetectorConstruction::DefineMaterials()
{ 
  // Define materials

  G4double A;  // atomic mass (mass of a mole)
  G4double Z;  // atomic number (mean number of protons)
  G4double d;  // density

  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  //nistManager->FindOrBuildMaterial("G4_TEFLON");
  nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
  nistManager->FindOrBuildMaterial("G4_POLYPROPYLENE"); 
  nistManager->FindOrBuildMaterial("G4_AIR");
  nistManager->FindOrBuildMaterial("G4_Al");

  // plexiglass, lucite 
  A = 1.01*g/mole;
  G4Element* elH = new G4Element ("Hydrogen","H",Z = 1.,A);
  A = 16.00*g/mole;
  G4Element* elO = new G4Element("Oxygen","O",Z = 8.,A);
  A = 12.011*g/mole;
  G4Element* elC = new G4Element("Carbon","C",Z = 6.,A);
  d = 1.19*g/cm3;
  G4Material* matplexiglass = new G4Material("Plexiglass",d,3);
  matplexiglass->AddElement(elH,0.08);
  matplexiglass->AddElement(elC,0.60);
  matplexiglass->AddElement(elO,0.32);

  //Volume of the He3 tubes
  G4double volume_He3 = (13*2.54*pi*1.18872*1.18872);

  // Build He3 gas
  G4int protons=2, neutrons=1, nucleons=protons+neutrons;
  G4double elements;
  G4double atomicMass_He3 = 3.016*g/mole; //molar mass
  G4Isotope* isoHe3 = new G4Isotope("He3", protons, nucleons, atomicMass_He3);
  G4Element* elHe3 = new G4Element("Helium3", "He3", 1);
  elHe3->AddIsotope(isoHe3, 100*perCent);
  G4double pressure_He3 = 11.02*atmosphere;
  G4double temperature = 293.15*kelvin;
  G4double molar_constant = Avogadro*k_Boltzmann;
  //G4double density = (atomicMass*pressure)/(temperature*molar_constant);
  G4double density_He3 = (atomicMass_He3*pressure_He3)/(temperature*molar_constant);

  G4Material* Helium3 = new G4Material("Helium3", density_He3, elements=1, kStateGas, temperature, pressure_He3);
  Helium3->AddElement(elHe3, 100*perCent);

  // Argon 
  G4double atomicMass_Ar = 39.948*g/mole;
  G4double pressure_Ar = 0.58*atmosphere;
  G4double density_Ar = (atomicMass_Ar*pressure_Ar)/(temperature*molar_constant);
  G4Element* elAr = new G4Element("Argon", "Ar", Z=18., atomicMass_Ar);
  G4Material* Argon = new G4Material("Argon"  , density_Ar, 1,  kStateGas, temperature, pressure_Ar);
  Argon->AddElement(elAr, 1);

/*  // Methane
  G4Material* methane  = new G4Material("methane",density, 2);
  methane->AddElement(elC,1);
  methane->AddElement(elH,4);
*/

  // 95% He3 + 4.95% Ar + 0.05% CH4,
  G4double atomicMass_He3Ar = ((0.95 *3.016) + (0.05* 39.948))*g/mole;
  G4double pressure_He3Ar = 11.60*atmosphere;
  G4double density_He3Ar = (atomicMass_He3Ar*pressure_He3Ar)/(temperature*molar_constant);

  G4Material* He3Ar = new G4Material("He3Ar"  , density_He3Ar, 2,  kStateGas, temperature, pressure_He3Ar);
  //He3Ar->AddMaterial( Helium3, 0.95 );
  //He3Ar->AddMaterial( Argon,   0.05 );
  He3Ar->AddMaterial( Helium3, (0.95 *3.016) / ((0.95 *3.016) + (0.05* 39.948)) );
  He3Ar->AddMaterial( Argon, (0.05* 39.948) / ((0.95 *3.016) + (0.05* 39.948)) );
  //He3ArCH4->AddMaterial( methane, 0.0005 );
  G4cout << "MY CHECK....!!!!  " << "pressure_He3  =  " << pressure_He3/atmosphere << G4endl;
  G4cout << "MY CHECK....!!!!  " << "pressure_Ar  =  " << pressure_Ar/atmosphere << G4endl;
  G4cout << "MY CHECK....!!!!  " << "pressure_He3Ar  =  " << pressure_He3Ar/atmosphere << G4endl;
  G4cout << G4endl;
  G4cout << "MY CHECK....!!!!  " << "density_He3  =  " << density_He3/(mg/cm3) << G4endl;
  G4cout << "MY CHECK....!!!!  " << "density_Ar  =  " << density_Ar/(mg/cm3) << G4endl;
  G4cout << "MY CHECK....!!!!  " << "density_He3Ar  =  " << density_He3Ar/(mg/cm3) << G4endl;
  G4cout << G4endl;
  G4cout << "MY CHECK....!!!!  " << "volume of the tube  =  " << volume_He3 << " cm3"<< G4endl;
  G4cout << G4endl;
  G4cout << "MY CHECK....!!!!  " << "mass_He3  =  " << (density_He3/(mg/cm3))*volume_He3 << " mg" << G4endl;
  G4cout << "MY CHECK....!!!!  " << "mass_Ar  =  " << (density_Ar/(mg/cm3))*volume_He3 << " mg"<<G4endl;
  G4cout << "MY CHECK....!!!!  " << "mass_He3Ar  =  " << (density_He3Ar/(mg/cm3))*volume_He3 << " mg"<< G4endl;
  
  
  // UHMW for Thermal Scattering of neutron
  //G4Element *C = new G4Element("C", "C", 2);
  G4Element *H = new G4Element("TS_H_of_Polyethylene", "H", 1., 1.0079*g/mole);
  G4Material *POLYMAT = new G4Material("POLYMAT", 0.94*g/cm3, 2, kStateSolid, 293.15*kelvin);
  POLYMAT->AddElement(H, 0.143711);
  POLYMAT->AddElement(elC, 0.856289);

  // UHMW (Ultra High Molecular Weight Polyethylene)
  d = 0.94*g/cm3;
  nistManager->BuildMaterialWithNewDensity("UHMWPE","G4_POLYETHYLENE",d);

  // POLYPROPYLENE
  d = 913.43685*mg/cm3; // Geant4 density: 900.000 mg/cm3
  nistManager->BuildMaterialWithNewDensity("POLYPROPYLENE","G4_POLYPROPYLENE",d);
  //G4NistManager* man = G4NistManager::Instance();
  //G4Material* matUHMW = man->BuildMaterialWithNewDensity("UHMW","G4_POLYETHYLENE",d);
  
  
/*  
  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);*/

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* He3DetectorConstruction::DefineVolumes()
{
  // Geometry parameters

  // Parameters: source holder
  //G4double innerRadius_sHold = 0.*cm;
  //G4double outerRadius_sHold = 16.*mm;
  G4double l_sHold = 4.3*cm;
  G4double w_sHold = 3.25*cm;
  G4double hz_sHold = 1.27*cm;

  G4double startAngle = 0.*deg;
  G4double spanningAngle = 360.*deg;

  // Outer layer of He3 tube
  G4double outerRadius_THe3 = 25.4/2*mm;
  G4double thickness_AlTube = (0.032*25.4)*mm;


  // Parameters: Inner layer of He3
  // Tubes filled with He3 gas
  G4double outerRadius_gHe3 = outerRadius_THe3 - thickness_AlTube;
  G4double innerRadius_gHe3 = 0.*mm;

  // This is the effective length (13") of He-3 tubes
  // divide by 2 because of G4 style (-33.02/2 to +33.02/2)
  G4double hz_gHe3 = 33.02/2*cm; 
  //G4double startAngle = 0.*deg;
  //G4double spanningAngle = 360.*deg;

  // Gap between UHMWPE and outer layer of He3 tube
  // Polypropylene U-Shaped Channel of Dimension,   
  // 2-1/2"(2.5" --> 6.35 Cm) Base,   
  // 2-1/4"(2.25" --> 5.715 Cm ) Legs X 3/16" Wall Thick
  // gapPEHe3T = Legs - Diameter of He3 Tubes
  G4double gapPEHe3T = 5.715*cm - 2*outerRadius_THe3; // Legs 2.25" (5.715 cm)
  
  // Rotation of the Tubes 
  // (both outer and inner filled with He3 gas)
  G4RotationMatrix* rotD1 = new G4RotationMatrix();
  G4RotationMatrix* rotD2 = new G4RotationMatrix();
  G4RotationMatrix* rotD3 = new G4RotationMatrix();
  rotD1->rotateZ(90.*deg);
  rotD2->rotateX(90.*deg);
  rotD3->rotateY(90.*deg);


  G4int nofLayers = 1;
  //G4double Thickness = 20.*mm;
  G4double gapThickness = 5*25.40*mm;  // Thickness of moderator (Varies: 1", 2", 2.5", 5" etc)
  G4double gapThicknessB = 50.8*mm;  // 2" Thickness behind (Underground) or below (SDSMT) He3 tubes
  G4double gapThicknessB2 = 4.7625*mm; //(0.1875)*25.4*mm polypropylene U-shaped channel for He3 tubes support
  G4double gapThicknessB2z = 52.3875*mm; //(2.0625)*25.4*mm polypropylene U-shaped channel for He3 tubes support
  //G4double gapThicknessSW = 50.8*mm; // 2" Thickness for side wall 1 and side wall 2
  G4double calorSizeY  = 30.48*cm;   // 12" 
  G4double calorSizeX  = 31.115*cm;  // 12.25"
  G4double calorSizeXY = 62.0*cm;
  //G4double calorSizeZ = 42.0*cm;

  // For side wall
  G4double sidePlateZ = (12.25 * 2.54)*cm; // 12.25"
  G4double sidePlateY = 50.80*mm; // 2" width for the 2-side wall
  G4double BsidePlateY = 25.40*mm; // 1" width for the back side wall
  G4double BsidePlateX = 50.0*cm; // Back side plate length

  auto layerThickness2 = 100.0*cm;
  auto layerThickness = 10.5*mm + gapThickness + 2*outerRadius_THe3 + gapPEHe3T;
  //auto calorThickness = nofLayers * layerThickness;
  auto calorThickness2 = nofLayers * layerThickness2;
  auto worldSizeXY = 5.0 * calorSizeX; 
  auto worldSizeZ  = 5.0 * (layerThickness2);
  
  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("G4_AIR");
  auto absorberMaterial = G4Material::GetMaterial("Helium3");
  //auto gapMaterial = G4Material::GetMaterial("G4_POLYETHYLENE");
  auto gapMaterial = G4Material::GetMaterial("POLYMAT");
  //auto gapMaterial = G4Material::GetMaterial("UHMWPE");
  auto gapMaterial2 = G4Material::GetMaterial("POLYPROPYLENE");
  auto holdMaterial = G4Material::GetMaterial("Plexiglass");
  auto gapMaterialT = G4Material::GetMaterial("G4_Al");

  
  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("He3DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
  // Volume defination:
  //     
  // World
  //
  auto worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                               
  // Calorimeter
  //  
  auto calorimeterS
    = new G4Box("Calorimeter",     // its name
                 calorSizeXY/2, calorSizeXY/2, calorThickness2/2); // its size
                         
  auto calorLV
    = new G4LogicalVolume(
                 calorimeterS,     // its solid
                 defaultMaterial,  // its material
                 "Calorimeter");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 calorLV,          // its logical volume                         
                 "Calorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                                 
  // Layer
  //
  auto layerS 
    = new G4Box("Layer",           // its name
                 calorSizeXY/2, calorSizeXY/2, layerThickness2/2.); // its size
                         
  auto layerLV
    = new G4LogicalVolume(
                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Layer");         // its name

  new G4PVReplica(
                 "Layer",          // its name
                 layerLV,          // its logical volume
                 calorLV,          // its mother
                 kZAxis,           // axis of replication
                 nofLayers,        // number of replica
                 layerThickness2);  // witdth of replica

  //
  // Source holder ( Plexiglass)
  //

  auto SBox 
    = new G4Box("Source",             // its name
                 l_sHold/2, w_sHold/2, hz_sHold/2); // its size
                         
  auto SBoxLV
    = new G4LogicalVolume(
                 SBox,             // its solid
                 holdMaterial, 
                 //defaultMaterial,     // its material
                 "Source");           // its name
                                   
  fHoldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -(gapThickness/2. + hz_sHold/2.)), // its position
                 SBoxLV,            // its logical volume                         
                 "Source",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

/*
  G4Tubs* sourceTube = 
    new G4Tubs("Source", innerRadius_sHold, outerRadius_sHold, hz_sHold, startAngle, spanningAngle);

  G4LogicalVolume* logicRing =                         
    new G4LogicalVolume(sourceTube,           //its solid
                        //holdMaterial, 
                        defaultMaterial,        //its material
                        "Source"); 

  fHoldPV
   = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -(gapThickness/2 + hz_sHold)), // its position
                 logicRing,            // its logical volume                         
                 "Source",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps */

  //                               
  // Gap
  // Here Gap material is UHMWPE as neutron moderator
  //
  auto gapS 
    = new G4Box("Gap",             // its name
                 calorSizeX/2, calorSizeY/2, gapThickness/2); // its size
                         
  auto gapLV
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "Gap");           // its name
                                   
  fGapPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0.), // its position
                 gapLV,            // its logical volume                         
                 "Gap",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //
  // Gap2 --> Below the He3 tube SDSMT setup 
              //OR
  //          Behind the He3 tubes for DD setup
  //  
  auto gapS2 
    = new G4Box("Gap2",             // its name
                 calorSizeX/2, calorSizeY/2, gapThicknessB/2); // its size
                         
  auto gapLV2
    = new G4LogicalVolume(
                 gapS2,             // its solid
                 gapMaterial,
                 //defaultMaterial,      // its material
                 "Gap2");           // its name
                                   
  fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., (103.188 + (12.7 + 4.76 + gapThicknessB/2))*mm), // 103.188 is the distance from the origen to center of the tube for 5" MT
                 //G4ThreeVector(0., 0., (gapThickness/2 + 2*outerRadius_THe3 + gapPEHe3T - 0.47625/2*cm + gapThicknessB/4.)), // 0.47625*cm -> wall thick of U -shape
                 gapLV2,            // its logical volume                         
                 "Gap2",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
//    
//(0.1875)*25.4*mm polypropylene U-shaped channel for He3 tubes support
// below the tube (make a single slab for now)   
  auto gapS2b 
    = new G4Box("Gap2b",             // its name
                 calorSizeX/2, calorSizeY/2, gapThicknessB2/2); // its size
                         
  auto gapLV2b
    = new G4LogicalVolume(
                 gapS2b,             // its solid
                 gapMaterial2,      // its material
                 "Gap2b");           // its name
                                   
  fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., (103.188 + (12.7 + gapThicknessB2/2))*mm), // 103.188 is the distance from the origen to center of the tube for 5" MT
                 //G4ThreeVector(0., 0., (gapThickness/2 + 2*outerRadius_THe3 + gapPEHe3T - 0.47625/2*cm + gapThicknessB/4.)), // 0.47625*cm -> wall thick of U -shape
                 gapLV2b,            // its logical volume                         
                 "Gap2b",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

//inner side wall of the U-channel (double thick wall)
  auto gapS2s 
    = new G4Box("Gap2s",             // its name
                 calorSizeX/2, 2*gapThicknessB2/2, gapThicknessB2z/2); // its size
                         
  auto gapLV2s
    = new G4LogicalVolume(
                 gapS2s,             // its solid
                 gapMaterial2,      // its material
                 "Gap2s");           // its name
                                   
  fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 31.7505*mm, gapThickness/2. + gapThicknessB2z/2.), // 103.188 is the distance from the origen to center of the tube for 5" MT
                 //G4ThreeVector(0., 0., (gapThickness/2 + 2*outerRadius_THe3 + gapPEHe3T - 0.47625/2*cm + gapThicknessB/4.)), // 0.47625*cm -> wall thick of U -shape
                 gapLV2s,            // its logical volume                         
                 "Gap2s",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

 // make a copy for other sides
                                   
  fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., (31.7505 + 63.5)*mm, gapThickness/2. + gapThicknessB2z/2.), // 103.188 is the distance from the origen to center of the tube for 5" MT
                 //G4ThreeVector(0., 0., (gapThickness/2 + 2*outerRadius_THe3 + gapPEHe3T - 0.47625/2*cm + gapThicknessB/4.)), // 0.47625*cm -> wall thick of U -shape
                 gapLV2s,            // its logical volume                         
                 "Gap2s",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 1,                // copy number
                 fCheckOverlaps);  // checking overlaps 

      fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., -31.75*mm, gapThickness/2. + gapThicknessB2z/2.), // 103.188 is the distance from the origen to center of the tube for 2.5" MT
                 //G4ThreeVector(0., 0., (gapThickness/2 + 2*outerRadius_THe3 + gapPEHe3T - 0.47625/2*cm + gapThicknessB/4.)), // 0.47625*cm -> wall thick of U -shape
                 gapLV2s,            // its logical volume                         
                 "Gap2s",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 2,                // copy number
                 fCheckOverlaps);  // checking overlaps 

      fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., -(31.75 + 63.5)*mm, gapThickness/2. + gapThicknessB2z/2.), // 103.188 is the distance from the origen to center of the tube for 2.5" MT
                 //G4ThreeVector(0., 0., (gapThickness/2 + 2*outerRadius_THe3 + gapPEHe3T - 0.47625/2*cm + gapThicknessB/4.)), // 0.47625*cm -> wall thick of U -shape
                 gapLV2s,            // its logical volume                         
                 "Gap2s",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 2,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //side wall of the U-channel (double thick wall)
  auto gapS2s2 
    = new G4Box("Gap2s2",             // its name
                 calorSizeX/2, gapThicknessB2/2, gapThicknessB2z/2); // its size
                         
  auto gapLV2s2
    = new G4LogicalVolume(
                 gapS2s2,             // its solid
                 gapMaterial2,      // its material
                 "Gap2s2");           // its name
                                   
  fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 152.4*mm, gapThickness/2. + gapThicknessB2z/2.), // 103.188 is the distance from the origen to center of the tube for 5" MT
                 //G4ThreeVector(0., 0., (gapThickness/2 + 2*outerRadius_THe3 + gapPEHe3T - 0.47625/2*cm + gapThicknessB/4.)), // 0.47625*cm -> wall thick of U -shape
                 gapLV2s2,            // its logical volume                         
                 "Gap2s2",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

// copy for the other side
  fGapPV2
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., -(152.4)*mm, gapThickness/2. + gapThicknessB2z/2.), // 103.188 is the distance from the origen to center of the tube for 5" MT
                 //G4ThreeVector(0., 0., (gapThickness/2 + 2*outerRadius_THe3 + gapPEHe3T - 0.47625/2*cm + gapThicknessB/4.)), // 0.47625*cm -> wall thick of U -shape
                 gapLV2s2,            // its logical volume                         
                 "Gap2s2",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 1,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //
  // For He3 tube outerlayer
  // Tube 1 (center tube) outerlayer
  // centre tube is off by 0.635cm(0.25")
  // Length of the tube (not symmetric)
  // 1.6"-1.25" = 0.35" active length extending out of moderator at NEAR END
  // 1.6"-0.95" = 0.65" active length extending out of moderator at FAR END
  // So the X-center is pushed by the 0.15" 
  // 103.188 is the distance from the origen to center of the tube for 2.5" MT  
  auto gapsT1 
    = new G4Tubs("T1", outerRadius_gHe3, outerRadius_THe3, hz_gHe3, startAngle, spanningAngle); 
                         
  auto gapLVT1
    = new G4LogicalVolume(
                 gapsT1,        // its solid
                 gapMaterialT, // its material
                 "T1");          // its name
                                   
  fGapPV
    = new G4PVPlacement(
                 rotD3,                // no rotation
                 G4ThreeVector((0.15 * 2.54)*cm, (0.5625*2.54)*cm, 103.188*mm), //// 103.188 is the distance from the origen to center of the tube for 5" MT
                 //G4ThreeVector((0.15 * 2.54)*cm, (0.5625*2.54)*cm, gapThickness/2. + outerRadius_THe3 + gapPEHe3T/2.), // its position centre is off by 0.635cm(0.25")
                 gapLVT1,       // its logical volume                         
                 "T1",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  // Tube 2 outerlayer
  auto gapsT2 
    = new G4Tubs("T2", outerRadius_gHe3, outerRadius_THe3, hz_gHe3, startAngle, spanningAngle); 
                         
  auto gapLVT2
    = new G4LogicalVolume(
                 gapsT2,        // its solid
                 gapMaterialT, // its material
                 "T2");          // its name
                                   
  fGapPV
    = new G4PVPlacement(
                 rotD3,                // no rotation
                 G4ThreeVector((0.15 * 2.54)*cm, -112.71*mm, 103.188*mm), //103.188 is the distance from the origen to center of the tube for 5" MT
                 gapLVT2,       // its logical volume                         
                 "T2",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps  
  
  // Tube 3 outerlayer
  auto gapsT3 
    = new G4Tubs("T3", outerRadius_gHe3, outerRadius_THe3, hz_gHe3, startAngle, spanningAngle); 
                         
  auto gapLVT3
    = new G4LogicalVolume(
                 gapsT3,        // its solid
                 gapMaterialT, // its material
                 "T3");          // its name
                                   
  fGapPV
    = new G4PVPlacement(
                 rotD3,                // no rotation
                 G4ThreeVector((0.15 * 2.54)*cm, 134.94*mm, 103.188*mm), // its position
                 gapLVT3,       // its logical volume                         
                 "T3",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                               
  // Absorber
  // Absorber is the He3 gas
  // For He3 gas Tube 1 (Center tube) 
  // Centre tube is off by 0.635cm (0.25")
  // Length of the tube (not symmetric)
  // 1.6"-1.25" = 0.35" active length extending out of moderator at NEAR END
  // 1.6"-0.95" = 0.65" active length extending out of moderator at FAR END
  // So the X-center is pushed by the 0.15" 
  // 52.388 is the distance from the origen to center of the tube for 1" MT  
  auto absorberS 
    = new G4Tubs("Abso", innerRadius_gHe3, outerRadius_gHe3, hz_gHe3, startAngle, spanningAngle); 
                         
  auto absorberLV
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "Abso");          // its name
                                   
  fAbsorberPV
    = new G4PVPlacement(
                 rotD3,                // no rotation
                 G4ThreeVector((0.15 * 2.54)*cm, (0.5625*2.54)*cm, 103.188*mm), // its position Centre tube is off by 0.635cm (0.25")
                 absorberLV,       // its logical volume                         
                 "Abso",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  // For He3 gas Tube 2
  auto absorberS2 
    = new G4Tubs("Abso2", innerRadius_gHe3, outerRadius_gHe3, hz_gHe3, startAngle, spanningAngle); 
                         
  auto absorberLV2
    = new G4LogicalVolume(
                 absorberS2,        // its solid
                 absorberMaterial, // its material
                 "Abso2");          // its name
                                   
  fAbsorberPV1
    = new G4PVPlacement(
                 rotD3,                // no rotation
                 G4ThreeVector((0.15 * 2.54)*cm, -(112.71)*mm, 103.188*mm), // its position
                 absorberLV2,       // its logical volume                         
                 "Abso2",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  // For He3 gas Tube 3
  auto absorberS3 
    = new G4Tubs("Abso3", innerRadius_gHe3, outerRadius_gHe3, hz_gHe3, startAngle, spanningAngle); 
                         
  auto absorberLV3
    = new G4LogicalVolume(
                 absorberS3,        // its solid
                 absorberMaterial, // its material
                 "Abso3");          // its name
                                   
  fAbsorberPV2
    = new G4PVPlacement(
                 rotD3,                // no rotation
                 G4ThreeVector((0.15 * 2.54)*cm, (134.94)*mm, 103.188*mm), // its position
                 absorberLV3,       // its logical volume                         
                 "Abso3",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //
  // For the 2 side walls and 1 Back side wall ( Just for the SDSMT setup not for Underground)
  // gapS3, gapS4 & gapS5 name is used but energy deposition is not tracked)
  //
  // Gap3 --> side wall 1

  auto gapS3 
    = new G4Box("Gap3",             // its name
                 calorSizeX/2., sidePlateY/2., sidePlateZ/2.); // its size
                         
  auto gapLV3
    = new G4LogicalVolume(
                 gapS3,             // its solid
                 gapMaterial,
                 //defaultMaterial,      // its material
                 "Gap3");           // its name
                                   
  fHoldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., calorSizeY/2. + 2.7*cm, -(layerThickness - gapThickness -(3.25*2.54)*cm)), // its position ( - gapThicknessSW)
                 //G4ThreeVector(0., (103.188 + (12.7 + 25.4+4.72))*mm,-(layerThickness - gapThickness -(2.75*2.54)*cm )),
                 gapLV3,            // its logical volume                         
                 "Gap3",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  // Gap4 --> side wall 2
  auto gapS4 
    = new G4Box("Gap4",             // its name
                 calorSizeX/2., sidePlateY/2., sidePlateZ/2.); // its size
                         
  auto gapLV4
    = new G4LogicalVolume(
                 gapS4,             // its solid
                 gapMaterial,
                 //defaultMaterial,      // its material
                 "Gap4");           // its name
                                   
  fHoldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., -(calorSizeY/2. + 2.7*cm), -(layerThickness -gapThickness - (3.25*2.54)*cm)), // its position ( - gapThicknessSW)
                 gapLV4,            // its logical volume                         
                 "Gap4",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  // Gap5 --> side wall 3 (Back side)
  auto gapS5 
    = new G4Box("Gap5",             // its name
                 BsidePlateX/2., BsidePlateY/2., sidePlateZ/2.); // its size
                         
  auto gapLV5
    = new G4LogicalVolume(
                 gapS5,             // its solid
                 gapMaterial,
                 //defaultMaterial,      // its material
                 "Gap5");           // its name
                                   
  fHoldPV
    = new G4PVPlacement(
                 rotD1,                // no rotation
                 G4ThreeVector(calorSizeY/2. + ((0.5+1.25)*2.54)*cm, 0., -(layerThickness - gapThickness - (3.25*2.54)*cm)), // its position ()gapThicknessB
                 gapLV5,            // its logical volume                         
                 "Gap5",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  
  // print parameters
  //
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << nofLayers << " layers of: [ "
    << gapThickness/mm << "mm of " << absorberMaterial->GetName() 
    << " + "
    << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;
  
  //                                        
  // Visualization attributes
  //

  //G4VisAttributes* grey_solid = new G4VisAttributes(G4Colour::Grey());
  //G4VisAttributes* grey_clear = new G4VisAttributes(G4Colour(G4Colour::Grey().GetRed(), G4Colour::Grey().GetGreen(), G4Colour::Grey().GetBlue(), 0.4));
  G4VisAttributes* blue_clear = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.4));
  G4VisAttributes* aqua_clear = new G4VisAttributes(G4Colour(0.0, 0.5, 1.0, 0.4));
  
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
  calorLV->SetVisAttributes (G4VisAttributes::GetInvisible());
  layerLV->SetVisAttributes (G4VisAttributes::GetInvisible());
  //layerLV->SetVisAttributes(G4Colour(0.0,0.0,0.0)); //Black color

  //logicRing->SetVisAttributes(G4VisAttributes(G4Colour::Yellow()));
  SBoxLV->SetVisAttributes(G4VisAttributes(G4Colour::Yellow()));

  absorberLV->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
  absorberLV2->SetVisAttributes(G4VisAttributes(G4Colour::Red()));
  absorberLV3->SetVisAttributes(G4VisAttributes(G4Colour::Red()));

  gapLVT1->SetVisAttributes(blue_clear);
  gapLVT2->SetVisAttributes(blue_clear);
  gapLVT3->SetVisAttributes(blue_clear);
  //gapLVT1->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));
  //gapLVT2->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));
  //gapLVT3->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));

  gapLV->SetVisAttributes(aqua_clear);
  //gapLV->SetVisAttributes(G4VisAttributes(G4Colour::Green()));

  //auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,2.0,3.0));
  //simpleBoxVisAtt->SetVisibility(true);
  //calorLV->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void He3DetectorConstruction::ConstructSDandField()
{ 
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
