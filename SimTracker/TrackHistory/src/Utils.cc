/**\
   The following connection between m_map and CMSProcess is done using the following twiki page:
   https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMCTruth
   
*/

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimTracker/TrackHistory/interface/CMSProcessTypes.h"
#include "SimTracker/TrackHistory/interface/Utils.h"

G4toCMSLegacyProcTypeMap::G4toCMSLegacyProcTypeMap()
{
  // Geant4Process -> CmsProcess
  m_map[0]   = CMS::Primary;              // "Primary"         -> "Primary"
  m_map[91]  = CMS::Unknown;              // "Transportation"  -> "Unknown"
  m_map[92]  = CMS::Unknown;              // "CoupleTrans"     -> "Unknown"
  m_map[1]   = CMS::Unknown;              // "CoulombScat"     -> "Unknown"
  m_map[2]   = CMS::Ionisation;           // "Ionisation"      -> "Ionisation" 
  m_map[3]   = CMS::Bremsstrahlung;       // "Brems"           -> "Bremsstrahlung"
  m_map[4]   = CMS::EPairProd;            // "eePairProd"      -> "EPairProd"
  m_map[5]   = CMS::Annihilation;         // "AnnihIn2Gamma"   -> "Annihilation"
  m_map[6]   = CMS::Annihilation;         // "AnnihToMuMu"     -> "Annihilation"
  m_map[7]   = CMS::Annihilation;         // "AnnihToHad"      -> "Annihilation"
  m_map[10]  = CMS::Unknown;              // "Msc"             -> "Unknown"
  m_map[12]  = CMS::Photon;               // "PhotoElectric"   -> "Photon"
  m_map[13]  = CMS::Compton;              // "Compton"         -> "Compton"
  m_map[14]  = CMS::Conversions;          // "ConvToee"        -> "Conversions"
  m_map[15]  = CMS::Conversions;          // "ConvToMuMu"      -> "Conversions"
  m_map[23]  = CMS::SynchrotronRadiation; // "SynchRad"        -> "SynchrotronRadiation"
  m_map[111] = CMS::Hadronic;             // "HadElastic"      -> "Hadronic"
  m_map[121] = CMS::Hadronic;             // "HadInelastic"    -> "Hadronic"
  m_map[131] = CMS::Neutron;              // "NeutronCapture"  -> "Neutron"
  m_map[141] = CMS::Neutron;              // "NeutronFission"  -> "Neutron"
  m_map[151] = CMS::Hadronic;             // "HadAtRest"       -> "Hadronic"
  m_map[201] = CMS::Decay;                // "Decay"           -> "Decay"
  m_map[202] = CMS::Decay;                // "DecayWSpin"      -> "Decay"
  m_map[203] = CMS::Decay;                // "DecayPiWSpin"    -> "Decay"
  m_map[204] = CMS::Decay;                // "DecayRadio"      -> "Decay"
  m_map[205] = CMS::Decay;                // "DecayUnKnown"    -> "Decay"
  m_map[206] = CMS::Decay;                // "DecayExt"        -> "Decay"
  m_map[301] = CMS::Unknown;              // "GFlash"          -> "Unknown"
  m_map[401] = CMS::Unknown;              // "StepLimiter"     -> "Unknown"
  m_map[403] = CMS::Unknown;              // "NeutronKiller"   -> "Unknown"
}

const unsigned int G4toCMSLegacyProcTypeMap::processId(unsigned int g4ProcessId) const
{
  MapType::const_iterator it = m_map.find(g4ProcessId);

  if( it == m_map.end() )
  {
    edm::LogError("UnknownProcessType") << "Encountered an unknown process type: " << g4ProcessId << ". Mapping it to 'Undefined'.";
    return CMS::Undefined;
  }
  else
    return it->second;
}
