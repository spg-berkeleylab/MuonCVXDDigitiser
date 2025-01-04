#include "MuonCVXDDigitiser.h"

// Standard
#include <iostream>
#include <algorithm>

// edm4hep
#include <edm4hep/MCParticle.h>
#include "BitField64.hxx"

// DD4hep
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DD4hepUnits.h"

// Random
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_math.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h" 
#include "CLHEP/Random/RandFlat.h" 

using CLHEP::RandGauss;
using CLHEP::RandPoisson;
using CLHEP::RandFlat;
using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::rec::ZPlanarData;
using dd4hep::rec::ZDiskPetalsData;
using dd4hep::rec::SurfaceManager;
using dd4hep::rec::SurfaceMap;
using dd4hep::rec::ISurface;
using dd4hep::rec::Vector2D;
using dd4hep::rec::Vector3D;

DECLARE_COMPONENT(MuonCVXDDigitiser)

MuonCVXDDigitiser::MuonCVXDDigitiser(const std::string& name, ISvcLocator* svcLoc) : MultiTransformer(name, svcLoc,
          { KeyValues("CollectionName", {"VertexBarrelCollection"}) },
          { KeyValues("SimHitLocCollectionName", {"VertexBarrel"}),
            KeyValues("OutputCollectionName", {"VTXTrackerHits"}),
            KeyValues("RelationColName", {"VTXTrackerHitRelations"}),
            KeyValues("RawHitsLinkColName", {"VTXRawHitRelations"}) }) {}


StatusCode MuonCVXDDigitiser::initialize() { 
    debug() << "   init called  " << endmsg;
    // Determine if we're handling barrel or endcap geometry
    if (m_subDetName.find("Barrel") != std::string::npos) {
      isBarrel=true;
    } else if (m_subDetName.find("Endcap") != std::string::npos) {
      isBarrel=false;
    } else {
      std::stringstream err;
      err << " Could not determine sub-detector type for: " << m_subDetName;
      throw Exception ( err.str() );
    }

    // Determine if vertex, inner tracker, or outer tracker
    if (m_subDetName.find("Vertex") != std::string::npos) {
      isVertex=true;
    } else if (m_subDetName.find("InnerTracker") != std::string::npos) {
      isInnerTracker=true;
    } else if (m_subDetName.find("OuterTracker") != std::string::npos) {
      isOuterTracker=true;
    } else {
      std::stringstream err;
      err << " Could not determine sub-detector type for: " << m_subDetName;
      throw Exception ( err.str() );
    }

    m_fluctuate = new MyG4UniversalFluctuationForSi();
    
    return StatusCode::SUCCESS;
}

void MuonCVXDDigitiser::LoadGeometry() const{
    Detector& theDetector = Detector::getInstance();
    DetElement subDetector = theDetector.detector(m_subDetName);
    std::vector<ZPlanarData::LayerLayout> barrelLayers;
    std::vector<ZDiskPetalsData::LayerLayout> endcapLayers;
    if (isBarrel) {
      ZPlanarData*  zPlanarData=nullptr;
      // Barrel-like geometry
      zPlanarData = subDetector.extension<ZPlanarData>();
      if (! zPlanarData) {
	    std::stringstream err;
            err << " Could not find surface of type ZPlanarData for subdetector: " << m_subDetName;
	    throw Exception ( err.str() );
      }
      barrelLayers = zPlanarData->layers;
      m_numberOfLayers  = barrelLayers.size();
    } else {
      ZDiskPetalsData *zDiskPetalData = nullptr;
      //Endcap-like geometry
      zDiskPetalData = subDetector.extension<ZDiskPetalsData>();
      if (! zDiskPetalData) {
	std::stringstream err;
        err << " Could not find surface of type ZDiskPetalsData for subdetector: " << m_subDetName;
	throw Exception ( err.str() );	
      }
      endcapLayers = zDiskPetalData->layers;
      m_numberOfLayers = endcapLayers.size();
    } 
    SurfaceManager& surfMan = *theDetector.extension<SurfaceManager>();
    m_map = surfMan.map( subDetector.name() ) ;
    if( ! m_map ) 
    {
      std::stringstream err;
      err << " Could not find surface map for detector: "
          << m_subDetName << " in SurfaceManager " ;
      throw Exception( err.str() );
    }
    m_laddersInLayer.resize(m_numberOfLayers);
#ifdef ZSEGMENTED
    m_sensorsPerLadder.resize(m_numberOfLayers);
#endif
    m_layerHalfPhi.resize(m_numberOfLayers);
    m_layerHalfThickness.resize(m_numberOfLayers);
    m_layerThickness.resize(m_numberOfLayers);
    m_layerRadius.resize(m_numberOfLayers);
    m_layerLadderLength.resize(m_numberOfLayers);
    m_layerLadderWidth.resize(m_numberOfLayers);
    m_layerLadderHalfWidth.resize(m_numberOfLayers);
    m_layerActiveSiOffset.resize(m_numberOfLayers);
    m_layerPhiOffset.resize(m_numberOfLayers);
    m_petalsInLayer.resize(m_numberOfLayers);
    m_layerPetalLength.resize(m_numberOfLayers);
    m_layerPetalInnerWidth.resize(m_numberOfLayers);
    m_layerPetalOuterWidth.resize(m_numberOfLayers);
    int curr_layer = 0;
    if (isBarrel) {
      for(ZPlanarData::LayerLayout z_layout : barrelLayers)
	{
	  // ALE: Geometry is in cm, convert all lenght in mm
	  // ALE: Geometry is in cm, convert all length to mm
	  m_laddersInLayer[curr_layer] = z_layout.ladderNumber;
	  m_layerHalfPhi[curr_layer] = M_PI / ((double)_laddersInLayer[curr_layer]) ;
	  m_layerThickness[curr_layer] = z_layout.thicknessSensitive * dd4hep::cm / dd4hep::mm ;
	  m_layerHalfThickness[curr_layer] = 0.5 * m_layerThickness[curr_layer];
	  m_layerRadius[curr_layer] = z_layout.distanceSensitive * dd4hep::cm / dd4hep::mm  + m_layerHalfThickness[curr_layer];
#ifdef ZSEGMENTED
	  m_sensorsPerLadder[curr_layer] = z_layout.sensorsPerLadder;
	  m_layerLadderLength[curr_layer] = z_layout.lengthSensor * z_layout.sensorsPerLadder * dd4hep::cm / dd4hep::mm ;
#else
	  m_layerLadderLength[curr_layer] = z_layout.lengthSensor * dd4hep::cm / dd4hep::mm ;
#endif
     if (!isVertex) { m_layerLadderLength[curr_layer] = 2 * z_layout.zHalfSensitive;}
	  m_layerLadderWidth[curr_layer] = z_layout.widthSensitive * dd4hep::cm / dd4hep::mm ;	  
	  m_layerLadderHalfWidth[curr_layer] = m_layerLadderWidth[curr_layer] / 2.;
	  m_layerActiveSiOffset[curr_layer] = - z_layout.offsetSensitive * dd4hep::cm / dd4hep::mm ;
	  m_layerPhiOffset[curr_layer] = z_layout.phi0;
	  
	  curr_layer++;
	}
    } else {
      // TODO.. fill info necessary for endcap
      for (ZDiskPetalsData::LayerLayout z_layout : endcapLayers)
	{
	  // see /opt/ilcsoft/muonc/lcgeo/v00-18-01-MC/detector/tracker/VertexEndcap_o1_v06_geo.cpp, L144+
	  // for how the endcap (VXD) geometry is built
	  // Note: petal-like structure but current geometry only defines a single sensitive element for the whole disk afaics.
	  // The structure is defined in /opt/ilcsoft/muonc/DD4hep/v01-25-01/DDRec/include/DDRec/DetectorData.h

      //_petalsInLayer[curr_layer] = z_layout.ladderNumber;
	  //_layerHalfPhi[curr_layer] = M_PI / ((double)_laddersInLayer[curr_layer]) ;
	  m_layerThickness[curr_layer] = z_layout.thicknessSensitive * dd4hep::cm / dd4hep::mm ;
	  m_layerHalfThickness[curr_layer] = 0.5 * m_layerThickness[curr_layer];

	  //_sensorsPerPetal[curr_layer] = z_layout.sensorsPerPetal;
      // CS: does the petal sensitive length include all petal sub-sensors?
	  m_layerPetalLength[curr_layer] = z_layout.lengthSensitive * dd4hep::cm / dd4hep::mm ;

	  m_layerPetalInnerWidth[curr_layer] = z_layout.widthInnerSensitive * dd4hep::cm / dd4hep::mm ;
          m_layerPetalOuterWidth[curr_layer] = z_layout.widthOuterSensitive * dd4hep::cm / dd4hep::mm ;
          m_petalsInLayer[curr_layer] = z_layout.petalNumber;
          if (m_layerPetalOuterWidth[curr_layer] == 0){
            float outerEndcapRadius = 112.0 * dd4hep::cm / dd4hep::mm ;// FIX find source
            m_layerPetalOuterWidth[curr_layer] = 2 * outerEndcapRadius * std::tan(M_PI/m_petalsInLayer[curr_layer]);
          }
	  //_layerActiveSiOffset[curr_layer] = - z_layout.offsetSensitive * dd4hep::cm / dd4hep::mm ;
	  //_layerPhiOffset[curr_layer] = z_layout.phi0;
	
	  curr_layer++;
	}
    }

    // temp fix for issue: z_layout.lengthSensor = 0 for the inner tracker barrel
    // manually hard code values for ladderlength
    /* if (m_layerLadderLength[0] == 0 && isInnerTracker && isBarrel){
        m_layerLadderLength = {963.2,963.2,1384.6};
    }
    if (m_layerLadderLength[0] == 0 && $$ isOuterTracker && isBarrel){
        m_layerLadderLength = {1264.2*2,1264.2*2,1264.2*2};
    } */

    PrintGeometryInfo();

    // Bins for charge discretization
    // FIXME: Will move to assign more dynamically 
    if (m_ChargeDigitizeNumBits == 3) m_DigitizedBins = {500, 786, 1100, 1451, 1854, 2390, 3326, 31973};
    
    //  Here is the updated version of the m_DigitizedBins w/ 4 bits
    if (m_ChargeDigitizeNumBits == 4) m_DigitizedBins = {500, 657, 862, 1132, 1487, 1952, 2563, 3366, 4420, 5804, 7621, 10008, 13142, 17257, 22660, 29756}; //{500, 639, 769, 910, 1057, 1213, 1379, 1559, 1743, 1945, 2193, 2484, 2849, 3427, 4675, 29756};    
    
    if (m_ChargeDigitizeNumBits == 5) m_DigitizedBins = {500, 573, 633, 698, 757, 821, 890, 963, 1032, 1104, 1179, 1260, 1337, 1421, 1505, 1600, 1685, 1777, 1875, 1982, 2097, 2220, 2352, 2511, 2679, 2866, 3107, 3429, 3880, 4618, 6287, 16039};
    if (m_ChargeDigitizeNumBits == 6) m_DigitizedBins = {500, 542, 572, 601, 629, 661, 692, 721, 750, 779, 812, 842, 877, 913, 946, 981, 1016, 1051, 1087, 1121, 1161, 1196, 1237, 1275, 1313, 1350, 1391, 1431, 1468, 1514, 1560, 1606, 1646, 1687, 1733, 1777, 1821, 1872, 1920, 1976, 2036, 2091, 2145, 2213, 2272, 2337, 2411, 2488, 2573, 2651, 2739, 2834, 2938, 3053, 3194, 3356, 3532, 3764, 4034, 4379, 4907, 5698, 6957, 9636};
    if (m_ChargeDigitizeNumBits == 8) m_DigitizedBins = {500, 511, 523, 533, 542, 550, 556, 564, 570, 577, 585, 592, 598, 603, 610, 617, 624, 630, 638, 646, 654, 661, 668, 676, 684, 691, 699, 705, 712, 719, 724, 731, 738, 745, 752, 760, 767, 772, 780, 787, 795, 802, 810, 818, 826, 832, 839, 847, 856, 865, 874, 881, 889, 898, 906, 916, 924, 930, 938, 945, 955, 965, 971, 978, 986, 995, 1004, 1012, 1019, 1027, 1036, 1044, 1053, 1062, 1071, 1079, 1088, 1096, 1104, 1112, 1121, 1131, 1139, 1149, 1158, 1168, 1175, 1184, 1193, 1203, 1211, 1221, 1233, 1241, 1249, 1259, 1268, 1277, 1286, 1294, 1303, 1313, 1321, 1330, 1338, 1348, 1357, 1368, 1378, 1387, 1395, 1406, 1417, 1426, 1434, 1445, 1452, 1460, 1470, 1480, 1492, 1503, 1514, 1525, 1536, 1550, 1560, 1570, 1580, 1592, 1604, 1614, 1623, 1634, 1644, 1653, 1662, 1673, 1684, 1695, 1707, 1717, 1727, 1737, 1747, 1759, 1769, 1780, 1790, 1800, 1812, 1823, 1835, 1846, 1860, 1873, 1885, 1897, 1907, 1918, 1931, 1943, 1958, 1971, 1987, 2000, 2014, 2026, 2041, 2056, 2068, 2080, 2095, 2108, 2119, 2131, 2147, 2162, 2180, 2195, 2213, 2224, 2238, 2256, 2269, 2284, 2300, 2314, 2332, 2351, 2366, 2383, 2401, 2421, 2440, 2458, 2475, 2496, 2519, 2538, 2559, 2581, 2601, 2618, 2636, 2658, 2681, 2703, 2722, 2742, 2767, 2791, 2811, 2836, 2857, 2884, 2913, 2938, 2967, 2995, 3023, 3052, 3086, 3119, 3153, 3188, 3221, 3270, 3304, 3342, 3390, 3428, 3473, 3515, 3556, 3611, 3691, 3742, 3801, 3857, 3928, 3999, 4069, 4141, 4220, 4325, 4417, 4518, 4655, 4789, 4965, 5141, 5359, 5548, 5770, 6017, 6311, 6584, 7024, 7492, 8060, 8740, 9738, 11450, 14878, 23973};

    // shift digitized bins for inner and outer tracker by factor of 2
    // this adjusts for the fact that the resolution is 2x worse for inner and outer tracker
    if (!isVertex) {
        debug() << "Subdetector is: " << m_subDetName << endmsg;
        float shift = 500.; // first bin
        float scalefactor = 2.; 
        for (int i = 0; i < m_DigitizedBins.size(); i++){
            m_DigitizedBins[i] = (m_DigitizedBins[i] - m_DigitizedBins[0]) * scalefactor + shift;
        }
    }
} 

std::tuple<edm4hep::SimTrackerHitCollection,
           edm4hep::TrackerHitPlaneCollection,
           edm4hep::TrackerHitSimTrackerHitCollection,
           edm4hep::TrackerHitSimTrackerHitCollection> MuonCVXDDigitiser::operator(
     const edm4hep::SimTrackerHitCollection& STHcol) const{ 
    //SP. few TODO items:
    // - include noisy pixels (calculate rate from gaussian with unit sigma integral x > m_electronicNoise / m_threshold )
    // - change logic in creating pixels from all SimTrkHits, then cluster them (incl. timing info)
    // - include threshold dispersion effects
    // - add digi parametrization for time measurement
    // - change position determination of cluster to analog cluster (w-avg of corner hits)
    if (!m_map){
        LoadGeometry() ;
    }
    edm4hep::SimTrackerHitCollection               STHLocCol;
    edm4hep::TrackerHitPlaneCollection             THcol;
    edm4hep::TrackerHitSimTrackerHitLinkCollection relCol;
    edm4hep::TrackerHitSimTrackerHitLinkCollection rawHitsCol;

    BitField64 cellID_coder("subdet:5,side:-2,layer:9,module:8,sensor:8");

    int nSimHits = STHcol.size();
    debug() << "Processing collection " << STHcol.getID()  << " with " <<  nSimHits  << " hits ... " << endmsg;
    for (int i=0; i < nSimHits; ++i) {
        edm4hep::SimTrackerHit simTrkHit = STHcol.at(i);
        InternalState intState;
        // use CellID to set layer and ladder numbers
        cellID_coder.setValue(simTrkHit.getCellID());
        intState.currentLayer = cellID_coder["layer"];
        intState.currentLadder = cellID_coder["module"];
        debug() << "Processing simHit #" << i << ", from layer=" << intState.currentLayer << ", module=" << intState.currentLadder << "\n"
                << "- EDep = " << simTrkHit.getEDep() *dd4hep::GeV / dd4hep::keV << " keV, path length = " 
                << simTrkHit.getPathLength() * 1000. << " um" << endmsg;
        float mcp_r = std::sqrt(simTrkHit.getPosition().x *simTrkHit.getPosition().x +simTrkHit.getPosition().y *simTrkHit.getPosition().y);
        float mcp_phi = std::atan(simTrkHit.getPosition().y /simTrkHit.getPosition().x);
        float mcp_theta = simTrkHit.getPosition().z == 0 ? 3.1416/2 : std::atan(mcp_r/simTrkHit.getPosition().z);
        debug() << "- Position (mm) x,y,z,t = " << simTrkHit.getPosition().x << ", " 
                                                << simTrkHit.getPosition().y << ", " 
                                                << simTrkHit.getPosition().z << ", " 
                                                << simTrkHit.getTime() << "\n" 
                << "- Position r(mm),phi,theta = " << mcp_r << ", " << mcp_phi << ", " << mcp_theta 
                << "\n- MC particle pdg = ";
        edm4hep::MCParticle mcp = simTrkHit.getParticle();
        if (&mcp) {
	    debug() << mcp.getPDG();
        } else {
	    debug()<< " N.A.";
        }
        debug() <<  "\n- MC particle p (GeV) = " << std::sqrt(simTrkHit.getMomentum().x*simTrkHit.getMomentum().x+simTrkHit.getMomentum().y*simTrkHit.getMomentum().y+simTrkHit.getMomentum().z*simTrkHit.getMomentum().z)
                << "\n- isSecondary = " << simTrkHit.isProducedBySecondary() << ", isOverlay = " << simTrkHit.isOverlay()
                << "\n- Quality = " << simTrkHit.getQuality() << endmsg;
        ProduceIonisationPoints( simTrkHit, &intState );
        if (intState.currentLayer == -1)
          continue;
        ProduceSignalPoints(&intState);
        SimTrackerHitImplVec simTrkHitVec;
        ProduceHits(simTrkHitVec, simTrkHit, &intState);
        if (m_PoissonSmearing) PoissonSmearer(simTrkHitVec);
        if (m_electronicEffects) GainSmearer(simTrkHitVec);
        ApplyThreshold(simTrkHitVec);
        if (m_DigitizeCharge) ChargeDigitizer(simTrkHitVec);
        if (m_timeSmearingSigma > 0) TimeSmearer(simTrkHitVec);
        if (m_DigitizeTime) TimeDigitizer(simTrkHitVec);
	    
        //**************************************************************************
        // Create reconstructed cluster object (TrackerHitImpl)
        //**************************************************************************
        edm4hep::MutableTrackerHitPlane *recoHit = ReconstructTrackerHit(simTrkHitVec, &THcol, &intState) {
        if (recoHit == nullptr)
          debug() << "Skip hit" << endmsg;
          continue;
        }       
        // hit's layer/ladder/petal position does not change
        const int cellid = simTrkHit.getCellID();
        recoHit->setCellID( cellid );
            
        edm4hep::Vector3d localPos;
        edm4hep::Vector3d localIdx;
        edm4hep::Vector3d localDir;
        FindLocalPosition(simTrkHit, localPos, localDir, &intState);
        localIdx.x = localPos.x / m_pixelSizeX;
        localIdx.y = localPos.y / m_pixelSizeY;
        float incidentPhi = std::atan(localDir.x / localDir.z);
        float incidentTheta = std::atan(localDir.y / localDir.z);

        // Debug messages to check if reconstruction went correctly
        // true global
        debug() << "- TRUE GLOBAL position (mm) x,y,z,t = " << simTrkHit.getPosition().x << ", " 
                                                            << simTrkHit.getPosition().y << ", " 
                                                            << simTrkHit.getPosition().z << ", " 
                                                            << simTrkHit.getTime() << "\n"
        // true local (compare two verions)
                << "- TRUE LOCAL position (localPos) (mm) x,y,z,t = " << localPos.x << ", " 
                                                                      << localPos.y << ", " 
                                                                      << localPos.z << "\n"
        // reco local 
                << "- RECO LOCAL position (mm) x,y,z,t = " << recoHit->getPosition().x << ", "
                                                           << recoHit->getPosition().y << ", "
                                                           << recoHit->getPosition().z << "\n" << endmsg;
            
        edm4hep::Vector3d xLab;
        TransformToLab( cellid, recoHit->getPosition(), xLab);
        recoHit->setPosition( xLab );

        // reco global
        debug() << "- RECO GLOBAL position (mm) x,y,z,t = " << recoHit->getPosition().x << ", " 
                                                            << recoHit->getPosition().y << ", " 
                                                            << recoHit->getPosition().z << endmsg;
            
        SurfaceMap::const_iterator sI = m_map->find( cellid );
        const dd4hep::rec::ISurface* surf = sI->second;
        dd4hep::rec::Vector3D u = surf->u() ;
        dd4hep::rec::Vector3D v = surf->v() ;
            
            
        edm4hep::Vector2f u_direction;
        //TODO HACK: Store incidence angle of particle instead!
        u_direction.a = u.theta();
        u_direction.b = u.phi();
        edm4hep::Vector2f v_direction;
        v_direction.a = v.theta();
        v_direction.b = v.phi();
        recoHit->setU( u_direction );
        recoHit->setV( v_direction );
            
        //**************************************************************************
        // Set Relation to SimTrackerHit
        //**************************************************************************    
        edm4hep::MutableTrackerHitSimTrackerHitLink rel = relCol.create();
        rel.setFrom(*recoHit);
        rel.setTo(simTrkHit);
        rel.setWeight( 1.0 );
        debug() << "Reconstructed pixel cluster:\n"
                << "- local position (x,y) = " << localPos.x << "(Idx: " << localIdx.x << "), " 
                                               << localPos.y << "(Idy: " << localIdx.y << ")\n"
                << "(reco local) - (true local) (x,y,z): " << localPos.x - intState.currentLocalPosition.x << ", " 
                                                           << localPos.t - intState.currentLocalPosition.y << ", " 
                                                           << localPos.z - intState.currentLocalPosition.z << "\n"
                << "- global position (x,y,z, t) = " << recoHit->getPosition().x << ", "
                                                     << recoHit->getPosition().y << ", "
                                                     << recoHit->getPosition().z << ", " 
                                                     << recoHit->getTime() << "\n"
                << "- (reco global (x,y,z,t)) - (true global) = " << recoHit->getPosition().x - simTrkHit.getPosition().x << ", "
                                                                  << recoHit->getPosition().y - simTrkHit.getPosition().y << ", "
                                                                  << recoHit->getPosition().z - simTrkHit.getPosition().z << ", "
                                                                  << recoHit->getTime() - simTrkHit.getTime() << "\n"
                << "- charge = " << recoHit->getEDep() << "(True: " << simTrkHit.getEDep() << ")\n"
                << "- incidence angles: theta = " << incidentTheta << ", phi = " << incidentPhi << endmsg;
        
        std::vector<edm4hep::SimTrackerHit*> rawHits;
        if (m_produceFullPattern != 0) {
          // Store all the fired points
          for (int iS = 0; iS < (int)simTrkHitVec.size(); ++iS) {
            edm4hep::MutableSimTrackerHit *sth = simTrkHitVec[iS];
            float charge = sth->getEDep();
            //store hits that are above threshold. In case of m_ChargeDiscretization, just check for a small non-zero value                
            if ( (m_DigitizeCharge and (charge >1.0)) or (charge > m_threshold) ) {
              edm4hep::MutableSimTrackerHit newsth = STHLocCol.create();
              // hit's layer/ladder position is the same for all fired points 
              newsth.setCellID( cellid );
              //Store local position in units of pixels instead
              const double *sLab;
              //TransformToLab(cellid0, sth->getPosition(), sLab);
              sLab = sth->getPosition();
              edm4hep::Vector3d pixelPos;
              pixelPos.x = sLab.x / m_pixelSizeX;
              pixelPos.y = sLab.y / m_pixelSizeY;
              newsth.setPosition(pixelPos);
              newsth.setEDep(charge); // in unit of electrons
              newsth.setTime(sth->getTime());
              newsth.setPathLength(simTrkHit.getPathLength());
              newsth.setMCParticle(simTrkHit.getParticle());
              newsth.setMomentum(simTrkHit.getMomentum());
              newsth.setProducedBySecondary(simTrkHit.isProducedBySecondary());
              newsth.setOverlay(simTrkHit.isOverlay());
              
              numOfPixels++;
              rawHits.push_back(&newsth);
            }
          }
        }
        debug() << "\n- number of pixels: " << rawHits.size()
                << "\n- MC particle p=" << std::sqrt(simTrkHit.getMomentum().x*simTrkHit.getMomentum().x+simTrkHit.getMomentum().y*simTrkHit.getMomentum().y+simTrkHit.getMomentum().z*simTrkHit.getMomentum().z)
                << "\n- isSecondary = " << simTrkHit.isProducedBySecondary() << ", isOverlay = " << simTrkHit.isOverlay()
                << "\n- List of constituents (pixels/strips):";
        for (size_t iH = 0; iH < rawHits.size(); ++iH) {
          edm4hep::MutableTrackerHitSimTrackerHitLink rawLink = rawHitsCol.create();
          rawLink.setFrom(*recoHit);
          rawLink.setTo(rawHits.at(iH));
          rawLink.setWeight(1. / rawHits.size());
          debug() << "  - " << iH << ": Edep (e-) = " << rawHits.at(iH)->getEDep() << ", t (ns) =" << rawHits.at(iH)->getTime();
        }
        debug() << "--------------------------------\n";
        for (int k=0; k < int(simTrkHitVec.size()); ++k) {
          SimTrackerHit *hit = simTrkHitVec[k];
          delete hit;
        }
        delete intState
    }
    debug() << "Number of produced hits: " << THcol.size()  << endmsg;
    
    return std::make_tuple(std::move(STHLocCol),
                           std::move(THcol),
                           std::move(relCol),
                           std::move(rawHitsCol) );
}

Status Code MuonCVXDDigitiser::finalize() {
    debug() << "   end called  " << endmsg;
    delete m_fluctuate;
}

/** Function calculates local coordinates of the sim hit 
 * in the given ladder and local momentum of particle. 
 * Also returns module number and ladder number.
 * Local coordinate system within the ladder 
 * is defined as following :  <br> 
 *    - x axis lies in the ladder plane and orthogonal to the beam axis <br>
 *    - y axis lies in the ladder plane and parallel to the beam axis <br>
 *    - z axis is perpendicular to the ladder plane <br>
 * 
 */
void MuonCVXDDigitiser::FindLocalPosition(edm4hep::SimTrackerHit &hit, 
                                          edm4hep::Vector3d &localPosition,
                                          edm4hep::Vector3d &localDirection,
                                          InternalState *intState) const{
    // Use SurfaceManager to calculate local coordinates
    const int cellID = hit.getCellID() ;
    debug() << "Cell ID of Sim Hit: " << cellID << endmsg;
    SurfaceMap::const_iterator sI = m_map->find( cellID ) ;
    const dd4hep::rec::ISurface* surf = sI->second ;
    Vector3D oldPos( hit.getPosition().x, hit.getPosition().y, hit.getPosition().z );
    // We need it?
    if ( ! surf->insideBounds( dd4hep::mm * oldPos ) ) {
        debug() << "  hit at " << oldPos
                << " is not on surface "
                << *surf
                << " distance: " << surf->distance(  dd4hep::mm * oldPos )
                << endmsg;
      intState->currentLayer = -1;
      return;
    }    
    
    
    Vector2D lv = surf->globalToLocal( dd4hep::mm * oldPos  ) ;
    // Store local position in mm
    localPosition.x = lv[0] / dd4hep::mm ;
    localPosition.y = lv[1] / dd4hep::mm ;
    // Add also z ccordinate
    Vector3D origin( surf->origin()[0], surf->origin()[1], surf->origin()[2]);
    localPosition.z = ( dd4hep::mm * oldPos - dd4hep::cm * origin ).dot( surf->normal() ) / dd4hep::mm;
    edm4hep::Vector3d Momentum;
    edm4hep::MCParticle mcp = hit.getParticle();
    for (int j = 0; j < 3; ++j) {
      if (&mcp) {
        Momentum[j] = mcp.getMomentum()[j] * dd4hep::GeV;
      } else {
        Momentum[j] = hit.getMomentum()[j];
      }
    }
    // as default put electron's mass
    intState->currentParticleMass = 0.510e-3 * dd4hep::GeV;
    if (&(hit.getParticle()))
        intState->currentParticleMass = std::max(hit.getParticle().getMass() * dd4hep::GeV, intState->currentParticleMass);
    intState->currentParticleMomentum = sqrt(pow(Momentum.x, 2)
                                   + pow(Momentum.y, 2)
                                   + pow(Momentum.z, 2));                   
                         
    localDirection.x = Momentum * surf->u();
    localDirection.y = Momentum * surf->v();
    localDirection.z = Momentum * surf->normal();
    if (isBarrel){
      intState->currentPhi = intState->currentLadder * 2.0 * m_layerHalfPhi[intState->currentLayer] + m_layerPhiOffset[intState->currentLayer];
    }
}

void MuonCVXDDigitiser::ProduceIonisationPoints(edm4hep::SimTrackerHit &hit, InternalState *intState) const{
    debug() << "Creating Ionization Points" << endmsg;
    edm4hep::Vector3d pos(0,0,0);
    edm4hep::Vector3d dir(0,0,0);
    edm4hep::Vector3d entry;
    edm4hep::Vector3d exit;
    // hit and pos are in mm
    FindLocalPosition(hit, pos, dir, &intState);
    if ( intState->currentLayer == -1)
      return;
 
    entry.z = -m_layerHalfThickness[intState->currentLayer]; 
    exit.z = m_layerHalfThickness[intState->currentLayer];
    // entry points: hit position is in middle of layer. ex: entry_x = x - (z distance to bottom of layer) * px/pz
    for (int i = 0; i < 2; ++i) {
        entry[i] = pos[i] + dir[i] * (entry.z - pos.z) / dir.z;
        exit[i]= pos[i] + dir[i] * (exit.z - pos.z) / dir.z;
    }
    for (int i = 0; i < 3; ++i) {
        intState->currentLocalPosition[i] = pos[i];
        intState->currentEntryPoint[i] = entry[i];
        intState->currentExitPoint[i] = exit[i];
    }
    debug() << "local position: " << intState->currentLocalPosition.x << ", "
                                  << intState->currentLocalPosition.y << ", "
                                  << intState->currentLocalPosition.z << endmsg;
    double tanx = dir.x / dir.z;
    double tany = dir.y / dir.z;  
    
    // trackLength is in mm -> limit length at 1cm
    double trackLength = std::min(m_maxTrkLen,
         m_layerThickness[intState->currentLayer] * sqrt(1.0 + pow(tanx, 2) + pow(tany, 2)));
  
    intState->numberOfSegments = ceil(trackLength / m_segmentLength );
    double dEmean = (dd4hep::keV * m_energyLoss * trackLength) / ((double)(intState->numberOfSegments));
    intState->ionisationPoints.resize(intState->numberOfSegments);
    debug() <<  "Track path length: " << trackLength << ", calculated dEmean * N_segment = " << dEmean << " * " << intState->numberOfSegments << " = " << dEmean*intState->numberOfSegments << endmsg;
    intState->eSum = 0.0;
    // TODO m_segmentLength may be different from segmentLength, is it ok?
    double segmentLength = trackLength / ((double)(intState->numberOfSegments));
    intState->segmentDepth = m_layerThickness[intState->currentLayer] / ((double)(intState->numberOfSegments));
    double z = -m_layerHalfThickness[intState->currentLayer] - 0.5 * intState->segmentDepth;
    
    double hcharge = ( hit.getEDep() / dd4hep::GeV ); 	
    debug() << "Number of ionization points: " << intState->numberOfSegments << ", G4 EDep = "  << hcharge << endmsg;
    for (int i = 0; i < intState->numberOfSegments; ++i) {
        z += intState->segmentDepth;
        double x = pos.x + tanx * (z - pos.z);
        double y = pos.y + tany * (z - pos.z);
        // momentum in MeV/c, mass in MeV, tmax (delta cut) in MeV, 
        // length in mm, meanLoss eloss in MeV.
        double de = m_fluctuate->SampleFluctuations(double(intState->currentParticleMomentum * dd4hep::keV / dd4hep::MeV),
                                                    double(intState->currentParticleMass * dd4hep::keV / dd4hep::MeV),
                                                    m_cutOnDeltaRays,
                                                    segmentLength,
                                                    double(dEmean / dd4hep::MeV)) * dd4hep::MeV;
        intState->eSum += de;
        IonisationPoint ipoint;
        ipoint.eloss = de;
        ipoint.x = x;
        ipoint.y = y;
        ipoint.z = z;
        intState->ionisationPoints[i] = ipoint;
        debug() << " " << i << ": z=" << z << ", eloss = " << de << "(total so far: " << intState->eSum << "), x=" << x << ", y=" << y << endmsg;
    }
   
    const double thr = m_deltaEne/m_electronsPerKeV * dd4hep::keV;
    while ( hcharge > intState->eSum + thr ) {
      // Add additional charge sampled from an 1 / n^2 distribution.
      // Adjust charge to match expectations
      const double       q = randomTail( thr, hcharge - intState->eSum );
      const unsigned int h = floor(RandFlat::shoot(0.0, (double)(intState->numberOfSegments) ));
      intState->ionisationPoints[h].eloss += q;
      intState->eSum += q;
    }
    debug() << "Padding each segment charge (1/n^2 pdf) until total below " << m_deltaEne << "e- threshold. New total energy: " << intState->eSum
            << "\nList of ionization points:";
    for (int i =0; i < intState->numberOfSegments; ++i) {
        debug() << "\n- " << i << ": E=" << intState->ionisationPoints[i].eloss 
                               << ", x=" << intState->ionisationPoints[i].x 
                               << ", y=" << intState->ionisationPoints[i].y 
                               << ", z=" << intState->ionisationPoints[i].z;
    }
    debug() << endmsg;
}

void MuonCVXDDigitiser::ProduceSignalPoints(InternalState *intState) const{
    intState->signalPoints.resize(intState->numberOfSegments);
    // run over ionisation points
    debug() << "Creating signal points" << endmsg;
    for (int i = 0; i < intState->numberOfSegments; ++i) {
        IonisationPoint ipoint = intState->ionisationPoints[i]; // still local coords
        double z = ipoint.z;
        double x = ipoint.x;
        double y = ipoint.y;
        double DistanceToPlane = m_layerHalfThickness[intState->currentLayer] - z;
        double xOnPlane = x + m_tanLorentzAngleX * DistanceToPlane;
        double yOnPlane = y + m_tanLorentzAngleY * DistanceToPlane;
        // For diffusion-coeffieint calculation, see e.g. https://www.slac.stanford.edu/econf/C060717/papers/L008.PDF
        // or directly Eq. 13 of https://cds.cern.ch/record/2161627/files/ieee-tns-07272141.pdf
        // diffusionCoefficient = sqrt(2*D / mu / V), where
        //  - D = 12 cm^2/s // diffusion constant
        //  - mu = 450 cm^2/s/V // mobility
        //  - V = 10-30 V // expected depletion voltage
        //  => m_diffusionCoefficient = 0.04-0.07
        // and diffusion sigma = m_diffusionCoefficient * DistanceToPlane
        // e.g. fot 50um diffusion sigma = 2.1 - 3.7 um
        //double DriftLength = DistanceToPlane * sqrt(1.0 + pow(m_tanLorentzAngleX, 2) 
        //                                                + pow(m_tanLorentzAngleY, 2));
        double SigmaDiff = DistanceToPlane * m_diffusionCoefficient;
        double SigmaX = SigmaDiff * sqrt(1.0 + pow(m_tanLorentzAngleX, 2));
        double SigmaY = SigmaDiff * sqrt(1.0 + pow(m_tanLorentzAngleY, 2));
        // energy is in keV
        double charge = (ipoint.eloss / dd4hep::keV) * m_electronsPerKeV;
        SignalPoint  spoint;
        spoint.x = xOnPlane;
        spoint.y = yOnPlane;
        spoint.sigmaX = SigmaX;
        spoint.sigmaY = SigmaY;
        spoint.charge = charge; // electrons x keV
        intState->signalPoints[i] = spoint;
        debug() << "- " << i << ": charge=" << charge 
                << ", x="<< xOnPlane << "(delta=" << xOnPlane - x << ")"
                << ", y="<< yOnPlane << "(delta=" << yOnPlane - y << ")"
                << ", sigmaDiff=" << SigmaDiff
                << ", sigmaX="<< SigmaX << ", sigmay=" << SigmaY << endmsg;
    }
}

void MuonCVXDDigitiser::ProduceHits(MutableSimTrackerHitVec &simTrkVec, edm4hep::SimTrackerHit &simHit, InternalState *intState) const{  
    simTrkVec.clear();
    std::map<int, edm4hep::SimTrackerHit*> hit_Dict;
    debug() << "Creating hits" << endmsg;
    for (int i=0; i < intState->numberOfSegments; ++i) {
        SignalPoint spoint = intState->signalPoints[i];
        double xCentre = spoint.x;
        double yCentre = spoint.y;
        double sigmaX = spoint.sigmaX;
        double sigmaY = spoint.sigmaY;
        double xLo = spoint.x - 3 * spoint.sigmaX;
        double xUp = spoint.x + 3 * spoint.sigmaX;
        double yLo = spoint.y - 3 * spoint.sigmaY;
        double yUp = spoint.y + 3 * spoint.sigmaY;
        
        int ixLo, ixUp, iyLo, iyUp;
        TransformXYToCellID(xLo, yLo, ixLo, iyLo, intState);
        TransformXYToCellID(xUp, yUp, ixUp, iyUp, intState);
        debug() << i << ": Pixel idx boundaries: ixLo=" << ixLo << ", iyLo=" << iyLo  
                                          <<  ", ixUp=" << ixUp << ", iyUp=" << iyUp << endmsg;
        for (int ix = ixLo; ix< ixUp + 1; ++ix) {   
            if ( (ix < 0) or (ix >= GetPixelsInaColumn(intState)) ) {
                debug() << "Pixels in a column: " << GetPixelsInaColumn(intState)
                        << "\nSkipping pixels with ix =" << ix << endmsg;
                continue;
            }
            for (int iy = iyLo; iy < iyUp + 1; ++iy) {
                if ( (iy < 0) or (iy >= GetPixelsInaRow(intState)) ) {
                    debug() << "Pixels in a row: " << GetPixelsInaRow(intState)
                            << "\nSkipping pixels with iy =" << iy << endmsg;
                    continue;
                }
                double xCurrent, yCurrent;
                TransformCellIDToXY(ix, iy, xCurrent, yCurrent, intState);
                
                gsl_sf_result result;
                /*int status = */gsl_sf_erf_Q_e((xCurrent - 0.5 * m_pixelSizeX - xCentre)/sigmaX, &result);
                double LowerBound = 1 - result.val;
                /*status = */gsl_sf_erf_Q_e((xCurrent + 0.5 * m_pixelSizeX - xCentre)/sigmaX, &result);
                double UpperBound = 1 - result.val;
                double integralX = UpperBound - LowerBound;
                /*status = */gsl_sf_erf_Q_e((yCurrent - 0.5 * m_pixelSizeY - yCentre)/sigmaY, &result);
                LowerBound = 1 - result.val;
                /*status = */gsl_sf_erf_Q_e((yCurrent + 0.5 * m_pixelSizeY - yCentre)/sigmaY, &result);
                UpperBound = 1 - result.val;
                double integralY = UpperBound - LowerBound;
                debug() << "Integral x=" << integralX << ", Integral y=" << integralY << ", signal pt charge=" << spoint.charge << endmsg;
                float totCharge = float(spoint.charge * integralX * integralY);
                int pixelID = GetPixelsInaRow(intState) * ix + iy;
              
                auto item = hit_Dict.find(pixelID);
                if (item == hit_Dict.end()) {
                    edm4hep::MutableSimTrackerHit tmp_hit = new edm4hep::MutableSimTrackerHit();
                    edm4hep::Vector3d pos(
                        xCurrent,
                        yCurrent,
                        m_layerHalfThickness[intState->currentLayer]
                    );
                    tmp_hit.setPosition(pos); // still in local coordinates
                    tmp_hit.setCellID(pixelID);                   // workaround: cellID used for pixel index
                    tmp_hit.setEDep(totCharge);
                    tmp_hit.setTime(simHit.getTime()); //usual true timing as starting point
                  
                    hit_Dict.emplace(pixelID, tmp_hit);
                    debug() << "Created new pixel hit at idx=" << ix << ", idy=" << iy << ", charge=" << totCharge << endsmg;
                } else {
                    float edep = item->second->getEDep();
                    edep += totCharge;
                    item->second->setEDep(edep);
                    //TODO: handle multiple times. For now not needed since all deposits arrive at the same true time.
                    debug() << "Updating pixel hit at idx=" << ix << ", idy=" << iy << ", total charge=" << edep << "(delta = " << totCharge << ")" << endmsg;
                }
            }
        }
    }
    debug() << "List of pixel hits created:" << endmsg; // still in local coords
    int idx = 0;
    for(auto item : hit_Dict) {
        simTrkVec.push_back(item.second);       
        debug() << idx++ << ": x=" << item.second->getPosition().x 
                         << ", y=" << item.second->getPosition().y
                         << ", z=" << item.second->getPosition().z
                         << ", EDep = " << item.second->getEDep() << endmsg;
    }
}

/**
 * Function that fluctuates charge (in units of electrons)
 * deposited on the fired pixels according to the Poisson
 * distribution...
 */
void MuonCVXDDigitiser::PoissonSmearer(MutableSimTrackerHitVec &simTrkVec) const{
    debug() << "Adding Poisson smear to charge" << endmsg;
    for (int ihit = 0; ihit < int(simTrkVec.size()); ++ihit) {
        edm4hep::MutableSimTrackerHit *hit = simTrkVec[ihit];
        float charge = hit->getEDep();
        float rng;
        if (charge > 1e+03) // assume Gaussian
        {
            rng = float(RandGauss::shoot(charge, sqrt(charge)));
        }
        else // assume Poisson
        {
            rng = float(RandPoisson::shoot(charge));
        }
        hit->setEDep(rng);
        debug() << ihit << ": x=" << hit->getPosition().x
                        << ", y=" << hit->getPosition().y
                        << ", z=" << hit->getPosition().z
                << ", charge = " << rng << "(delta = " << charge-rng << ")" << endmsg;
    }  
}

/**
 * Simulation of electronic noise.
 */
void MuonCVXDDigitiser::GainSmearer(MutableSimTrackerHitVec &simTrkVec) const{
    debug() << "Adding FE noise smear to charge" << endmsg;
    for (int i = 0; i < (int)simTrkVec.size(); ++i) {
        double Noise = RandGauss::shoot(0., m_electronicNoise);
        edm4hep::MutableSimTrackerHit *hit = simTrkVec[i];
        hit->setEDep(hit->getEDep() + float(Noise));
        debug() << i << ": x=" << hit->getPosition().x 
                     << ", y=" << hit->getPosition().y
                     << ", z=" << hit->getPosition().z 
                << ", charge = " << hit->getEDep() << "(delta = " << Noise << ")" << endmsg;
    }
}

/**
 * Apply threshold.  
 * Sets the charge to 0 if less than the threshold
 * Smears the threshold by a Gaussian if sigma > 0
 */
void MuonCVXDDigitiser::ApplyThreshold(MutableSimTrackerHitVec &simTrkVec) {
   debug() << "Applying threshold" << endmsg;
   float actualThreshold = m_threshold;
   
   for (int i = 0; i < (int)simTrkVec.size(); ++i) {
     edm4hep::MutableSimTrackerHit *hit = simTrkVec[i];
     
     double smear = 0;
     float origCharge = hit->getEDep();
     if (m_thresholdSmearSigma > 0) smear = RandGauss::shoot(0., m_thresholdSmearSigma);
     actualThreshold = actualThreshold + smear;
     if (hit->getEDep() <= actualThreshold) hit->setEDep(0.0);
     
     debug() << i << ": x=" << hit->getPosition().x
                  << ", y=" << hit->getPosition().y
                  << ", z=" << hit->getPosition().z
             << ", new charge = " << hit->getEDep() << ", previous charge = " << origCharge
             << " smeared threshold = " << actualThreshold << "(delta = " << smear << ")" << endmsg;
   }
}

/**
 * Digitizes the charge.
 * Discretization based on number of bits and bin width scheme.
 */
void MuonCVXDDigitiser::ChargeDigitizer(MutableSimTrackerHitVec &simTrkVec) const{
  debug() << "Charge discretization" << endmsg;
  
  float minThreshold = m_threshold;
  float maxThreshold = m_chargeMax;
  //int split = 0.3; -- future use
  int numBins = pow(2, m_ChargeDigitizeNumBits)-1;
  double discCharge = -999; 
  for (int i = 0; i < (int)simTrkVec.size(); ++i) {
    edm4hep::MutableSimTrackerHit *hit = simTrkVec[i];
    float origCharge =	hit->getEDep();
    discCharge = origCharge;
     
    switch(m_ChargeDigitizeBinning) {
        case 0: { // uniform binning
            if (origCharge < 1.0) break;
            float binWidth = (maxThreshold-minThreshold)/(numBins);
            if (origCharge < binWidth) discCharge = (minThreshold+binWidth)/2;
            else if (origCharge > maxThreshold) discCharge = (maxThreshold-binWidth/2);
            else discCharge = ((ceil((origCharge-binWidth)/binWidth)*binWidth)*2+binWidth)/2;
            break;
        }
        case 1: { // variable binning
            if (origCharge < 1.0) break;
            int binVal=-1;
            for(unsigned int idx = 0; idx < m_DigitizedBins.size()-1; idx++) {
                if (m_DigitizedBins[idx+1] > origCharge) {
                    binVal = idx;
                    break;
                }
            }
            if (binVal < 0) discCharge = (m_DigitizedBins[m_DigitizedBins.size()-2] + m_DigitizedBins[m_DigitizedBins.size()-1]) / 2;
            else discCharge = (m_DigitizedBins[binVal] + m_DigitizedBins[binVal+1]) / 2;
            break;
        }
    }
    hit->setEDep(discCharge);
    debug() << i << ": x=" << hit->getPosition().x
                 << ", y=" << hit->getPosition().y
                 << ", z=" << hit->getPosition().z
                 << ", new charge = " << hit->getEDep() << ", previous charge = " << origCharge
                 << ", number of bits = " << m_ChargeDigitizeNumBits
                 << ", binning scheme = " << m_ChargeDigitizeBinning << endmsg;
    }
}

/**
 * Apply effective measurement resolution.
 * TODO: Right now assuming completely uncorrelated resolution across pixels, will need to divide into:
 * - correlated across pixels, uncorrelated across clusters
 * - correlated within the event, un-correlate 
*/
void MuonCVXDDigitiser::TimeSmearer(MutableSimTrackerHitVec &simTrkVec) const{
    debug() << "Adding resolution effect to timing measurements" << endmsg;
    for (int i = 0; i < (int)simTrkVec.size(); ++i) {
        float delta = RandGauss::shoot(0., m_timeSmearingSigma);
        edm4hep::MutableSimTrackerHit *hit = simTrkVec[i];
        hit->setTime(hit->getTime() + delta);
        debug() << i << ": x=" << hit->getPosition().x
                     << ", y=" << hit->getPosition().y
                     << ", z=" << hit->getPosition().z
                     << ", time = " << hit->getTime() << "(delta = " << delta << ")" << endmsg;
    }
}

/**
 * Digitizes the time information.
 * Discretization based on number of bits and bin width scheme.
 */
void MuonCVXDDigitiser::TimeDigitizer(MutableSimTrackerHitVec &simTrkVec) const{
    debug() << "Time discretization" << endmsg;
  
    static const int numBins = pow(2, m_TimeDigitizeNumBits)-1;
    double discTime;
    for (int i = 0; i < (int)simTrkVec.size(); ++i) {
        edm4hep::MutableSimTrackerHit *hit = simTrkVec[i];
        float origTime = hit->getTime();
        discTime = origTime;
     
        switch(m_TimeDigitizeBinning) {
            case 0: // uniform binning
                static const float binWidth = m_timeMax/numBins;
                if (origTime < binWidth) discTime = binWidth/2;
                else if (origTime > m_timeMax) discTime = m_timeMax - binWidth/2;
                else discTime = ((ceil((origTime-binWidth)/binWidth)*binWidth)*2+binWidth)/2;
                break;
            default:
                error() << "Invalid setting for pixel time digitization binning. Retaining original time." << endmsg;
        }
        hit->setTime(discTime);
        debug() << i << ": x=" << hit->getPosition().x
                     << ", y=" << hit->getPosition().y
                     << ", z=" << hit->getPosition().z
                << ", new time = " << hit->getTime() << ", previous time = " << origTime << endmsg;
    } //end loop over pixel cells
}

/**
 * Emulates reconstruction of Tracker Hit 
 * Tracker hit position is reconstructed as weighted average of edge pixels.
 * The position is corrected for Lorentz shift.
 * Time is the arithmetic average of constituents.
 */
edm4hep::MutableTrackerHitPlane *MuonCVXDDigitiser::ReconstructTrackerHit(MutableSimTrackerHitVec &simTrkVec, 
                                                                          edm4hep::TrackerHitPlaneCollection *THcol, 
                                                                          InternalState *intState) const{
    edm4hep::Vector3d pos(0, 0, 0);

    double minX = 99999999;
    double minY = 99999999;
    double maxX = -99999999;
    double maxY = -99999999;
      
    double charge = 0; // total cluster charge
    unsigned int size = 0; // number of pixels in the cluster
    unsigned int edge_size_minx = 0; //number of pixels at the lower edge of cluster in x direction
    unsigned int edge_size_miny = 0; //number of pixels at the lower edge of cluster in y direction
    unsigned int edge_size_maxx = 0; //number of pixels at the upper edge of cluster in x direction
    unsigned int edge_size_maxy = 0; //number of pixels at the upper edge of cluster in y direction

    debug() << "Creating reconstructed cluster" << endmsg;
    double time = 0; //average time

    /* Get extreme positions, currently only implemented for barrel */
    /* Calculate the mean */
    for (size_t iHit=0; iHit < simTrkVec.size(); ++iHit) {
        edm4hep::SimTrackerHit *hit = simTrkVec[iHit];
        //check for non-zero value (pixels below threshold have already been set to zero)
        if (hit->getEDep() < 1.0) continue;
	
        size += 1;
        time += hit->getTime();
        charge += hit->getEDep();
        debug() << iHit << ": Averaging position, x=" << hit->getPosition().x
                                            << ", y=" << hit->getPosition().y
                        << ", weight(EDep)=" << hit->getEDep() << endmsg;

        // calculate min x, min y, max x, max y
        if (hit->getPosition().x < minX) {
            minX = hit->getPosition().x;
        }
        if (hit->getPosition().y < minY) {
            minY = hit->getPosition().y;
        }
        if (hit->getPosition().x > maxX) {
            maxX = hit->getPosition().x;
        }
        if (hit->getPosition().y > maxY) {
            maxY = hit->getPosition().y;
        }
    }

    // Loop over all pixel hits again, find the pixels on the 4 extreme edges
    for (size_t iHit=0; iHit < simTrkVec.size(); ++iHit) {
        edm4hep::SimTrackerHit *hit = simTrkVec[iHit];
        if (hit->getEDep() < 1.0) continue; // ignore pixels below threshold

        if (hit->getPosition()[0] == minX) edge_size_minx += 1;
        if (hit->getPosition()[1] == minY) edge_size_miny += 1;
        if (hit->getPosition()[0] == maxX) edge_size_maxx += 1;
        if (hit->getPosition()[1] == maxY) edge_size_maxy += 1;
    }

    /* Calculate mean x and y by weighted ave: 
    x_reco = ( (x_max * max edge size) + (x_min * min edge size) ) / (min edge size + max edge size) */
    pos.x = ((minX * edge_size_minx) + (maxX * edge_size_maxx))/(edge_size_minx + edge_size_maxx);
    pos.y = ((minY * edge_size_miny) + (maxY * edge_size_maxy))/(edge_size_miny + edge_size_maxy);

    if ( not (charge > 0.) ) return nullptr;

    edm4hep::MutableTrackerHitPlane *recoHit = &(THcol->create())
    recoHit->setEDep((charge / m_electronsPerKeV) * dd4hep::keV);

    debug() << "Edge sizes, minx, maxx, miny, maxy: " << edge_size_minx << ", "
                                                      << edge_size_maxx << ", "
                                                      << edge_size_miny << ", "
                                                      << edge_size_maxy <<
            << "\nPosition: x = " << pos.x << " + " << m_layerHalfThickness[intState->currentLayer] * m_tanLorentzAngleX << "(LA-correction)";
    pos.x -= m_layerHalfThickness[intState->currentLayer] * m_tanLorentzAngleX;
    debug() << " = " << pos.x
            << "\n; y = " << pos.y << " + " << m_layerHalfThickness[intState->currentLayer] * m_tanLorentzAngleY << "(LA-correction)";
    pos.y -= m_layerHalfThickness[intState->currentLayer] * m_tanLorentzAngleY;
    debug() << " = " << pos.y;

    recoHit->setPosition(pos);
    recoHit->setdU( m_pixelSizeX / sqrt(12) );
    recoHit->setdV( m_pixelSizeY / sqrt(12) );
    time /= size;
    recoHit->setTime(time);
    debug() << "\ntime (ns) = " << time << endmsg;
          
    return recoHit;
}

/** Function transforms local coordinates in the ladder
 * into global coordinates
 */
void MuonCVXDDigitiser::TransformToLab(const int cellID, edm4hep::Vector3d xLoc, edm4hep::Vector3d xLab) const{
    // Use SurfaceManager to calculate global coordinates
    debug() << "Cell ID of Hit (used for transforming to lab coords)" << cellID << endmsg;
    SurfaceMap::const_iterator sI = m_map->find( cellID ) ;
    const dd4hep::rec::ISurface* surf = sI->second ;
    Vector2D oldPos( xLoc.x * dd4hep::mm, xLoc.y * dd4hep::mm );
    Vector3D lv = surf->localToGlobal( oldPos ) ;
    // Store local position in mm
    for ( int i = 0; i < 3; i++ )
      xLab[i] = lv[i] / dd4hep::mm;
}

/**
 * Function calculates position in pixel matrix based on the 
 * local coordinates of point in the ladder.
 */
void MuonCVXDDigitiser::TransformXYToCellID(double x, double y, int & ix, int & iy, InternalState *intState) const{
    int layer = intState->currentLayer;
    // Shift all of L/2 so that all numbers are positive
    if (isBarrel){
        double yInLadder = y + m_layerLadderLength[layer] / 2;
        iy = int(yInLadder / m_pixelSizeY);
        double xInLadder = x + m_layerLadderHalfWidth[layer];
        ix = int(xInLadder / m_pixelSizeX);
    } else {
        double yInPetal = y + m_layerPetalLength[layer] / 2;
        iy = int(yInPetal / m_pixelSizeY);
        //double localwidth = (m_layerPetalOuterWidth[layer] - m_layerPetalInnerWidth[layer])/(m_layerPetalLength[layer]) * yInPetal;
        double xInPetal = x + m_layerPetalOuterWidth[layer]/2;
        ix = int(xInPetal / m_pixelSizeX);
    }
}

/**
 Function calculates position in the local frame 
 based on the index of pixel in the ladder.
*/
void MuonCVXDDigitiser::TransformCellIDToXY(int ix, int iy, double & x, double & y, InternalState *intState) const{
    int layer = intState->currentLayer;
    // Put the point in the cell center
    if (isBarrel){
        y = ((0.5 + double(iy)) * m_pixelSizeY) - m_layerLadderLength[layer] / 2;
        x = ((0.5 + double(ix)) * m_pixelSizeX) - m_layerLadderHalfWidth[layer];
    } else {
        y = ((0.5 + double(iy)) * m_pixelSizeY) - m_layerPetalLength[layer] / 2;
        x = ((0.5 + double(ix)) * m_pixelSizeX) - m_layerPetalOuterWidth[layer] /2;
    }
}

int MuonCVXDDigitiser::GetPixelsInaColumn(InternalState *intState) const{//SP: why columns!?! I would have guess row..
    if (isBarrel){
        return ceil(m_layerLadderWidth[intState->currentLayer] / m_pixelSizeX);
    } else {
        return ceil(m_layerPetalOuterWidth[intState->currentLayer]/ m_pixelSizeX);
    }
}

int MuonCVXDDigitiser::GetPixelsInaRow(InternalState *intState) const{
    if (isBarrel){
        return ceil(m_layerLadderLength[intState->currentLayer] / m_pixelSizeY);
    } else {
        return ceil(m_layerPetalLength[intState->currentLayer] / m_pixelSizeY);
    }
}

void MuonCVXDDigitiser::PrintGeometryInfo() const{
    message() << "Number of layers: " << m_numberOfLayers
              << "\nPixel size X: " << m_pixelSizeX
              << "\nPixel size Y: " << m_pixelSizeY
              << "\nElectrons per KeV: " << m_electronsPerKeV
            //<< "\nSegment depth: " << m_segmentDepth;
    for (int i = 0; i < m_numberOfLayers; ++i) {
        message() << "Layer " << i
                  << "  Number of ladders: " << m_laddersInLayer[i]
                  << "  Radius: " << m_layerRadius[i]
                  << "  Ladder length: " << m_layerLadderLength[i]
                  << "  Ladder width: "<< m_layerLadderWidth[i]
                  << "  Ladder half width: " << m_layerLadderHalfWidth[i]
                  << "  Phi offset: " << m_layerPhiOffset[i]
                  << "  Active Si offset: " << m_layerActiveSiOffset[i]
                  << "  Half phi: " << m_layerHalfPhi[i]
                  << "  Thickness: " << m_layerThickness[i]
                  << "  Half thickness: " << m_layerHalfThickness[i]
                  //TODO: organize which get printed based on barrel/endcap
                  << "  Petal length: " << m_layerPetalLength[i]
                  << "  Petal Inner Width: " << m_layerPetalInnerWidth[i] 
                  << "  Petal Outer Width: " << m_layerPetalOuterWidth[i];
    }
    message() << endmsg;
}

//=============================================================================
// Sample charge from 1 / n^2 distribution.
//=============================================================================
double MuonCVXDDigitiser::randomTail( const double qmin, const double qmax ) const{
    const double offset = 1. / qmax;
    const double range  = ( 1. / qmin ) - offset;
    const double u      = offset + RandFlat::shoot() * range;
    return 1. / u;
}
