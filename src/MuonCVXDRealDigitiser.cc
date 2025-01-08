#include "MuonCVXDRealDigitiser.h"
#include <iostream>
#include <algorithm>
#include <math.h>

// edm4hep
#include <edm4hep/MCParticle.h>
#include "BitField64.hxx"
#include <GaudiKernel/ITHistSvc.h>

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

// Helpers
#include "DetElemSlidingWindow.h"
#include "TrivialSensor.h"
#include "HKBaseSensor.h"

// ROOT
#include <TFile.h>

using CLHEP::RandGauss;
using CLHEP::RandPoisson;
using CLHEP::RandFlat;

using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::rec::ZPlanarData;
using dd4hep::rec::SurfaceManager;
using dd4hep::rec::SurfaceMap;
using dd4hep::rec::ISurface;
using dd4hep::rec::Vector2D;
using dd4hep::rec::Vector3D;

DECLARE_COMPONENT(MuonCVXDRealDigitiser)

MuonCVXDRealDigitiser::MuonCVXDRealDigitiser(const std::string& name, ISvcLocator* svcLoc) : MultiTransformer(name, svcLoc,
          { KeyValues("CollectionName", {"VXDCollection"}) },
          { KeyValues("OutputCollectionName", {"VTXTrackerHits"}),
            KeyValues("RelationColName", {"VTXTrackerHitRelations"}) }) {}


StatusCode MuonCVXDRealDigitiser::initialize() {
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << "   init called  " << endmsg;



    if (m_create_stats) {
        ITHistSvc* histSvc{nullptr};
        StatusCode sc1 = service("THistSvc", histSvc);
        if ( sc1.isFailure() ) {
            log << MSG::ERROR << "Could not locate HistSvc" << endmsg;
            return StatusCode::FAILURE;
        }

        double max_histox = std::max(m_pixelSizeX, m_pixelSizeY) * 10;
        signal_dHisto = new TH1F("SignalHitDistance", "Signal Hit offset", 1000, 0., max_histox);
        bib_dHisto = new TH1F("BIBHitDistance", "BIB Hit offset", 1000, 0., max_histox);
        signal_cSizeHisto = new TH1F("SignalClusterSize", "Signal Cluster Size", 1000, 0., 50);
        signal_xSizeHisto = new TH1F("SignalClusterSizeinX", "Signal Cluster Size in x", 1000, 0., 20);
        signal_ySizeHisto = new TH1F("SignalClusterSizeinY", "Signal Cluster Size in y", 1000, 0., 20);
        signal_zSizeHisto = new TH1F("SignalClusterSizeinZ", "Signal Cluster Size in z", 1000, 0., 20);
        signal_eDepHisto = new TH1F("SignalClustereDep", "Signal Cluster Energy (MeV)", 1000, 0., 10e-1);
        bib_cSizeHisto = new TH1F("BIBClusterSize", "BIB Cluster Size", 1000, 0., 50);
        bib_xSizeHisto = new TH1F("BIBClusterSizeinX", "BIB Cluster Size in x", 1000, 0., 20);
        bib_ySizeHisto = new TH1F("BIBClusterSizeinY", "BIB Cluster Size in y", 1000, 0., 20);
        bib_zSizeHisto = new TH1F("BIBClusterSizeinZ", "BIB Cluster Size in z", 1000, 0., 20);
        bib_eDepHisto = new TH1F("BIBClusterSizeinZ", "BIB Cluster Energy (MeV)", 1000, 0., 10e-1);
 
        (void)histSvc->regHist("/histos/signal/hit_distance", signal_dHisto);
        (void)histSvc->regHist("/histos/signal/cluster_size", signal_cSizeHisto);
        (void)histSvc->regHist("/histos/signal/cluster_size_x", signal_xSizeHisto);
        (void)histSvc->regHist("/histos/signal/cluster_size_y", signal_ySizeHisto);
        (void)histSvc->regHist("/histos/signal/cluster_size_z", signal_zSizeHisto);
        (void)histSvc->regHist("/histos/signal/cluster_E_dep", signal_eDepHisto);

        (void)histSvc->regHist("/histos/bib/hit_distance", bib_dHisto);
        (void)histSvc->regHist("/histos/bib/cluster_size", bib_cSizeHisto);
        (void)histSvc->regHist("/histos/bib/cluster_size_x", bib_xSizeHisto);
        (void)histSvc->regHist("/histos/bib/cluster_size_y", bib_ySizeHisto);
        (void)histSvc->regHist("/histos/bib/cluster_size_z", bib_zSizeHisto);
        (void)histSvc->regHist("/histos/bib/cluster_E_dep", bib_eDepHisto);
    }

    return LoadGeometry();
}


StatusCode MuonCVXDRealDigitiser::LoadGeometry() { 
    Detector& theDetector = Detector::getInstance();
    DetElement vxBarrel = theDetector.detector(m_subDetName);              // TODO check missing barrel
    ZPlanarData&  zPlanarData = *vxBarrel.extension<ZPlanarData>();       // TODO check missing extension
    std::vector<ZPlanarData::LayerLayout> vx_layers = zPlanarData.layers;
    m_numberOfLayers = vx_layers.size();

    SurfaceManager& surfMan = *theDetector.extension<SurfaceManager>();
    m_map = surfMan.map( vxBarrel.name() );
    if( ! m_map ) 
    {
      MsgStream log(msgSvc(), name());
      log << MSG::ERROR << " Could not find surface map for detector: "
          << m_subDetName << " in SurfaceManager " << endmsg;
      return StatusCode::FAILURE;
    }

    m_barrelID = vxBarrel.id();
    m_laddersInLayer.resize(m_numberOfLayers);
    m_sensorsPerLadder.resize(m_numberOfLayers);
    m_layerHalfPhi.resize(m_numberOfLayers);
    m_layerHalfThickness.resize(m_numberOfLayers);
    m_layerThickness.resize(m_numberOfLayers);
    m_layerRadius.resize(m_numberOfLayers);
    m_layerLadderLength.resize(m_numberOfLayers);
    m_layerLadderWidth.resize(m_numberOfLayers);
    m_layerLadderHalfWidth.resize(m_numberOfLayers);
    m_layerActiveSiOffset.resize(m_numberOfLayers);
    m_layerPhiOffset.resize(m_numberOfLayers);

    int curr_layer = 0;
    for(ZPlanarData::LayerLayout z_layout : vx_layers) {
        // ALE: Geometry is in cm, convert all lenght in mm
        m_laddersInLayer[curr_layer] = z_layout.ladderNumber;

        m_layerHalfPhi[curr_layer] = M_PI / ((double)m_laddersInLayer[curr_layer]) ;

        m_layerThickness[curr_layer] = z_layout.thicknessSensitive * dd4hep::cm / dd4hep::mm ;

        m_layerHalfThickness[curr_layer] = 0.5 * m_layerThickness[curr_layer];

        m_layerRadius[curr_layer] = z_layout.distanceSensitive * dd4hep::cm / dd4hep::mm  + m_layerHalfThickness[curr_layer];

        m_sensorsPerLadder[curr_layer] = z_layout.sensorsPerLadder;

        m_layerLadderLength[curr_layer] = z_layout.lengthSensor * z_layout.sensorsPerLadder * dd4hep::cm / dd4hep::mm ;

        m_layerLadderWidth[curr_layer] = z_layout.widthSensitive * dd4hep::cm / dd4hep::mm ;

        m_layerLadderHalfWidth[curr_layer] = m_layerLadderWidth[curr_layer] / 2.;

        m_layerActiveSiOffset[curr_layer] = - z_layout.offsetSensitive * dd4hep::cm / dd4hep::mm ;

        m_layerPhiOffset[curr_layer] = z_layout.phi0;

        curr_layer++;
    }

    PrintGeometryInfo();
    return StatusCode::SUCCESS;
} 


std::tuple<edm4hep::TrackerHitPlaneCollection,
           edm4hep::TrackerHitSimTrackerHitLinkCollection> MuonCVXDRealDigitiser::operator()(
     const edm4hep::SimTrackerHitCollection& STHcol) const{
    MsgStream log(msgSvc(), name());

    edm4hep::TrackerHitPlaneCollection             THcol;
    edm4hep::TrackerHitSimTrackerHitLinkCollection relCol;

    std::string encoder_str = "subdet:5,side:-2,layer:9,module:8,sensor:8";
    BitField64 cellID_coder(encoder_str);

    if (STHcol.size() == 0) {
        log << MSG::INFO << "Number of produced hits: " << THcol.size()  << endmsg;
        return std::make_tuple( std::move(THcol), std::move(relCol) );
    }

    std::size_t RELHISTOSIZE { 10 };
    vector<std::size_t> relHisto {};
    relHisto.assign(RELHISTOSIZE, 0);

    HitTemporalIndexes t_index { STHcol };

    for (int layer = 0; layer < m_numberOfLayers; layer++) {
#pragma omp parallel for
        for (int ladder = 0; ladder < m_laddersInLayer[layer]; ladder++) {
            int num_segment_x = 1;
            int nun_segment_y = m_sensorsPerLadder[layer];

            float m_time = t_index.GetMinTime(layer, ladder);
            if (m_time == HitTemporalIndexes::MAXTIME) {
                if ( msgLevel(MSG::DEBUG) )
#pragma omp critical
                {
                    log << MSG::DEBUG << "Undefined min time for layer " << layer
                        << " ladder " << ladder << endmsg;
                }
                continue;
            }
            //clock time centered at 0
            float nw = floor(fabs(m_time) / m_window_size);
            float start_time = (m_time >= 0) ? nw * m_window_size : -1 * (nw + 1) * m_window_size;

            AbstractSensor* sensor = nullptr;
            if (m_sensor_type == 1) {
                sensor = new TrivialSensor(layer, ladder, num_segment_x, nun_segment_y,
                                          m_layerLadderLength[layer], m_layerLadderWidth[layer],
                                          m_layerThickness[layer], m_pixelSizeX, m_pixelSizeY,
                                          encoder_str, m_barrelID, m_threshold,
                                          start_time, m_window_size);
            } else {
                sensor = new HKBaseSensor(layer, ladder, num_segment_x, nun_segment_y, 
                                          m_layerLadderLength[layer], m_layerLadderWidth[layer],
                                          m_layerThickness[layer], m_pixelSizeX, m_pixelSizeY,
                                          encoder_str, m_barrelID, m_threshold, m_fe_slope,
                                          start_time, m_window_size);
            }

            if (sensor->GetStatus() != MatrixStatus::ok and msgLevel(MSG::ERROR)) {
                if (sensor->GetStatus() == MatrixStatus::pixel_number_error)
#pragma omp critical
                {
                    log << MSG::ERROR << "Pixel number error for layer " << layer
                        << " ladder " << ladder << endmsg;
                } else
#pragma omp critical
                {
                    log << MSG::ERROR << "Segment number error for layer " << layer
                        << " ladder " << ladder << endmsg;
                }
                continue;
            }

            DetElemSlidingWindow t_window {
                t_index, *sensor,
                m_window_size, start_time,
                m_tanLorentzAngleX, m_tanLorentzAngleY,
                m_cutOnDeltaRays,
                m_diffusionCoefficient,
                m_electronsPerKeV,
                m_segmentLength,
                m_energyLoss,
                3.0,
                m_electronicNoise,
                m_maxTrkLen,
                m_deltaEne,
                m_map
            };

            vector<std::size_t> histo_buffer {};

            while(t_window.active()) {
                t_window.process(msgSvc());

                SegmentDigiHitList hit_buffer {};
                sensor->buildHits(hit_buffer, service<IMessageSvc>("MessageSvc"));
                if (hit_buffer.size() == 0) continue;

                std::vector<edm4hep::MutableTrackerHitPlane*> reco_buffer;
                reco_buffer.assign(hit_buffer.size(), nullptr);

                std::vector<edm4hep::MutableTrackerHitSimTrackerHitLink*> rel_buffer;
                histo_buffer.assign(RELHISTOSIZE, 0);

                int idx = 0;
                for (SegmentDigiHit& digiHit : hit_buffer) {
                    edm4hep::MutableTrackerHitPlane *recoHit = new edm4hep::MutableTrackerHitPlane();
                    recoHit->setEDep((digiHit.charge / m_electronsPerKeV) * dd4hep::keV);

                    bool sig = false;
                    double minx = 999;
                    double maxx = -999;
                    double miny = 999;
                    double maxy = -999;
                    double minz = 999;
                    double maxz = -999;

                    edm4hep::Vector3d loc_pos(
                        digiHit.x - m_layerHalfThickness[layer] * m_tanLorentzAngleX,
                        digiHit.y - m_layerHalfThickness[layer] * m_tanLorentzAngleY,
                        0
                    );

                    recoHit->setCellID(digiHit.cellID);

                    SurfaceMap::const_iterator sI = m_map->find(digiHit.cellID);
                    const ISurface* surf = sI->second;

                    // See DetElemSlidingWindow::StoreSignalPoints
                    cellID_coder.setValue(recoHit->getCellID());
                    int segment_id = cellID_coder["sensor"];
                    float s_offset = sensor->GetSensorCols() * sensor->GetPixelSizeY();
                    s_offset *= (float(segment_id) + 0.5);
                    s_offset -= sensor->GetHalfLength();

                    Vector2D oldPos(loc_pos.x * dd4hep::mm, (loc_pos.y - s_offset)* dd4hep::mm);
                    Vector3D lv = surf->localToGlobal(oldPos);

                    edm4hep::Vector3d xLab;
                    xLab.x = lv[0] / dd4hep::mm;
                    xLab.y = lv[1] / dd4hep::mm;
                    xLab.z = lv[2] / dd4hep::mm;

                    recoHit->setPosition(xLab);

                    recoHit->setTime(digiHit.time);

                    Vector3D u = surf->u() ;
                    Vector3D v = surf->v() ;

                    edm4hep::Vector2f u_direction( u.theta(), u.phi() );
                    edm4hep::Vector2f v_direction( v.theta(), v.phi() );

                    recoHit->setU( u_direction ) ;
                    recoHit->setV( v_direction ) ;

                    // ALE Does this make sense??? TO CHECK
                    recoHit->setDu( m_pixelSizeX / sqrt(12) );
                    recoHit->setDv( m_pixelSizeY / sqrt(12) );  

                    //All the sim-hits are registered for a given reco-hit
                    for (edm4hep::SimTrackerHit *st_item : digiHit.sim_hits) {
                      if (m_create_stats) {
                        if (!st_item->isOverlay()) sig = true;
                        if ( st_item->getPosition().x < minx ) minx = st_item->getPosition().x;
                        else if ( st_item->getPosition().x > maxx ) maxx = st_item->getPosition().x;
                        if ( st_item->getPosition().y < miny ) miny = st_item->getPosition().y;
                        else if ( st_item->getPosition().y > maxy ) maxy = st_item->getPosition().y;
                        if ( st_item->getPosition().z < minz ) minz = st_item->getPosition().z;
                        else if ( st_item->getPosition().z > maxz ) maxz = st_item->getPosition().z;
                      }
                        //recoHit->rawHits().push_back( st_item );
                        edm4hep::MutableTrackerHitSimTrackerHitLink* t_rel = new edm4hep::MutableTrackerHitSimTrackerHitLink();
                        t_rel->setFrom( *recoHit );
                        t_rel->setTo( *st_item );
                        t_rel->setWeight( 1.0 );
                        rel_buffer.push_back(t_rel);
                    }

                    if (digiHit.sim_hits.size() < RELHISTOSIZE) {
                        histo_buffer[digiHit.sim_hits.size()]++;
                    }

                    reco_buffer[idx] = recoHit;
                    idx++;

                    if (m_create_stats) {
                      // cluster size histograms
                      if ( !sig ) {
                        //bib_cSizeHisto->Fill(recoHit->getRawHits().size());
                        bib_cSizeHisto->Fill(digiHit.size);
                        bib_xSizeHisto->Fill(maxx-minx);
                        bib_ySizeHisto->Fill(maxy-miny);
                        bib_zSizeHisto->Fill(maxz-minz);
                        bib_eDepHisto->Fill(1000*recoHit->getEDep());
                      } else {
                        //signal_cSizeHisto->Fill(recoHit->getRawHits().size());
                        signal_cSizeHisto->Fill(digiHit.size);
                        signal_xSizeHisto->Fill(maxx-minx);
                        signal_ySizeHisto->Fill(maxy-miny);
                        signal_zSizeHisto->Fill(maxz-minz);
                        signal_eDepHisto->Fill(1000*recoHit->getEDep());
                      }
                   }
                }

                if (reco_buffer.size() > 0)
#pragma omp critical               
                {
                    for(edm4hep::MutableTrackerHitPlane* recoHit : reco_buffer) {
                        if (recoHit == nullptr) continue;
                        THcol.push_back( *recoHit );

                        if ( msgLevel(MSG::DEBUG) ) {
                            cellID_coder.setValue(recoHit->getCellID());
                            log << MSG::DEBUG << "Reconstructed pixel cluster for " 
                                << sensor->GetLayer() << ":" << sensor->GetLadder() 
                                << ":" << cellID_coder["sensor"] << std::endl
                                << "- global position (x,y,z,t) = " 
                                << recoHit->getPosition().x 
                                << ", " << recoHit->getPosition().y
                                << ", " << recoHit->getPosition().z
                                << ", " << recoHit->getTime() << std::endl
                                << "- charge = " << recoHit->getEDep() << endmsg;
                        }
                    }

                    for (edm4hep::MutableTrackerHitSimTrackerHitLink* rel_item : rel_buffer) {
                        relCol.push_back( *rel_item );
                    }
                    for (std::size_t k = 0; k < RELHISTOSIZE; k++) {
                        relHisto[k] += histo_buffer[k];
                    }
                }
            }

            delete sensor;
        }
    }
    log << MSG::INFO << "Number of produced hits: " << THcol.size()  << endmsg;
    int count = 0;
    log << MSG::DEBUG << "Hit relation histogram:" << std::endl;
    for (std::size_t k = 0; k < RELHISTOSIZE; k++) {
        log << MSG::DEBUG << k << " " << relHisto[k] << std::endl;
        count += relHisto[k];
    }
    log << MSG::DEBUG << "> " << THcol.size() - count << endmsg;

    if (m_create_stats) {
        for (int i = 0; i < relCol.size(); ++i) {
            edm4hep::TrackerHitSimTrackerHitLink hitRel = relCol.at(i);

            float tmpf = 0.0;
            tmpf += pow(hitRel.getFrom().getPosition().x - hitRel.getTo().getPosition().x, 2);
            tmpf += pow(hitRel.getFrom().getPosition().y - hitRel.getTo().getPosition().y, 2);
            tmpf += pow(hitRel.getFrom().getPosition().z - hitRel.getTo().getPosition().z, 2);
            
            if (hitRel.getTo().isOverlay()) {
                bib_dHisto->Fill(sqrt(tmpf));
            } else {
                signal_dHisto->Fill(sqrt(tmpf));
            }
        }
    }

    return std::make_tuple( std::move(THcol), std::move(relCol) );
}

StatusCode MuonCVXDRealDigitiser::finalize() {
  return StatusCode::SUCCESS;
}

void MuonCVXDRealDigitiser::PrintGeometryInfo() {
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "Number of layers: " << m_numberOfLayers << std::endl
        << "Pixel size X: " << m_pixelSizeX << std::endl
        << "Pixel size Y: " << m_pixelSizeY << std::endl
        << "Electrons per KeV: " << m_electronsPerKeV << std::endl;
    for (int i = 0; i < m_numberOfLayers; ++i) {
        log << MSG::INFO << "Layer " << i << std::endl
            << "  Number of ladders: " << m_laddersInLayer[i] << std::endl
            << "  Radius: " << m_layerRadius[i] << std::endl
            << "  Ladder length: " << m_layerLadderLength[i] << std::endl
            << "  Ladder width: "<< m_layerLadderWidth[i] << std::endl
            << "  Ladder half width: " << m_layerLadderHalfWidth[i] << std::endl
            << "  Phi offset: " << m_layerPhiOffset[i] << std::endl
            << "  Active Si offset: " << m_layerActiveSiOffset[i] << std::endl
            << "  Half phi: " << m_layerHalfPhi[i] << std::endl
            << "  Thickness: " << m_layerThickness[i] << std::endl
            << "  Half thickness: " << m_layerHalfThickness[i] << std::endl;
    }
    log << MSG::INFO << endmsg;
}

