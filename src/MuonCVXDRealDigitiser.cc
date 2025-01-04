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
          { KeyValues("SimHitLocCollectionName", {"VertexBarrel"}),
            KeyValues("OutputCollectionName", {"VTXTrackerHits"}),
            KeyValues("RelationColName", {"VTXTrackerHitRelations"}),
            KeyValues("RawHitsLinkColName", {"VTXRawHitRelations"}) }) {}


StatusCode MuonCVXDRealDigitiser::initialize() {
    debug() << "   init called  " << endmsg;

    create_stats = stat_filename.compare(std::string { "None" }) != 0;

    if (create_stats) {
        ITHistSvc* histSvc{nullptr};
        StatusCode sc1 = service("THistSvc", histSvc);
        if ( sc1.isFailure() ) {
            error() << "Could not locate HistSvc" << endmsg;
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

    return StatusCode::SUCCESS;
}


void MuonCVXDRealDigitiser::LoadGeometry() const{ 
    Detector& theDetector = Detector::getInstance();
    DetElement vxBarrel = theDetector.detector(m_subDetName);              // TODO check missing barrel
    ZPlanarData&  zPlanarData = *vxBarrel.extension<ZPlanarData>();       // TODO check missing extension
    std::vector<ZPlanarData::LayerLayout> vx_layers = zPlanarData.layers;
    m_numberOfLayers  = vx_layers.size();

    SurfaceManager& surfMan = *theDetector.extension<SurfaceManager>();
    m_map = surfMan.map( vxBarrel.name() );
    if( ! m_map ) 
    {
      std::stringstream err;
      err << " Could not find surface map for detector: "
          << m_subDetName << " in SurfaceManager ";
      throw Exception( err.str() ) ;
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
} 


std::tuple<edm4hep::SimTrackerHitCollection,
                       edm4hep::TrackerHitPlaneCollection,
                       edm4hep::TrackerHitSimTrackerHitCollection,
                       edm4hep::TrackerHitSimTrackerHitCollection> operator(
                 const edm4hep::SimTrackerHitCollection& STHcol) const{
    if ( !m_map ) {
        LoadGeometry();
    }
    edm4hep::SimTrackerHitCollection               STHLocCol;
    edm4hep::TrackerHitPlaneCollection             THcol;
    edm4hep::TrackerHitSimTrackerHitLinkCollection relCol;
    edm4hep::TrackerHitSimTrackerHitLinkCollection rawHitsCol;

    BitField64 cellID_coder("subdet:5,side:-2,layer:9,module:8,sensor:8");

    if (STHcol == nullptr or STHcol.size() == 0) {
        message() << "Number of produced hits: " << THcol.size()  << endmsg;
        return std::make_tuple(std::move(STHLocCol),
                               std::move(THcol),
                               std::move(relCol),
                               std::move(rawHitsCol) );
    }

    std::size_t RELHISTOSIZE { 10 };
    vector<std::size_t> relHisto {};
    relHisto.assign(RELHISTOSIZE, 0);

    HitTemporalIndexes t_index { STHcol };

    for (int layer = 0; layer < _numberOfLayers; layer++)
    {
#pragma omp parallel for
        for (int ladder = 0; ladder < m_laddersInLayer[layer]; ladder++)
        {
            int num_segment_x = 1;
            int nun_segment_y = m_sensorsPerLadder[layer];

            float m_time = t_index.GetMinTime(layer, ladder);
            if (m_time == HitTemporalIndexes::MAXTIME)
            {
                if ( msgLevel(MSG::DEBUG) )
#pragma omp critical
                {
                    debug() << "Undefined min time for layer " << layer
                            << " ladder " << ladder << endmsg;
                }
                continue;
            }
            //clock time centered at 0
            float nw = floor(fabs(m_time) / m_window_size);
            float start_time = (m_time >= 0) ? nw * m_window_size : -1 * (nw + 1) * m_window_size;

            AbstractSensor* sensor = nullptr;
            if (sensor_type == 1) {
                sensor = new TrivialSensor(layer, ladder, num_segment_x, nun_segment_y,
                                          m_layerLadderLength[layer], m_layerLadderWidth[layer],
                                          m_layerThickness[layer], m_pixelSizeX, m_pixelSizeY,
                                          encoder_str, _barrelID, m_threshold,
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
                    error() << "Pixel number error for layer " << layer
                            << " ladder " << ladder << endmsg;
                } else
#pragma omp critical
                {
                    error() << "Segment number error for layer " << layer
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
                t_window.process();

                SegmentDigiHitList hit_buffer {};
                sensor->buildHits(hit_buffer);
                if (hit_buffer.size() == 0) continue;

                std::vector<edm4hep::MutableTrackerHitPlane*> reco_buffer;
                reco_buffer.assign(hit_buffer.size(), nullptr);

                std::vector<edm4hep::TrackerHitSimTrackerHitLink*> rel_buffer;
                histo_buffer.assign(RELHISTOSIZE, 0);

                std::vector<edm4hep::SimTrackerHit*> rawHits_buffer;
                std::ccvector<edm4hep::TrackerHitSimTrackerHitLink*> rawHitsLink_buffer;

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
                    cellID_coder.setValue(recoHit.getcellID())
                    int segment_id = cellID_coder["sensor"];
                    float s_offset = sensor->GetSensorCols() * sensor->GetPixelSizeY();
                    s_offset *= (float(segment_id) + 0.5);
                    s_offset -= sensor->GetHalfLength();

                    Vector2D oldPos(loc_pos.x * dd4hep::mm, (loc_pos.y - s_offset)* dd4hep::mm);
                    Vector3D lv = surf->localToGlobal(oldPos);

                    edm4hep::Vector3d xLab;
                    for ( int i = 0; i < 3; i++ )
                    {
                        xLab[i] = lv[i] / dd4hep::mm;
                    }

                    recoHit->setPosition(xLab);

                    recoHit->setTime(digiHit.time);

                    Vector3D u = surf->u() ;
                    Vector3D v = surf->v() ;

                    edm4hep::Vector2f u_direction( u.theta(), u.phi() );
                    edm4hep::Vector2f v_direction( v.theta(), v.phi() );

                    recoHit->setU( u_direction ) ;
                    recoHit->setV( v_direction ) ;

                    // ALE Does this make sense??? TO CHECK
                    recoHit->setdU( m_pixelSizeX / sqrt(12) );
                    recoHit->setdV( m_pixelSizeY / sqrt(12) );  

                    //All the sim-hits are registered for a given reco-hit
                    for (edm4hep::SimTrackerHit *st_item : digiHit.sim_hits) {
                      if (create_stats) {
                        if (!st_item->isOverlay()) sig = true;
                        if ( st_item->getPosition().x < minx ) minx = st_item->getPosition().x;
                        else if ( st_item->getPosition().x > maxx ) maxx = st_item->getPosition().x;
                        if ( st_item->getPosition().y < miny ) miny = st_item->getPosition().y;
                        else if ( st_item->getPosition().y > maxy ) maxy = st_item->getPosition().y;
                        if ( st_item->getPosition().z < minz ) minz = st_item->getPosition().z;
                        else if ( st_item->getPosition().z > maxz ) maxz = st_item->getPosition().z;
                      }
                        //recoHit->rawHits().push_back( st_item );
                        rawHits_buffer.push_back( st_item );
                        edm4hep::MutableTrackerHitSimTrackerHitLink* t_rel = new edm4hep::MutableTrackerHitSimTrackerHitLink();
                        t_rel->setFrom(recoHit);
                        t_rel->setTo(st_item);
                        t_rel->setWeight( 1.0 );
                        rel_buffer.push_back(t_rel);

                        edm4hep::MutableTrackerHitSimTrackerHitLink* rawLink = new edm4hep::MutableTrackerHitSimTrackerHitLink();
                        rawLink->setFrom(recoHit);
                        rawLink->setTo(st_item)
                    }

                    if (digiHit.sim_hits.size() < RELHISTOSIZE)
                    {
                        histo_buffer[digiHit.sim_hits.size()]++;
                    }

                    reco_buffer[idx] = recoHit;
                    idx++;

                    if (create_stats)
                    {
                      // cluster size histograms
                      if ( !sig )
                      {
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
                    for(TrackerHitPlaneImpl* recoHit : reco_buffer)
                    {
                        if (recoHit == nullptr) continue;
                        THcol->addElement(recoHit);

                        if (streamlog::out.write<streamlog::DEBUG7>())
                        {
                            streamlog::out() << "Reconstructed pixel cluster for " 
                                             << sensor->GetLayer() << ":" << sensor->GetLadder() 
                                             << ":" << cellid_decoder(recoHit)["sensor"] << std::endl
                                             << "- global position (x,y,z,t) = " << recoHit->getPosition()[0] 
                                             << ", " << recoHit->getPosition()[1] 
                                             << ", " << recoHit->getPosition()[2] 
                                             << ", " << recoHit->getTime() << std::endl
                                             << "- charge = " << recoHit->getEDep() << std::endl;
                        }
                    }

                    for (LCRelationImpl* rel_item : rel_buffer)
                    {
                        relCol->addElement(rel_item);
                    }
                    for (std::size_t k = 0; k < RELHISTOSIZE; k++)
                    {
                        relHisto[k] += histo_buffer[k];
                    }
                }
            }

            delete sensor;
        }
    }
    streamlog_out(MESSAGE) << "Number of produced hits: " << THcol->getNumberOfElements()  << std::endl;
    int count = 0;
      streamlog_out(DEBUG) << "Hit relation histogram:" << std::endl;
      for (std::size_t k = 0; k < RELHISTOSIZE; k++)
      {
        streamlog_out(DEBUG) << k << " " << relHisto[k] << std::endl;
        count += relHisto[k];
      }
      streamlog_out(DEBUG) << "> " << THcol->getNumberOfElements() - count << std::endl;

    if (create_stats)
    {
        for (int i = 0; i < relCol->getNumberOfElements(); ++i)
        {
            LCRelationImpl* hitRel = dynamic_cast<LCRelationImpl*>(relCol->getElementAt(i));
            TrackerHitPlaneImpl* recoHit = dynamic_cast<TrackerHitPlaneImpl*>(hitRel->getFrom());
            SimTrackerHit* simTrkHit = dynamic_cast<SimTrackerHit*>(hitRel->getTo());

            float tmpf = pow(recoHit->getPosition()[0] - simTrkHit->getPosition()[0], 2);
            tmpf += pow(recoHit->getPosition()[1] - simTrkHit->getPosition()[1], 2);
            tmpf += pow(recoHit->getPosition()[2] - simTrkHit->getPosition()[2], 2);
            
            if (simTrkHit->isOverlay())
            {
                bib_dHisto->Fill(sqrt(tmpf));
            }
            else
            {
                signal_dHisto->Fill(sqrt(tmpf));
            }
        }
    }
}

void MuonCVXDRealDigitiser::check(LCEvent *evt)
{}

void MuonCVXDRealDigitiser::end()
{
    streamlog_out(DEBUG) << "   end called  " << std::endl;

    if (create_stats)
    {
        TFile statFile = TFile(stat_filename.c_str(), "recreate");
        statFile.WriteObject(signal_dHisto, "Signal offset");
        statFile.WriteObject(bib_dHisto, "BIB offset");
        statFile.WriteObject(signal_cSizeHisto, "Signal cluster size");
        statFile.WriteObject(signal_xSizeHisto, "Signal cluster size in x");
        statFile.WriteObject(signal_ySizeHisto, "Signal cluster size in y");
        statFile.WriteObject(signal_zSizeHisto, "Signal cluster size in z");
        statFile.WriteObject(signal_eDepHisto, "Signal energy");
        statFile.WriteObject(bib_cSizeHisto, "BIB cluster size");
        statFile.WriteObject(bib_xSizeHisto, "BIB cluster size in x");
        statFile.WriteObject(bib_ySizeHisto, "BIB cluster size in y");
        statFile.WriteObject(bib_zSizeHisto, "BIB cluster size in z");
        statFile.WriteObject(bib_eDepHisto, "BIB energy");
        statFile.Flush();
        statFile.Close();
    }

    if (signal_dHisto != nullptr) delete(signal_dHisto);
    if (bib_dHisto != nullptr) delete(bib_dHisto);
    if (signal_cSizeHisto != nullptr) delete(signal_cSizeHisto);
    if (signal_xSizeHisto != nullptr) delete(signal_xSizeHisto);
    if (signal_ySizeHisto != nullptr) delete(signal_ySizeHisto);
    if (signal_zSizeHisto != nullptr) delete(signal_zSizeHisto);
    if (signal_eDepHisto != nullptr) delete(signal_eDepHisto);
    if (bib_cSizeHisto != nullptr) delete(bib_cSizeHisto);
    if (bib_xSizeHisto != nullptr) delete(bib_xSizeHisto);
    if (bib_ySizeHisto != nullptr) delete(bib_ySizeHisto);
    if (bib_zSizeHisto != nullptr) delete(bib_zSizeHisto);
    if (bib_eDepHisto != nullptr) delete(bib_eDepHisto);

}

void MuonCVXDRealDigitiser::PrintGeometryInfo()
{
    streamlog_out(MESSAGE) << "Number of layers: " << _numberOfLayers << std::endl;
    streamlog_out(MESSAGE) << "Pixel size X: " << _pixelSizeX << std::endl;
    streamlog_out(MESSAGE) << "Pixel size Y: " << _pixelSizeY << std::endl;
    streamlog_out(MESSAGE) << "Electrons per KeV: " << _electronsPerKeV << std::endl;
    for (int i = 0; i < _numberOfLayers; ++i) 
    {
        streamlog_out(MESSAGE) << "Layer " << i << std::endl;
        streamlog_out(MESSAGE) << "  Number of ladders: " << _laddersInLayer[i] << std::endl;
        streamlog_out(MESSAGE) << "  Radius: " << _layerRadius[i] << std::endl;
        streamlog_out(MESSAGE) << "  Ladder length: " << _layerLadderLength[i] << std::endl;
        streamlog_out(MESSAGE) << "  Ladder width: "<< _layerLadderWidth[i] << std::endl;
        streamlog_out(MESSAGE) << "  Ladder half width: " << _layerLadderHalfWidth[i] << std::endl;
        streamlog_out(MESSAGE) << "  Phi offset: " << _layerPhiOffset[i] << std::endl;
        streamlog_out(MESSAGE) << "  Active Si offset: " << _layerActiveSiOffset[i] << std::endl;
        streamlog_out(MESSAGE) << "  Half phi: " << _layerHalfPhi[i] << std::endl;
        streamlog_out(MESSAGE) << "  Thickness: " << _layerThickness[i] << std::endl;
        streamlog_out(MESSAGE) << "  Half thickness: " << _layerHalfThickness[i] << std::endl;
    }
}

