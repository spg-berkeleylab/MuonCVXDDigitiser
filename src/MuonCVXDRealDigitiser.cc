#include "MuonCVXDRealDigitiser.h"
#include <iostream>
#include <algorithm>
#include <math.h>

#include <EVENT/LCIO.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>

#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include "UTIL/LCTrackerConf.h"

#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>

#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DD4hepUnits.h"

#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_math.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h" 
#include "CLHEP/Random/RandFlat.h"

#include "DetElemSlidingWindow.h"
#include "TrivialSensor.h"
#include "HKBaseSensor.h"
    
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

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

MuonCVXDRealDigitiser aMuonCVXDRealDigitiser ;

MuonCVXDRealDigitiser::MuonCVXDRealDigitiser() :
    Processor("MuonCVXDRealDigitiser")
{}


void MuonCVXDRealDigitiser::init()
{ 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;
    _totEntries = 0;

    create_stats = stat_filename.compare(std::string { "None" }) != 0;

    if (create_stats)
    {
        double max_histox = std::max(_pixelSizeX, _pixelSizeY) * 10;
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
    }
}


void MuonCVXDRealDigitiser::processRunHeader(LCRunHeader* run)
{ 
    _nRun++ ;

    Detector& theDetector = Detector::getInstance();
    DetElement vxBarrel = theDetector.detector(_subDetName);              // TODO check missing barrel
    ZPlanarData&  zPlanarData = *vxBarrel.extension<ZPlanarData>();       // TODO check missing extension
    std::vector<ZPlanarData::LayerLayout> vx_layers = zPlanarData.layers;
    _numberOfLayers  = vx_layers.size();

    SurfaceManager& surfMan = *theDetector.extension<SurfaceManager>();
    _map = surfMan.map( vxBarrel.name() ) ;
    if( ! _map ) 
    {
      std::stringstream err  ; err << " Could not find surface map for detector: "
                                 << _subDetName << " in SurfaceManager " ;
      throw Exception( err.str() ) ;
    }

    _barrelID = vxBarrel.id();
    _laddersInLayer.resize(_numberOfLayers);
    _sensorsPerLadder.resize(_numberOfLayers);
    _layerHalfPhi.resize(_numberOfLayers);
    _layerHalfThickness.resize(_numberOfLayers);
    _layerThickness.resize(_numberOfLayers);
    _layerRadius.resize(_numberOfLayers);
    _layerLadderLength.resize(_numberOfLayers);
    _layerLadderWidth.resize(_numberOfLayers);
    _layerLadderHalfWidth.resize(_numberOfLayers);
    _layerActiveSiOffset.resize(_numberOfLayers);
    _layerPhiOffset.resize(_numberOfLayers);

    int curr_layer = 0;
    for(ZPlanarData::LayerLayout z_layout : vx_layers)
    {
        // ALE: Geometry is in cm, convert all lenght in mm
        _laddersInLayer[curr_layer] = z_layout.ladderNumber;

        _layerHalfPhi[curr_layer] = M_PI / ((double)_laddersInLayer[curr_layer]) ;

        _layerThickness[curr_layer] = z_layout.thicknessSensitive * dd4hep::cm / dd4hep::mm ;

        _layerHalfThickness[curr_layer] = 0.5 * _layerThickness[curr_layer];

        _layerRadius[curr_layer] = z_layout.distanceSensitive * dd4hep::cm / dd4hep::mm  + _layerHalfThickness[curr_layer];

        _sensorsPerLadder[curr_layer] = z_layout.sensorsPerLadder;

        _layerLadderLength[curr_layer] = z_layout.lengthSensor * z_layout.sensorsPerLadder * dd4hep::cm / dd4hep::mm ;

        _layerLadderWidth[curr_layer] = z_layout.widthSensitive * dd4hep::cm / dd4hep::mm ;

        _layerLadderHalfWidth[curr_layer] = _layerLadderWidth[curr_layer] / 2.;

        _layerActiveSiOffset[curr_layer] = - z_layout.offsetSensitive * dd4hep::cm / dd4hep::mm ;

        _layerPhiOffset[curr_layer] = z_layout.phi0;

        curr_layer++;
    }

    PrintGeometryInfo();
} 


void MuonCVXDRealDigitiser::processEvent(LCEvent * evt)
{ 
    LCCollectionVec *THcol = new LCCollectionVec(LCIO::TRACKERHITPLANE);
    CellIDEncoder<TrackerHitPlaneImpl> cellid_encoder(lcio::LCTrackerCellID::encoding_string(), THcol);

    LCCollectionVec* relCol = new LCCollectionVec(LCIO::LCRELATION);
    // to store the weights
    LCFlagImpl lcFlag { 0 };
    lcFlag.setBit(LCIO::LCREL_WEIGHTED);
    relCol->setFlag(lcFlag.getFlag());

    evt->addCollection(THcol, _outputCollectionName.c_str());
    evt->addCollection(relCol, _colVTXRelation.c_str());

    LCCollection* STHcol = nullptr;
    try
    {
        STHcol = evt->getCollection(_colName);
        streamlog_out( DEBUG9 ) << "Processing collection " << _colName  << " with "
                                <<  STHcol->getNumberOfElements()  << " hits ... " << std::endl ;
    }
    catch( lcio::DataNotAvailableException ex )
    {
        streamlog_out(ERROR) << _colName << " collection not available" << std::endl;
        STHcol = nullptr;
    }

    if (STHcol == nullptr or STHcol->getNumberOfElements() == 0)
    {
        streamlog_out(MESSAGE) << "Number of produced hits: " << THcol->getNumberOfElements()  << std::endl;
        return;
    }

    std::string encoder_str { STHcol->getParameters().getStringVal(lcio::LCIO::CellIDEncoding) };
    CellIDDecoder<TrackerHitPlaneImpl> cellid_decoder { encoder_str };

    std::size_t RELHISTOSIZE { 10 };
    vector<std::size_t> relHisto {};
    relHisto.assign(RELHISTOSIZE, 0);

    HitTemporalIndexes t_index { STHcol };

    for (int layer = 0; layer < _numberOfLayers; layer++)
    {
#pragma omp parallel for
        for (int ladder = 0; ladder < _laddersInLayer[layer]; ladder++)
        {
            int num_segment_x = 1;
            int nun_segment_y = _sensorsPerLadder[layer];

            float m_time = t_index.GetMinTime(layer, ladder);
            if (m_time == HitTemporalIndexes::MAXTIME)
            {
                if (streamlog::out.write<streamlog::DEBUG6>())
#pragma omp critical
                {
                    streamlog::out() << "Undefined min time for layer " << layer
                        << " ladder " << ladder << std::endl;
                }
                continue;
            }
            //clock time centered at 0
            float nw = floor(fabs(m_time) / _window_size);
            float start_time = (m_time >= 0) ? nw * _window_size : -1 * (nw + 1) * _window_size;

            AbstractSensor* sensor = nullptr;
            if (sensor_type == 1)
            {
                sensor = new TrivialSensor(layer, ladder, num_segment_x, nun_segment_y,
                                            _layerLadderLength[layer], _layerLadderWidth[layer],
                                            _layerThickness[layer], _pixelSizeX, _pixelSizeY,
                                            encoder_str, _barrelID, _threshold,
                                            start_time, _window_size);
            }
            else
            {
                sensor = new HKBaseSensor(layer, ladder, num_segment_x, nun_segment_y, 
                                            _layerLadderLength[layer], _layerLadderWidth[layer],
                                            _layerThickness[layer], _pixelSizeX, _pixelSizeY,
                                            encoder_str, _barrelID, _threshold, _fe_slope,
                                            start_time, _window_size);
            }

            if (sensor->GetStatus() != MatrixStatus::ok and streamlog::out.write<streamlog::ERROR>())
            {
                if (sensor->GetStatus() == MatrixStatus::pixel_number_error)
#pragma omp critical
                {
                    streamlog::out() << "Pixel number error for layer " << layer
                                        << " ladder " << ladder << std::endl;
                }
                else
#pragma omp critical
                {
                    streamlog::out() << "Segment number error for layer " << layer
                                        << " ladder " << ladder << std::endl;
                }
                continue;
            }

            DetElemSlidingWindow t_window {
                t_index, *sensor,
                _window_size, start_time,
                _tanLorentzAngleX, _tanLorentzAngleY,
                _cutOnDeltaRays,
                _diffusionCoefficient,
                _electronsPerKeV,
                _segmentLength,
                _energyLoss,
                3.0,
                _electronicNoise,
                _maxTrkLen,
                _deltaEne,
                _map
            };

            vector<std::size_t> histo_buffer {};

            while(t_window.active())
            {
                t_window.process();

                SegmentDigiHitList hit_buffer {};
                sensor->buildHits(hit_buffer);
                if (hit_buffer.size() == 0) continue;

                vector<TrackerHitPlaneImpl*> reco_buffer;
                reco_buffer.assign(hit_buffer.size(), nullptr);

                vector<LCRelationImpl*> rel_buffer;
                histo_buffer.assign(RELHISTOSIZE, 0);

                int idx = 0;
                for (SegmentDigiHit& digiHit : hit_buffer)
                {
                    TrackerHitPlaneImpl *recoHit = new TrackerHitPlaneImpl();
                    recoHit->setEDep((digiHit.charge / _electronsPerKeV) * dd4hep::keV);

                    bool sig = false;
                    double minx = 999;
                    double maxx = -999;
                    double miny = 999;
                    double maxy = -999;
                    double minz = 999;
                    double maxz = -999;

                    double loc_pos[3] = { 
                        digiHit.x - _layerHalfThickness[layer] * _tanLorentzAngleX,
                        digiHit.y - _layerHalfThickness[layer] * _tanLorentzAngleY,
                        0
                    };

                    recoHit->setCellID0(digiHit.cellID0);
                    recoHit->setCellID1(0);

                    SurfaceMap::const_iterator sI = _map->find(digiHit.cellID0);
                    const ISurface* surf = sI->second;

                    // See DetElemSlidingWindow::StoreSignalPoints
                    int segment_id = cellid_decoder(recoHit)["sensor"];
                    float s_offset = sensor->GetSensorCols() * sensor->GetPixelSizeY();
                    s_offset *= (float(segment_id) + 0.5);
                    s_offset -= sensor->GetHalfLength();

                    Vector2D oldPos(loc_pos[0] * dd4hep::mm, (loc_pos[1] - s_offset)* dd4hep::mm);
                    Vector3D lv = surf->localToGlobal(oldPos);

                    double xLab[3];
                    for ( int i = 0; i < 3; i++ )
                    {
                        xLab[i] = lv[i] / dd4hep::mm;
                    }

                    recoHit->setPosition(xLab);

                    recoHit->setTime(digiHit.time);

                    Vector3D u = surf->u() ;
                    Vector3D v = surf->v() ;

                    float u_direction[2] = { u.theta(), u.phi() };
                    float v_direction[2] = { v.theta(), v.phi() };

                    recoHit->setU( u_direction ) ;
                    recoHit->setV( v_direction ) ;

                    // ALE Does this make sense??? TO CHECK
                    recoHit->setdU( _pixelSizeX / sqrt(12) );
                    recoHit->setdV( _pixelSizeY / sqrt(12) );  

                    //All the sim-hits are registered for a given reco-hit
                    for (SimTrackerHit* st_item : digiHit.sim_hits)
                    {
                      if (create_stats)
                      {
                        if (!st_item->isOverlay()) sig = true;
                        if ( st_item->getPosition()[0] < minx ) minx = st_item->getPosition()[0];
                        else if ( st_item->getPosition()[0] > maxx ) maxx = st_item->getPosition()[0];
                        if ( st_item->getPosition()[1] < miny ) miny = st_item->getPosition()[1];
                        else if ( st_item->getPosition()[1] > maxy ) maxy = st_item->getPosition()[1];
                        if ( st_item->getPosition()[2] < minz ) minz = st_item->getPosition()[2];
                        else if ( st_item->getPosition()[2] > maxz ) maxz = st_item->getPosition()[2];
                      }
                        recoHit->rawHits().push_back( st_item );
                        LCRelationImpl* t_rel = new LCRelationImpl {};
                        t_rel->setFrom(recoHit);
                        t_rel->setTo(st_item);
                        t_rel->setWeight( 1.0 );
                        rel_buffer.push_back(t_rel);
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

