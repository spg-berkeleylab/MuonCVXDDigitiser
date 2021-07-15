#include "MuonCVXDRealDigitiser.h"
#include <iostream>
#include <algorithm>

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
#include "HKBaseSensor.h"
    
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

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
    Processor("MuonCVXDRealDigitiser"),
    _nRun(0),
    _nEvt(0),
    _totEntries(0),
    _barrelID(0)
{
    _description = "MuonCVXDRealDigitiser should create VTX TrackerHits from SimTrackerHits";

    registerInputCollection(LCIO::SIMTRACKERHIT,
                           "CollectionName", 
                           "Name of the SimTrackerHit collection",
                           _colName,
                           std::string("VXDCollection"));

    registerOutputCollection(LCIO::TRACKERHITPLANE,
                            "OutputCollectionName", 
                            "Name of the output TrackerHit collection",
                            _outputCollectionName,
                            std::string("VTXTrackerHits"));

    registerOutputCollection(LCIO::LCRELATION,
                            "RelationColName", 
                            "Name of the output VTX trackerhit relation collection",
                            _colVTXRelation,
                            std::string("VTXTrackerHitRelations"));

    registerProcessorParameter("SubDetectorName", 
                               "Name of Vertex detector",
                               _subDetName,
                               std::string("VertexBarrel"));

    registerProcessorParameter("TanLorentz",
                               "Tangent of Lorentz Angle",
                               _tanLorentzAngleX,
                               (double)0.8);

    registerProcessorParameter("TanLorentzY",
                               "Tangent of Lorentz Angle along Y",
                               _tanLorentzAngleY,
                               (double)0);

    registerProcessorParameter("CutOnDeltaRays",
                               "Cut on delta-ray energy (MeV)",
                               _cutOnDeltaRays,
                               (double)0.030);

    // For diffusion-coeffieint calculation, see e.g. https://www.slac.stanford.edu/econf/C060717/papers/L008.PDF
    // or directly Eq. 13 of https://cds.cern.ch/record/2161627/files/ieee-tns-07272141.pdf
    // diffusionCoefficient = sqrt(2*D / mu / V), where
    //  - D = 12 cm^2/s // diffusion constant
    //  - mu = 450 cm^2/s/V // mobility
    //  - V = 10-30 V // expected depletion voltage
    //  => _diffusionCoefficient = 0.04-0.07
    registerProcessorParameter("DiffusionCoefficient",
                               "Diffusion coefficient, sqrt(D / mu / V).",
                               _diffusionCoefficient,
                               (double)0.07);

    registerProcessorParameter("PixelSizeX",
                               "Pixel Size X",
                               _pixelSizeX,
                               (double)0.025);

    registerProcessorParameter("PixelSizeY",
                               "Pixel Size Y",
                               _pixelSizeY,
                               (double)0.025);

    registerProcessorParameter("ElectronsPerKeV",
                               "Electrons per keV",
                               _electronsPerKeV,
                               (double)270.3);

    registerProcessorParameter("Threshold",
                               "Cell Threshold in electrons",
                               _threshold,
                               200.);

    registerProcessorParameter("SegmentLength",
                               "Segment Length in mm",
                               _segmentLength,
                               double(0.005));

    registerProcessorParameter("PoissonSmearing",
                               "Apply Poisson smearing of electrons collected on pixels",
                               _PoissonSmearing,
                               1);

    registerProcessorParameter("ElectronicEffects",
                               "Apply Electronic Effects",
                               _electronicEffects,
                               int(1));

    registerProcessorParameter("ElectronicNoise",
                               "electronic noise in electrons",
                               _electronicNoise,
                               100.);

    registerProcessorParameter("StoreFiredPixels",
                               "Store fired pixels",
                               _produceFullPattern,
                               int(0));

    registerProcessorParameter("EnergyLoss",
                               "Energy Loss keV/mm",
                               _energyLoss,
                               double(280.0));
                               
    registerProcessorParameter("MaxEnergyDelta",
                               "Max delta in energy between G4 prediction and random sampling for each hit in electrons",
                               _deltaEne,
                               100.0);                               

    registerProcessorParameter("MaxTrackLength",
                               "Maximum values for track length (in mm)",
                               _maxTrkLen,
                               10.0); 

    registerProcessorParameter("WindowSize",
                               "Window size",
                               _window_size,
                               (float)0.1);

    registerProcessorParameter("RD53Aslope",
                               "ADC slope for chip RD53A",
                               _fe_slope,
                               (float)0.1);
}


void MuonCVXDRealDigitiser::init()
{ 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;
    _totEntries = 0;
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

    if( STHcol == nullptr ) return;
    std::string encoder_str { STHcol->getParameters().getStringVal(lcio::LCIO::CellIDEncoding) };
    CellIDDecoder<TrackerHitPlaneImpl> cellid_decoder { encoder_str };

    LCCollectionVec *THcol = new LCCollectionVec(LCIO::TRACKERHITPLANE);

    HitTemporalIndexes t_index{STHcol};

    for (int layer = 0; layer < _numberOfLayers; layer++)
    {
#pragma omp parallel for
        for (int ladder = 0; ladder < _laddersInLayer[layer]; ladder++)
        {
            int num_segment_x = 1;
            int nun_segment_y = _sensorsPerLadder[layer];

            float start_time = t_index.GetMinTime(layer, ladder);
            if (start_time == HitTemporalIndexes::MAXTIME)
            {
                if (streamlog::out.write<streamlog::DEBUG6>())
#pragma omp critical
                {
                    streamlog::out() << "Undefined min time for layer " << layer
                        << " ladder " << ladder << std::endl;
                }
                continue;
            }
            start_time -= _window_size / 2;

            HKBaseSensor sensor {
                layer, ladder,
                num_segment_x, nun_segment_y,
                _layerLadderLength[layer],
                _layerLadderWidth[layer],
                _layerThickness[layer],
                _pixelSizeX, _pixelSizeY,
                encoder_str, _barrelID,
                _threshold, _fe_slope,
                start_time, _window_size
            };
            
            if (sensor.GetStatus() != MatrixStatus::ok and streamlog::out.write<streamlog::ERROR>())
            {
                if (sensor.GetStatus() == MatrixStatus::pixel_number_error)
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
                t_index, sensor,
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

            while(t_window.active())
            {
                t_window.process();

                SegmentDigiHitList hit_buffer {};
                sensor.buildHits(hit_buffer);
                if (hit_buffer.size() == 0) continue;

                vector<TrackerHitPlaneImpl*> reco_buffer;
                reco_buffer.assign(hit_buffer.size(), nullptr);

                int idx = 0;
                for (SegmentDigiHit digiHit : hit_buffer)
                {
                    TrackerHitPlaneImpl *recoHit = new TrackerHitPlaneImpl();
                    recoHit->setEDep((digiHit.charge / _electronsPerKeV) * dd4hep::keV);

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
                    float s_offset = sensor.GetSensorCols() * sensor.GetPixelSizeY();
                    s_offset *= (float(segment_id) + 0.5);
                    s_offset -= sensor.GetHalfLength();

                    Vector2D oldPos(loc_pos[0] * dd4hep::mm, loc_pos[1] * dd4hep::mm - s_offset);
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

                    reco_buffer[idx] = recoHit;
                    idx++;
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
                                             << sensor.GetLayer() << ":" << sensor.GetLadder() 
                                             << ":" << cellid_decoder(recoHit)["sensor"] << std::endl
                                             << "- global position (x,y,z,t) = " << recoHit->getPosition()[0] 
                                             << ", " << recoHit->getPosition()[1] 
                                             << ", " << recoHit->getPosition()[2] 
                                             << ", " << recoHit->getTime() << std::endl
                                             << "- charge = " << recoHit->getEDep() << std::endl;
                        }
                    }
                }
            }
        }
    }

    evt->addCollection(THcol, _outputCollectionName.c_str());
}

void MuonCVXDRealDigitiser::check(LCEvent *evt)
{}

void MuonCVXDRealDigitiser::end()
{
    streamlog_out(DEBUG) << "   end called  " << std::endl;
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



