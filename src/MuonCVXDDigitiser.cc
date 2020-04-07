#include "MuonCVXDDigitiser.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"

#include "gsl/gsl_math.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"


using namespace lcio ;
using namespace marlin ;

using dd4hep::Detector;
using dd4hep::DetElement;
using dd4hep::rec::ZPlanarData;

MuonCVXDDigitiser aMuonCVXDDigitiser ;

MuonCVXDDigitiser::MuonCVXDDigitiser() : Processor("MuonCVXDDigitiser")
{
    _description = "MuonCVXDDigitiser should create VTX TrackerHits from SimTrackerHits";

    registerInputCollection(LCIO::SIMTRACKERHIT,
                           "CollectionName", 
                           "Name of the SimTrackerHit collection",
                           _colName,
                           std::string("VXDCollection"));

    registerOutputCollection(LCIO::TRACKERHIT,
                            "OutputCollectionName", 
                            "Name of the output TrackerHit collection",
                            _outputCollectionName,
                            std::string("VTXTrackerHits"));

    registerOutputCollection(LCIO::LCRELATION,
                            "RelationColName", 
                            "Name of the output VTX trackerhit relation collection",
                            _colVTXRelation,
                            std::string("VTXRelation"));

    registerProcessorParameter("SubDetectorName", 
                               "Name of Vertex detector",
                               _subDetName,
                               std::string("VertexBarrel"));

    registerProcessorParameter("TanLorentz",
                               "Tangent of Lorentz Angle",
                               _tanLorentzAngle,
                               (double)0.8);

    registerProcessorParameter("CutOnDeltaRays",
                               "Cut on delta-ray energy (MeV)",
                               _cutOnDeltaRays,
                               (double)0.030);

    registerProcessorParameter("Diffusion",
                               "Diffusion coefficient (in mm) for layer thickness",
                               _diffusionCoefficient,
                               (double)0.002);

    registerProcessorParameter("PixelSizeX",
                               "Pixel Size X",
                               _pixelSizeX,
                               (double)0.025);

    registerProcessorParameter("PixelSizeY",
                               "Pixel Size Y",
                               _pixelSizeY,
                               (double)0.025);

    registerProcessorParameter("Debug",
                               "Debug option",
                               _debug,
                               int(0));

    registerProcessorParameter("ElectronsPerKeV",
                               "Electrons per keV",
                               _electronsPerKeV,
                               (double)270.3);

    std::vector<float> bkgdHitsInLayer;
    bkgdHitsInLayer.push_back(34400.);
    bkgdHitsInLayer.push_back(23900.);
    bkgdHitsInLayer.push_back(9600.);
    bkgdHitsInLayer.push_back(5500.);
    bkgdHitsInLayer.push_back(3100.);    
    registerProcessorParameter("BackgroundHitsPerLayer",
                               "Background Hits per Layer",
                               _bkgdHitsInLayer,
                               bkgdHitsInLayer);

    registerProcessorParameter("SegmentLength",
                               "Segment Length",
                               _segmentLength,
                               double(0.005));

    registerProcessorParameter("WidthOfCluster",
                               "Width of cluster",
                               _widthOfCluster,
                               double(3.0));

    registerProcessorParameter("Threshold",
                               "Cell Threshold in electrons",
                               _threshold,
                               200.);

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

    registerProcessorParameter("UseMCPMomentum",
                               "Use Particle Momentum",
                               _useMCPMomentum,
                               int(1));

    registerProcessorParameter("EnergyLoss",
                               "Energy Loss keV/mm",
                               _energyLoss,
                               double(280.0));

    registerProcessorParameter("RemoveDRayPixels",
                               "Remove D-Ray Pixels",
                               _removeDrays,
                               int(1));

    registerProcessorParameter("GenerateBackground",
                               "Generate Background",
                               _generateBackground,
                               int(0));
}



void MuonCVXDDigitiser::init()
{ 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;
    _totEntries = 0;
    _fluctuate = new MyG4UniversalFluctuationForSi();
}


void MuonCVXDDigitiser::processRunHeader( LCRunHeader* run)
{ 
    _nRun++ ;

    Detector& theDetector = Detector::getInstance();
    DetElement vxBarrel = theDetector.detector(_subDetName);
    ZPlanarData&  zPlanarData = *vxBarrel.extension<ZPlanarData>();
    _numberOfLayers  = zPlanarData.layers.size();

    _laddersInLayer.resize(_numberOfLayers);
    _layerHalfPhi.resize(_numberOfLayers);
    _layerHalfThickness.resize(_numberOfLayers);
    _layerThickness.resize(_numberOfLayers);
    _layerRadius.resize(_numberOfLayers);
    _layerLadderLength.resize(_numberOfLayers);
    _layerLadderWidth.resize(_numberOfLayers);
    _layerLadderHalfWidth.resize(_numberOfLayers);
    _layerActiveSiOffset.resize(_numberOfLayers);
} 



void MuonCVXDDigitiser::processEvent( LCEvent * evt )
{ 
    LCCollection * STHcol = NULL;
    try
    {
        STHcol = evt->getCollection(_colName);
    }
    catch( lcio::DataNotAvailableException ex )
    {
        streamlog_out(WARNING) << _colName << " collection not available" << std::endl;
        STHcol = NULL;
    }

    if( STHcol != NULL )
    {
        // TODO missing implementation
    }

    streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
        << "   in run:  " << evt->getRunNumber() << std::endl ;

    _nEvt ++ ;
}

void MuonCVXDDigitiser::check( LCEvent * evt )
{}


void MuonCVXDDigitiser::end()
{}

