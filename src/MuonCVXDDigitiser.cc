#include "MuonCVXDDigitiser.h"
#include <iostream>
#include <algorithm>

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

MuonCVXDDigitiser::MuonCVXDDigitiser() :
    Processor("MuonCVXDDigitiser"),
    _nRun(0),
    _nEvt(0),
    _totEntries(0),
    _currentLayer(0)
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
    DetElement vxBarrel = theDetector.detector(_subDetName);              // TODO check missing barrel
    ZPlanarData&  zPlanarData = *vxBarrel.extension<ZPlanarData>();       // TODO check missing extension
    std::vector<ZPlanarData::LayerLayout> vx_layers = zPlanarData.layers;
    _numberOfLayers  = vx_layers.size();

    _laddersInLayer.resize(_numberOfLayers);
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
        _laddersInLayer[curr_layer] = z_layout.ladderNumber;

        _layerHalfPhi[curr_layer] = M_PI/((double)_laddersInLayer[curr_layer]);   // TODO investigate

        _layerThickness[curr_layer] = z_layout.thicknessSensitive;

        _layerHalfThickness[curr_layer] = 0.5 * _layerThickness[curr_layer];

        _layerRadius[curr_layer] = z_layout.distanceSensitive + 0.5 * _layerThickness[curr_layer];

        _layerLadderLength[curr_layer] = 2 * z_layout.zHalfSensitive;

        _layerLadderWidth[curr_layer] = z_layout.widthSensitive;

        _layerLadderHalfWidth[curr_layer] = _layerLadderWidth[curr_layer]/2.;

        _layerActiveSiOffset[curr_layer] = - z_layout.offsetSensitive;

        //_layerLadderGap[curr_layer] = laddergaps[curr_layer];

        _layerPhiOffset[curr_layer] = z_layout.phi0;

        curr_layer++;
    }
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
        LCCollectionVec *THcol = new LCCollectionVec(LCIO::TRACKERHIT);
        LCCollectionVec *STHLocCol = NULL;
        if (_produceFullPattern != 0)
        {
            STHLocCol = new LCCollectionVec(LCIO::SIMTRACKERHIT);
        }

        for (int i=0; i < STHcol->getNumberOfElements(); ++i)
        {
            SimTrackerHit * simTrkHit = 
                dynamic_cast<SimTrackerHit*>(STHcol->getElementAt(i));

            ProduceIonisationPoints( simTrkHit );      

            if (_currentLayer < 0 || _currentLayer >= int(_layerRadius.size())) continue;

            ProduceSignalPoints( );

            SimTrackerHitImplVec simTrkHitVec;

            if (_PoissonSmearing != 0) PoissonSmearer( simTrkHitVec );

            if (_electronicEffects != 0) GainSmearer( simTrkHitVec );

            TrackerHitImpl *recoHit = ReconstructTrackerHit( simTrkHitVec );
            if (recoHit == nullptr)
            {
                streamlog_out(DEBUG) << "Number of pixels above threshold = 0 " << std::endl;
                continue;
            }

            TrackerHitToLab( recoHit );
            if (_debug != 0) PrintInfo(simTrkHit, recoHit);

            recoHit->rawHits().push_back(simTrkHit);
            if (_produceFullPattern == 0)
            {
                recoHit->rawHits().push_back(simTrkHit);
            }
            else
            {
                int nSimHits = int( simTrkHitVec.size() );
                for (int iS=0;iS<nSimHits;++iS)
                {
                    SimTrackerHitImpl * sth = simTrkHitVec[iS];
                    float charge = sth->getEDep();
                    if ( charge >_threshold)
                    {
                        SimTrackerHitImpl * newsth = new SimTrackerHitImpl();
                        double spos[3];
                        double sLab[3];
                        for (int iC=0;iC<3;++iC) 
                        spos[iC] = sth->getPosition()[iC];
                        TransformToLab(spos,sLab);
                        newsth->setPosition(sLab);
                        newsth->setEDep(charge);
                        STHLocCol->addElement(newsth);
                        recoHit->rawHits().push_back( newsth );
                    }
                }
            }

            float pointResoRPhi=0.004;
            float pointResoZ=0.004;
            float covMat[TRKHITNCOVMATRIX] = {
                0.,0.,
                pointResoRPhi * pointResoRPhi,
                0.,0.,
                pointResoZ*pointResoZ
            };
            recoHit->setCovMatrix(covMat);      

            recoHit->setType(100+simTrkHit->getCellID0());
            THcol->addElement( recoHit );
            for (int k=0; k < int(simTrkHitVec.size()); ++k)
            {
                SimTrackerHit * hit = simTrkHitVec[k];
                delete hit;
            }     
        }
    }

    streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
        << "   in run:  " << evt->getRunNumber() << std::endl ;

    _nEvt ++ ;
}

void MuonCVXDDigitiser::check(LCEvent *evt)
{}

void MuonCVXDDigitiser::end()
{
    streamlog_out(DEBUG) << "   end called  " << std::endl;
    delete _fluctuate;
}

/** Function calculates local coordinates of the sim hit 
 * in the given ladder and local momentum of particle. 
 * Also returns module number and ladder 
 * number.
 * Local coordinate system within the ladder 
 * is defined as following :  <br> 
 *    - x axis lies in the ladder plane and orthogonal to the beam axis <br>
 *    - y axis is perpendicular to the ladder plane <br>
 *    - z axis lies in the ladder plane and parallel to the beam axis <br>
 * 
 *    Encoding of modules: <br>
 *    ======================  <br>
 *    - 0 = left endcap <br>
 *    - 1 = left ladder in the barrel <br>
 *    - 2 = right ladder in the barrel <br>
 *    - 3 = right endcap <br>
 * 
 */
void MuonCVXDDigitiser::FindLocalPosition(SimTrackerHit *hit, 
                                          double *localPosition,
                                          double *localDirection)
{
    double xLab[3] = {
        hit->getPosition()[0],
        hit->getPosition()[1],
        hit->getPosition()[2]
    };

    int layer = -1;

    double RXY = sqrt(xLab[0] * xLab[0] + xLab[1] * xLab[1]);

    layer = hit->getCellID0() - 1;                          
    _currentLayer = layer;

    if (layer < 0 || layer > _numberOfLayers) 
    return;

    int module = (xLab[2] < 0.0 ) ? 1 : 2;
    _currentModule = module;

    double Momentum[3];
    if (hit->getMCParticle())
    {
        for (int j=0; j<3; ++j)
            Momentum[j] = hit->getMCParticle()->getMomentum()[j];
    }
    else
    {
        for (int j=0; j<3; ++j)
            Momentum[j] = hit->getMomentum()[j];
    }

    _currentParticleMass = 0;
    if (hit->getMCParticle())
        _currentParticleMass   = hit->getMCParticle()->getMass();
    if (_currentParticleMass < 0.510e-3)
        _currentParticleMass = 0.510e-3;  

    _currentParticleMomentum = 0.0;
    for (int i=0; i<3; ++i)
        _currentParticleMomentum += Momentum[i] * Momentum[i];
    _currentParticleMomentum = sqrt(_currentParticleMomentum);

    double PXY = sqrt(Momentum[0] * Momentum[0] + Momentum[1] * Momentum[1]);

    double PhiInLab = (double)atan2(xLab[1],xLab[0]);
    if (PhiInLab < 0.0) PhiInLab += M_2_PI;
    double PhiInLabMom = atan2(Momentum[1],Momentum[0]);
    if (PhiInLabMom < 0.0) PhiInLabMom += M_2_PI;
    double Radius = _layerRadius[layer];

    double Phi0 = _layerPhiOffset[layer];

    int nLadders = _laddersInLayer[layer];

    double dPhi = 2.0*_layerHalfPhi[layer];

    double PhiLadder=0;
    double PhiInLocal=0;

    if (nLadders > 2) // laddered structure
    {
        int iLadder=0;
        for (int ic = 0; ic < nLadders; ++ic)
        {
            PhiLadder = double(ic)*dPhi + Phi0;
            PhiInLocal = PhiInLab - PhiLadder;
            if (RXY*cos(PhiInLocal)-Radius > -_layerThickness[layer] && 
                RXY*cos(PhiInLocal)-Radius < _layerThickness[layer])
            {
                iLadder = ic;
                break;
            }
        }

        double PhiLocalMom = PhiInLabMom - PhiLadder;
        localPosition[0] = RXY * sin(PhiInLocal);
        localPosition[1] = xLab[2];
        localPosition[2] = RXY * cos(PhiInLocal) - Radius;
        localDirection[0]=PXY * sin(PhiLocalMom);
        localDirection[1]=Momentum[2];
        localDirection[2]=PXY * cos(PhiLocalMom);
        _currentPhi = PhiLadder;
    }  
    else // cyllindrical structure
    {
        localPosition[0]=0.0;
        localPosition[1]=xLab[2];
        localPosition[2]=RXY - Radius;
        double PhiLocalMom = PhiInLabMom - PhiInLab;
        localDirection[0]=PXY * sin(PhiLocalMom);
        localDirection[1]=Momentum[2];
        localDirection[2]=PXY * cos(PhiLocalMom); 
        _currentPhi = PhiInLab;
    }
}

void MuonCVXDDigitiser::ProduceIonisationPoints(SimTrackerHit *hit)
{
    double pos[3];
    double dir[3];
    double entry[3];
    double exit[3];

    FindLocalPosition( hit, pos, dir);

    if (_currentLayer < 0 || _currentLayer > _numberOfLayers) 
    return;

    entry[2] = -_layerHalfThickness[_currentLayer]; 
    exit[2] = _layerHalfThickness[_currentLayer];

    for (int i=0; i<2; ++i) {
        entry[i]=pos[i]+dir[i]*(entry[2]-pos[2])/dir[2];
        exit[i]=pos[i]+dir[i]*(exit[2]-pos[2])/dir[2];
    }

    for (int i=0; i<3; ++i) {
        _currentLocalPosition[i] = pos[i];
        _currentEntryPoint[i] = entry[i];
        _currentExitPoint[i] = exit[i];
    }

    double tanx = dir[0]/dir[2];
    double tany = dir[1]/dir[2];  
    double trackLength = std::min(1.0e+3,
        _layerThickness[_currentLayer] * sqrt(1.0 + tanx * tanx + tany * tany));
    double dEmean = 1e-6*_energyLoss * trackLength;  

    _numberOfSegments = int(trackLength / _segmentLength) + 1;
    dEmean = dEmean / ((double)_numberOfSegments);
    _ionisationPoints.resize(_numberOfSegments);

    _eSum = 0.0;

    double segmentLength = trackLength/((double)_numberOfSegments);
    _segmentDepth = _layerThickness[_currentLayer]/((double)_numberOfSegments);

    for (int i=0; i<_numberOfSegments; ++i)
    {
        double z = -_layerHalfThickness[_currentLayer] + ((double)(i) + 0.5) * _segmentDepth;
        double x = pos[0] + dir[0] * (z - pos[2]) / dir[2];
        double y = pos[1] + dir[1] * (z - pos[2]) / dir[2];
        IonisationPoint ipoint;
        double de = _fluctuate->SampleFluctuations(double(1000. * _currentParticleMomentum),
                                                   double(1000. * _currentParticleMass),
                                                   _cutOnDeltaRays,
                                                   segmentLength,
                                                   double(1000.*dEmean)) / 1000.;
        _eSum = _eSum + de;
        ipoint.eloss = de;
        ipoint.x = x;
        ipoint.y = y;
        ipoint.z = z;
        _ionisationPoints[i] = ipoint;
    }
}

void MuonCVXDDigitiser::ProduceSignalPoints()
{}

void MuonCVXDDigitiser::PoissonSmearer(SimTrackerHitImplVec &simTrkVec)
{}

void MuonCVXDDigitiser::GainSmearer(SimTrackerHitImplVec &simTrkVec)
{}

TrackerHitImpl *MuonCVXDDigitiser::ReconstructTrackerHit(SimTrackerHitImplVec &simTrkVec)
{
    return nullptr;
}

void MuonCVXDDigitiser::TrackerHitToLab(TrackerHitImpl *recoHit)
{}

void MuonCVXDDigitiser::TransformToLab(double *xLoc, double *xLab)
{}

void MuonCVXDDigitiser::PrintInfo(SimTrackerHit *simTrkHit, TrackerHitImpl *recoHit)
{}


