#include "MuonCVXDDigitiser.h"
#include <iostream>
#include <algorithm>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"

#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_math.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h" 

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using CLHEP::RandGauss;
using CLHEP::RandPoisson;

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


void MuonCVXDDigitiser::processRunHeader(LCRunHeader* run)
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

        _layerHalfPhi[curr_layer] = M_PI / ((double)_laddersInLayer[curr_layer]);   // TODO investigate

        _layerThickness[curr_layer] = z_layout.thicknessSensitive;

        _layerHalfThickness[curr_layer] = 0.5 * _layerThickness[curr_layer];

        _layerRadius[curr_layer] = z_layout.distanceSensitive + _layerHalfThickness[curr_layer];

        _layerLadderLength[curr_layer] = z_layout.lengthSensor;

        _layerLadderWidth[curr_layer] = z_layout.widthSensitive;

        _layerLadderHalfWidth[curr_layer] = _layerLadderWidth[curr_layer] / 2.;

        _layerActiveSiOffset[curr_layer] = - z_layout.offsetSensitive;

        //_layerLadderGap[curr_layer] = laddergaps[curr_layer];

        _layerPhiOffset[curr_layer] = z_layout.phi0;

        curr_layer++;
    }

    PrintGeometryInfo();
} 



void MuonCVXDDigitiser::processEvent(LCEvent * evt)
{ 
    LCCollection * STHcol = nullptr;
    try
    {
        STHcol = evt->getCollection(_colName);
    }
    catch( lcio::DataNotAvailableException ex )
    {
        streamlog_out(WARNING) << _colName << " collection not available" << std::endl;
        STHcol = nullptr;
    }

    if( STHcol != nullptr )
    {
        LCCollectionVec *THcol = new LCCollectionVec(LCIO::TRACKERHIT);
        LCCollectionVec *STHLocCol = nullptr;
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

            ProduceSignalPoints();

            SimTrackerHitImplVec simTrkHitVec;
            ProduceHits(simTrkHitVec);

            if (_PoissonSmearing != 0) PoissonSmearer(simTrkHitVec);

            if (_electronicEffects != 0) GainSmearer(simTrkHitVec);

            TrackerHitImpl *recoHit = ReconstructTrackerHit(simTrkHitVec);
            if (recoHit == nullptr)
            {
                streamlog_out(DEBUG) << "Number of pixels above threshold = 0 " << std::endl;
                continue;
            }

            TrackerHitToLab(recoHit);
            if (_debug != 0) PrintInfo(simTrkHit, recoHit);

            recoHit->rawHits().push_back(simTrkHit);
            if (_produceFullPattern == 0)
            {
                recoHit->rawHits().push_back(simTrkHit);
            }
            else
            {
                for (int iS = 0; iS < simTrkHitVec.size(); ++iS)
                {
                    SimTrackerHitImpl *sth = simTrkHitVec[iS];
                    float charge = sth->getEDep();
                    if (charge >_threshold)
                    {
                        SimTrackerHitImpl *newsth = new SimTrackerHitImpl();
                        double spos[3];
                        double sLab[3];
                        for (int iC = 0; iC < 3; ++iC) 
                            spos[iC] = sth->getPosition()[iC];
                        TransformToLab(spos,sLab);
                        newsth->setPosition(sLab);
                        newsth->setEDep(charge);
                        STHLocCol->addElement(newsth);
                        recoHit->rawHits().push_back(newsth);
                    }
                }
            }

            // TODO what is this?
            float pointResoRPhi = 0.004;
            float pointResoZ = 0.004;
            float covMat[TRKHITNCOVMATRIX] = {
                0., 0., pointResoRPhi * pointResoRPhi,
                0., 0., pointResoZ * pointResoZ
            };
            recoHit->setCovMatrix(covMat);

            recoHit->setType(100 + simTrkHit->getCellID0());
            THcol->addElement(recoHit);

            for (int k=0; k < int(simTrkHitVec.size()); ++k)
            {
                SimTrackerHit *hit = simTrkHitVec[k];
                delete hit;
            }

            evt->addCollection(THcol, _outputCollectionName.c_str());
            if (_produceFullPattern != 0)
            {
                evt->addCollection(STHLocCol, "VTXPixels");
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
 *    - 1 = left ladder in the barrel <br>
 *    - 2 = right ladder in the barrel <br>
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

    double RXY = sqrt(pow(xLab[0], 2) + pow(xLab[1], 2));

    _currentLayer = -1;
    for (int i = 0; i < _numberOfLayers; ++i)
    {
        if (_laddersInLayer[i] > 2) // laddered structure
        {
            double xmin = _layerRadius[i] - _layerHalfThickness[i];
            double xmax = _layerRadius[i] + _layerHalfThickness[i];
            xmax /= (double)cos(_layerHalfPhi[i]);
            if (RXY > xmin && RXY < xmax)
            {
                _currentLayer = i;
                break;
            }
        }
        else // cyllindrical structure
        {}
    }
    if (_currentLayer < 0) return;

    _currentModule = (xLab[2] < 0.0 ) ? 1 : 2;

    double Momentum[3];
    if (hit->getMCParticle())
    {
        for (int j = 0; j < 3; ++j)
            Momentum[j] = hit->getMCParticle()->getMomentum()[j];
    }
    else
    {
        for (int j = 0; j < 3; ++j)
            Momentum[j] = hit->getMomentum()[j];
    }

    _currentParticleMass = 0.510e-3;
    if (hit->getMCParticle())
        _currentParticleMass = std::max(hit->getMCParticle()->getMass(), _currentParticleMass);

    _currentParticleMomentum = sqrt(pow(Momentum[0], 2) + pow(Momentum[1], 2) 
                                    + pow(Momentum[2], 2));

    double PXY = sqrt(pow(Momentum[0], 2) + pow(Momentum[1], 2));

    double PhiInLab = atan2(xLab[1], xLab[0]);
    if (PhiInLab < 0.0) PhiInLab += M_2_PI;

    double PhiInLabMom = atan2(Momentum[1], Momentum[0]);
    if (PhiInLabMom < 0.0) PhiInLabMom += M_2_PI;

    double Radius = _layerRadius[_currentLayer];

    double Phi0 = _layerPhiOffset[_currentLayer];

    int nLadders = _laddersInLayer[_currentLayer];


    if (nLadders > 2) // laddered structure
    {
        bool ladder_missing = true;

        for (int ic = 0; ic < nLadders; ++ic)
        {
            double PhiLadder = double(ic) * 2.0 * _layerHalfPhi[_currentLayer] + Phi0;
            double PhiInLocal = PhiInLab - PhiLadder;
            double delta_radius = RXY * cos(PhiInLocal) - Radius;

            if (abs(delta_radius) <= _layerHalfThickness[_currentLayer])
            {
                double PhiLocalMom = PhiInLabMom - PhiLadder;
                localPosition[0] = RXY * sin(PhiInLocal);
                localPosition[1] = xLab[2];
                localPosition[2] = delta_radius;
                localDirection[0]= PXY * sin(PhiLocalMom);
                localDirection[1]= Momentum[2];
                localDirection[2]= PXY * cos(PhiLocalMom);
                _currentPhi = PhiLadder;
                ladder_missing = false;
                break;
            }
        }
        if (ladder_missing)
            streamlog_out(DEBUG) << "Hit out of ladder" << std::endl;
    }  
    else // cyllindrical structure TODO is it necessary?
    {
        localPosition[0]= 0.0;
        localPosition[1]= xLab[2];
        localPosition[2]= RXY - Radius;
        double PhiLocalMom = PhiInLabMom - PhiInLab;
        localDirection[0]= PXY * sin(PhiLocalMom);
        localDirection[1]= Momentum[2];
        localDirection[2]= PXY * cos(PhiLocalMom);
        _currentPhi = PhiInLab;
    }
}

void MuonCVXDDigitiser::ProduceIonisationPoints(SimTrackerHit *hit)
{
    double pos[3];
    double dir[3];
    double entry[3];
    double exit[3];

    FindLocalPosition(hit, pos, dir);

    if (_currentLayer < 0 || _currentLayer > _numberOfLayers) 
    return;

    entry[2] = -_layerHalfThickness[_currentLayer]; 
    exit[2] = _layerHalfThickness[_currentLayer];

    for (int i = 0; i < 2; ++i) {
        entry[i] = pos[i] + dir[i] * (entry[2] - pos[2]) / dir[2];
        exit[i]= pos[i] + dir[i] * (exit[2] - pos[2]) / dir[2];
    }

    for (int i = 0; i < 3; ++i) {
        _currentLocalPosition[i] = pos[i];
        _currentEntryPoint[i] = entry[i];
        _currentExitPoint[i] = exit[i];
    }

    double tanx = dir[0] / dir[2];
    double tany = dir[1] / dir[2];  
    double trackLength = std::min(1.0e+3,
        _layerThickness[_currentLayer] * sqrt(1.0 + pow(tanx, 2) + pow(tany, 2)));
    _numberOfSegments = ceil(trackLength / _segmentLength);

    double dEmean = (1e-6 * _energyLoss * trackLength) / ((double)_numberOfSegments);

    _ionisationPoints.resize(_numberOfSegments);

    _eSum = 0.0;

    // TODO _segmentLength may be different from segmentLength, is it ok?
    double segmentLength = trackLength / ((double)_numberOfSegments);
    _segmentDepth = _layerThickness[_currentLayer] / ((double)_numberOfSegments);

    for (int i = 0; i < _numberOfSegments; ++i)
    {
        double z = -_layerHalfThickness[_currentLayer] + ((double)(i) + 0.5) * _segmentDepth;
        double x = pos[0] + dir[0] * (z - pos[2]) / dir[2];
        double y = pos[1] + dir[1] * (z - pos[2]) / dir[2];
        double de = _fluctuate->SampleFluctuations(double(1000. * _currentParticleMomentum),
                                                   double(1000. * _currentParticleMass),
                                                   _cutOnDeltaRays,
                                                   segmentLength,
                                                   double(1000. * dEmean)) / 1000.;
        _eSum = _eSum + de;

        IonisationPoint ipoint;
        ipoint.eloss = de;
        ipoint.x = x;
        ipoint.y = y;
        ipoint.z = z;
        _ionisationPoints[i] = ipoint;
    }
}

void MuonCVXDDigitiser::ProduceSignalPoints()
{
    _signalPoints.resize(_numberOfSegments);

    // run over ionisation points
    for (int i = 0; i < _numberOfSegments; ++i)
    {
        IonisationPoint ipoint = _ionisationPoints[i];
        double z = ipoint.z;
        double x = ipoint.x;
        double y = ipoint.y;
        double DistanceToPlane = _layerHalfThickness[_currentLayer] - z;
        double xOnPlane = x + _tanLorentzAngleX * DistanceToPlane;
        double yOnPlane = y + _tanLorentzAngleY * DistanceToPlane;
        double DriftLength = DistanceToPlane * sqrt(1.0 + pow(_tanLorentzAngleX, 2) 
                                                        + pow(_tanLorentzAngleY, 2));

        double SigmaDiff = sqrt(DriftLength / _layerThickness[_currentLayer])
                           * _diffusionCoefficient;

        double SigmaX = SigmaDiff * sqrt(1.0 + pow(_tanLorentzAngleX, 2));
        double SigmaY = SigmaDiff * sqrt(1.0 + pow(_tanLorentzAngleY, 2));

        double charge = 1.0e+6 * ipoint.eloss * _electronsPerKeV;

        SignalPoint  spoint;
        spoint.x = xOnPlane;
        spoint.y = yOnPlane;
        spoint.sigmaX = SigmaX;
        spoint.sigmaY = SigmaY;
        spoint.charge = charge;
        _signalPoints[i] = spoint;
    }
}

void MuonCVXDDigitiser::ProduceHits(SimTrackerHitImplVec &simTrkVec)
{  
    simTrkVec.clear();

    for (int i=0; i<_numberOfSegments; ++i)
    {
        SignalPoint spoint = _signalPoints[i];
        double xCentre = spoint.x;
        double yCentre = spoint.y;
        double sigmaX = spoint.sigmaX;
        double sigmaY = spoint.sigmaY;
        double xLo = spoint.x - _widthOfCluster * spoint.sigmaX;
        double xUp = spoint.x + _widthOfCluster * spoint.sigmaX;
        double yLo = spoint.y - _widthOfCluster * spoint.sigmaY;
        double yUp = spoint.y + _widthOfCluster * spoint.sigmaY;
        
        int ixLo, ixUp, iyLo, iyUp;

        TransformXYToCellID(xLo, yLo, ixLo, iyLo);
        TransformXYToCellID(xUp, yUp, ixUp, iyUp);

        for (int ix = ixLo; ix< ixUp + 1; ++ix)
        {
            if (ix < 0) continue;

            for (int iy = iyLo; iy < iyUp + 1; ++iy)
            {
                if (iy < 0) continue;

                double xCurrent, yCurrent;
                TransformCellIDToXY(ix, iy, xCurrent, yCurrent);
                
                gsl_sf_result result;
                int status = gsl_sf_erf_Q_e((xCurrent - 0.5 * _pixelSizeX - xCentre)/sigmaX, &result);
                double LowerBound = 1 - result.val;

                status = gsl_sf_erf_Q_e((xCurrent + 0.5 * _pixelSizeX - xCentre)/sigmaX, &result);
                double UpperBound = 1 - result.val;
                double integralX = UpperBound - LowerBound;

                status = gsl_sf_erf_Q_e((yCurrent - 0.5 * _pixelSizeY - yCentre)/sigmaY, &result);
                LowerBound = 1 - result.val;

                status = gsl_sf_erf_Q_e((yCurrent + 0.5 * _pixelSizeY - yCentre)/sigmaY, &result);
                UpperBound = 1 - result.val;
                double integralY = UpperBound - LowerBound;

                float totCharge = float(spoint.charge * integralX * integralY);

                int iexist = 0;
                int cellID = PX_PER_ROW * ix + iy;                 // TODO check cellID
                SimTrackerHitImpl *existingHit = nullptr;
                for (int iHits = 0; iHits < int(simTrkVec.size()); ++iHits)
                {
                    existingHit = simTrkVec[iHits];
                    int cellid = existingHit->getCellID0();
                    if (cellid == cellID)
                    {
                        iexist = 1;
                        break;
                    }
                }
                if (iexist == 1)
                {
                    float edep = existingHit->getEDep();
                    edep += totCharge;
                    existingHit->setEDep( edep );
                }
                else 
                {
                    SimTrackerHitImpl * hit = new SimTrackerHitImpl();
                    double pos[3] = {
                        xCurrent,
                        yCurrent,
                        _layerHalfThickness[_currentLayer]
                    };
                    hit->setPosition(pos);
                    hit->setCellID0(cellID);
                    hit->setEDep(totCharge);
                    simTrkVec.push_back(hit);
                }
            }
        }
    }
}

/**
 * Function that fluctuates charge (in units of electrons)
 * deposited on the fired pixels according to the Poisson
 * distribution...
 */
void MuonCVXDDigitiser::PoissonSmearer(SimTrackerHitImplVec &simTrkVec)
{
    for (int ihit = 0; ihit < int(simTrkVec.size()); ++ihit)
    {
        SimTrackerHitImpl *hit = simTrkVec[ihit];
        float charge = hit->getEDep();
        float rng;
        if (charge > 1000.) // assume Gaussian
        {
            rng = float(RandGauss::shoot(charge, sqrt(charge)));
        }
        else // assume Poisson
        {
            rng = float(RandPoisson::shoot(charge));
        }
        hit->setEDep(rng);
    }  
}

/**
 * Simulation of electronic noise.
 */
void MuonCVXDDigitiser::GainSmearer(SimTrackerHitImplVec &simTrkVec)
{
    for (int i = 0; i < simTrkVec.size(); ++i)
    {
        double Noise = RandGauss::shoot(0., _electronicNoise);
        SimTrackerHitImpl *hit = simTrkVec[i];
        hit->setEDep(hit->getEDep() + float(Noise));
    }
}

/**
 * Emulates reconstruction of Tracker Hit 
 * Tracker hit position is reconstructed as center-of-gravity 
 * of cluster of fired cells. Position is corrected for Lorentz shift.
 * TODO check the track reco algorithm (is it the one we need?)
 */
TrackerHitImpl *MuonCVXDDigitiser::ReconstructTrackerHit(SimTrackerHitImplVec &simTrkVec)
{
    double pos[3] = {0, 0, 0};
    double charge = 0;

    /* Simple center-of-gravity
    for (int iHit=0; iHit < int(simTrkVec.size()); ++iHit)
    {
        SimTrackerHit *hit = simTrkVec[iHit];
        if (hit->getEDep() <= _threshold) continue;

        charge += hit->getEDep();
        pos[0] += hit->getEDep() * hit->getPosition()[0];
        pos[1] += hit->getEDep() * hit->getPosition()[1];
    }

    if (charge > 0.)
    {
        TrackerHitImpl *recoHit = new TrackerHitImpl();
        recoHit->setEDep(charge);

        pos[0] /= charge;
        pos[0] -= _layerHalfThickness[_currentLayer] * _tanLorentzAngleX;

        pos[1] /= charge;
        pos[1] -= _layerHalfThickness[_currentLayer] * _tanLorentzAngleY;

        recoHit->setPosition(pos);
        return recoHit;
    }
    */

    /* Partial histogram */
    int nPixels = 0;
    int ixmin =  1000000;
    int ixmax = -1000000;
    int iymin =  1000000;
    int iymax = -1000000;

    for (int iHit=0; iHit < int(simTrkVec.size()); ++iHit)
    {
        SimTrackerHit *hit = simTrkVec[iHit];
        if (hit->getEDep() <= _threshold) continue;

        nPixels++;
        charge += hit->getEDep();
        int cellID = hit->getCellID0();
        int ix = cellID / PX_PER_ROW ;
        int iy = cellID - PX_PER_ROW * ix;      

        if (ix > ixmax) ixmax = ix;
        if (ix < ixmin) ixmin = ix;
        if (iy > iymax) iymax = iy;
        if (iy < iymin) iymin = iy;
    }

    if (charge > 0. && nPixels > 0)
    {
        TrackerHitImpl *recoHit = new TrackerHitImpl();
        recoHit->setEDep(charge);

        double _amplX[20];
        double _amplY[20];

        for (int k = 0; k < 20; ++k)
        {
            _amplY[k] = 0.0;
            _amplX[k] = 0.0;
        }

        for (int iHit = 0; iHit < int(simTrkVec.size()); ++iHit)
        {
            SimTrackerHit *hit = simTrkVec[iHit];
            if (hit->getEDep() <= _threshold) continue;

            int cellID = hit->getCellID0();
            int ix = cellID / PX_PER_ROW ;
            int iy = cellID - PX_PER_ROW * ix;
            if ((iy - iymin) < 20)
                _amplY[iy - iymin] = _amplY[iy - iymin] + hit->getEDep();
            if ((ix - ixmin) < 20) 
                _amplX[ix - ixmin] = _amplX[ix - ixmin] + hit->getEDep();        

        }

        double aXCentre = 0;
        double aYCentre = 0;
        for (int i = ixmin + 1; i < ixmax; ++i)
        {
            aXCentre += _amplX[i - ixmin];
        }
        for (int i = iymin + 1; i < iymax; ++i)
        {
            aYCentre += _amplY[i - iymin];
        }
        aXCentre = aXCentre / std::max(1, ixmax - ixmin - 1);
        aYCentre = aYCentre / std::max(1, iymax - iymin - 1);

        double aTot = 0;
        for (int i = ixmin; i < ixmax + 1; ++i)
        {
            double xx, yy;
            aTot += _amplX[i - ixmin];
            TransformCellIDToXY(i, 2, xx, yy);
            if (i != ixmin && i != ixmax)
            {
                pos[0] += xx * aXCentre;
            }
            else
            {
                pos[0] += xx * _amplX[i - ixmin];
            }
        }
        pos[0] /= aTot;
        pos[0] -= _layerHalfThickness[_currentLayer] * _tanLorentzAngleX;

        aTot = 0;
        for (int i = iymin; i < iymax + 1 ;++i)
        {
            double xx, yy;
            TransformCellIDToXY(i, i, xx, yy);
            aTot += _amplY[i - iymin];
            if (i != iymin && i != iymax)
            {
                pos[1] += yy * aYCentre;
            }
            else
            {
                pos[1] += yy * _amplY[i-iymin];
            }
        }
        pos[1] /= aTot;
        pos[1] -= _layerHalfThickness[_currentLayer] * _tanLorentzAngleY;

        recoHit->setPosition(pos);
        return recoHit;
    }

    return nullptr;
}

void MuonCVXDDigitiser::TrackerHitToLab(TrackerHitImpl *recoHit)
{
    double pos[3];
    for (int i = 0; i < 3; ++i)
    {
        pos[i] = recoHit->getPosition()[i];
    }
    double xLab[3];
    TransformToLab(pos, xLab);

    recoHit->setPosition(xLab);
}

/** Function transforms local coordinates in the ladder
 * into global coordinates
 */
void MuonCVXDDigitiser::TransformToLab(double *xLoc, double *xLab)
{
    int nLadders = _laddersInLayer[_currentLayer];
    double Radius = _layerRadius[_currentLayer];

    if (nLadders > 2 ) // laddered structure
    {
        double baseLine = Radius + xLoc[2];
        double PhiInLab = _currentPhi + atan2(xLoc[0], baseLine);
        double RXY = sqrt(pow(baseLine, 2) + pow(xLoc[0], 2));
        xLab[2] = xLoc[1];
        xLab[0] = RXY * cos(PhiInLab);
        xLab[1] = RXY * sin(PhiInLab);
    }
    else // cyllindrical structure
    {
        double baseLine = Radius + xLoc[2];
        double PhiInLab = _currentPhi + xLoc[0] / baseLine;
        xLab[0] = baseLine * cos(PhiInLab);
        xLab[1] = baseLine * sin(PhiInLab);
        xLab[2] = xLoc[1];    
    } 
}

/**
 * Function calculates position in pixel matrix based on the 
 * local coordinates of point in the ladder.
 */
void MuonCVXDDigitiser::TransformXYToCellID(double x, double y, int & ix, int & iy)
{
    int layer = _currentLayer;
    int nladders = _laddersInLayer[layer];
    //double ladderGap = _layerLadderGap[layer];
    double ladderGap = 0;                             // TODO investigate
    double ladderLength = _layerLadderLength[layer];

    double yInLadder = (y < 0.0) ? y + ladderLength : y - ladderGap;

    iy = (yInLadder < 0.0) ? -1 : int(yInLadder / _pixelSizeY);

    double xInLadder = x;
    if (nladders > 2) // laddered structure
    {
        xInLadder += _layerLadderHalfWidth[layer] + _layerActiveSiOffset[layer];
    }
    else // cyllindrical structure  TODO is it necessary
    {
        xInLadder += (_layerRadius[layer] + _layerHalfThickness[layer]) * _currentPhi;
    }

    ix = (xInLadder < 0.0) ? -1 : int(xInLadder / _pixelSizeX);
}

/**
 Function calculates position in the local frame 
 based on the index of pixel in the ladder.
*/
void MuonCVXDDigitiser::TransformCellIDToXY(int ix, int iy, double & x, double & y)
{
    int layer = _currentLayer;
    int nladders = _laddersInLayer[layer];
    //double ladderGap = _layerLadderGap[layer];
    double ladderGap = 0;                             // TODO investigate
    double ladderLength = _layerLadderLength[layer];

    y = (0.5 + double(iy)) * _pixelSizeY;

    if (_currentModule == 1) 
        y -= ladderLength;
    else
        y += ladderGap;

    x = (0.5 + double(ix)) * _pixelSizeX;

    if (nladders > 2) // laddered structure
        x -= _layerLadderHalfWidth[layer] + _layerActiveSiOffset[layer];
    else // cyllindrical structure TODO is it necessary
        x -= (_layerRadius[layer] + _layerHalfThickness[layer]) * _currentPhi;
}

void MuonCVXDDigitiser::PrintGeometryInfo()
{
    streamlog_out(MESSAGE) << "Number of layers: " << _numberOfLayers << std::endl;
    streamlog_out(MESSAGE) << "Pixel size X: " << _pixelSizeX << std::endl;
    streamlog_out(MESSAGE) << "Pixel size Y: " << _pixelSizeY << std::endl;
    streamlog_out(MESSAGE) << "Electrons per KeV: " << _electronsPerKeV << std::endl;
    streamlog_out(MESSAGE) << "Segment depth: " << _segmentDepth << std::endl;
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

void MuonCVXDDigitiser::PrintInfo(SimTrackerHit *simTrkHit, TrackerHitImpl *recoHit)
{
    // TODO TBD
}



