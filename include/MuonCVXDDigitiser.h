#ifndef MuonCVXDDigitiser_h
#define MuonCVXDDigitiser_h 1

#include <string>
#include <vector>

#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/LCIO.h"
#include "IMPL/TrackerHitImpl.h"
#include "IMPL/SimTrackerHitImpl.h"
#include <IMPL/LCCollectionVec.h>
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"
#include "MyG4UniversalFluctuationForSi.h"

// TODO check the following value
#define PX_PER_ROW 100000

using marlin::Processor;

struct IonisationPoint
{
    double x;
    double y;
    double z;
    double eloss;
};

struct SignalPoint
{
    double x;
    double y;
    double sigmaX;
    double sigmaY;
    double charge;
};

typedef std::vector<SimTrackerHitImpl*> SimTrackerHitImplVec;
typedef std::vector<IonisationPoint> IonisationPointVec;
typedef std::vector<SignalPoint> SignalPointVec;

/**  Digitizer for Simulated Hits in the Vertex Detector. <br>
 * Digitization follows the procedure adopted in the CMS software package.
 * Processor produces an output collection of the Tracker Hits. 
 * Collection has name "VTXTrackerHits".
 * @param CollectionName name of input SimTrackerHit collection <br>
 * (default parameter value : "vxd01_VXD", taken from Mokka)
 * @param TanLorentz tangent of the Lorentz angle <br>
 * (default parameter value : 0.8) <br>
 * @param CutOnDeltaRays cut on the energy of delta-electrons (in MeV) <br>
 * (default parameter value : 0.03) <br>
 * @param Diffusion diffusion coefficient for the nominal active layer thickness (in mm) <br>
 * (default parameter value : 0.002) <br>
 * @param LayerThickness thickness of the active Silicon layer (in mm) <br>
 * (default parameter value : 0.03744) <br>
 * @param PixelSizeX pixel size along direction perpendicular to beam axis (in mm) <br>
 * (default value : 0.025) <br>
 * @param PixelSizeY pixel size along beam axis (in mm) <br>
 * (default value : 0.025) <br>
 * @param ElectronsPerMeV number of electrons produced per MeV of deposited energy <br>
 * (default parameter value : 270.3) <br>
 * @param Threshold threshold on charge deposited on one pixel (in electons) <br>
 * (default parameter value : 200.0) <br>
 * @param LaddersInLayer vector of integers, numbers of phi-ladders in each layer <br>
 * (default parameter values : 8, 8, 12, 16, 20; taken from Mokka database for VXD model vxd01) <br>
 * @param LadderRadius vector of doubles, radii of layers (in mm) <br>
 * (default parameter values : 15.301, 26.301, 38.301, 49.301, 60.301; taken from Mokka database 
 * for VXD model vxd01) <br>
 * @param ActiveLadderOffset vector of doubles, offset of active Si layer along phi-angle (in mm) for each layer <br>
 * (default parameter values : 1.455, 1.39866, 2.57163, 3.59295, 4.42245) <br>
 * @param LadderHalfWidth vector of doubles, half-width of the ladders in each layer (in mm)<br>
 * (default parameter values : 6.5, 11.0, 11.0, 11.0, 11.0; taken from Mokka database for VXD model vxd01) <br>
 * @param LadderGaps vector of doubles, gaps between two symmetric ladders (+/-z) in mm<br>
 * (default parameter values : 0.0, 0.04, 0.04, 0.04, 0.04; taken from Mokka database) <br>
 * @param PhiOffset vector of doubles, offset in phi angle for starting ladder in each layer <br>
 * (default parameter values : 0, 0, 0, 0, 0; taken from Mokka database for VXD model vxd01) <br>
 * @param SegmentLength segment length along track path which is used to subdivide track into segments (in mm).
 * The number of track subsegments is calculated as int(TrackLengthWithinActiveLayer/SegmentLength)+1 <br>
 * (default parameter value : 0.005) <br>
 * @param WidthOfCluster defines width in Gaussian sigmas to perform charge integration for 
 * a given pixel <br>
 * (default parameter value : 3.0) <br>
 * @param ElectronicEffects flag to switch on gaussian smearing of signal (electronic noise) <br>
 * (default parameter value : 1) <br>
 * @param ElectronicNoise electronic noise in electrons <br>
 * (default parameter value : 100) <br>
 * @param GenerateBackground flag to switch on pseudo-generation of pair background hits
 * Background hits are uniformly generated in cosQ and Phi <br>
 * (default parameter value : 0)
 * @param BackgroundHitsPerLayer expected mean value of background hits accumulated in each layer 
 * over readout time. This number is calculated as Number of background hits per bunch crossing times
 * number of bunch crossings over integration time. <br>
 * (default values : 34400 23900 9600 5500 3100 corresponding to 100 overlaid bunch crossings) <br>
 * @param Debug boolean variable, if set to 1, printout is activated <br>
 * (default parameter value : 0) <br>
 * <br>
 */
class MuonCVXDDigitiser : public Processor
{  
public:
  
    virtual Processor*  newProcessor() { return new MuonCVXDDigitiser ; }

    MuonCVXDDigitiser();

    /** Called at the begin of the job before anything is read.
    * Use to initialize the processor, e.g. book histograms.
    */
    virtual void init();

    /** Called for every run.
    */
    virtual void processRunHeader( LCRunHeader* run );

    /** Called for every event - the working horse.
    */
    virtual void processEvent( LCEvent * evt ); 

    virtual void check( LCEvent * evt );

    /** Called after data processing for clean up.
    */
    virtual void end();

protected:

    int _nRun;
    int _nEvt;
    int _debug;
    int _totEntries;
    std::string _subDetName;

    // input/output collections
    std::string _colName;
    std::string _outputCollectionName;
    std::string _colVTXRelation;

    // processor parameters
    double _tanLorentzAngleX;
    double _tanLorentzAngleY;
    double _cutOnDeltaRays;
    double _diffusionCoefficient;
    double _pixelSizeX;
    double _pixelSizeY;
    double _electronsPerKeV;
    double _segmentLength;
    double _widthOfCluster;
    double _threshold;
    double _electronicNoise;
    double _energyLoss;
    int _PoissonSmearing;
    int _electronicEffects;
    int _produceFullPattern;
    int _useMCPMomentum;
    int _removeDrays;
    int _generateBackground;
    std::vector<float> _bkgdHitsInLayer;  // TODO is it necessary

    MyG4UniversalFluctuationForSi *_fluctuate;

    // geometry
    int _numberOfLayers;
    std::vector<int>   _laddersInLayer{};
    std::vector<float> _layerRadius{};
    std::vector<float> _layerThickness{};
    std::vector<float> _layerHalfThickness{};
    std::vector<float> _layerLadderLength{};
    std::vector<float> _layerLadderHalfWidth{};
    std::vector<float> _layerPhiOffset{};
    std::vector<float> _layerActiveSiOffset{};
    std::vector<float> _layerHalfPhi{};
    //std::vector<float> _layerLadderGap{};
    std::vector<float> _layerLadderWidth{};
    const dd4hep::rec::SurfaceMap* _map ;


    // internal state
    int _currentLayer;
    int _currentModule;
    int _currentLadder;
    int _numberOfSegments;
    double _currentParticleMass;
    double _currentParticleMomentum;
    double _currentPhi;
    double _eSum;
    double _segmentDepth;
    double _currentLocalPosition[3];
    double _currentEntryPoint[3];
    double _currentExitPoint[3];
    IonisationPointVec _ionisationPoints;
    SignalPointVec _signalPoints;

    void ProduceIonisationPoints(SimTrackerHit *hit);
    void ProduceSignalPoints();
    void ProduceHits(SimTrackerHitImplVec &simTrkVec);
    void PoissonSmearer(SimTrackerHitImplVec &simTrkVec);
    void GainSmearer(SimTrackerHitImplVec &simTrkVec);
    TrackerHitImpl *ReconstructTrackerHit(SimTrackerHitImplVec &simTrkVec);
    void TrackerHitToLab(TrackerHitImpl *recoHit);
    void TransformToLab(double *xLoc, double *xLab);
    void FindLocalPosition(SimTrackerHit *hit, double *localPosition, double *localDirection);
    void TransformXYToCellID(double x, double y, int & ix, int & iy);
    void TransformCellIDToXY(int ix, int iy, double & x, double & y);

    void PrintGeometryInfo();
    void PrintInfo(SimTrackerHit *simTrkHit, TrackerHitImpl *recoHit);
};

#endif



