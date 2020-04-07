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

#include "MyG4UniversalFluctuationForSi.h"


using namespace lcio ;
using namespace marlin ;


/**  Digitizer for Simulated Hits in the Vertex Detector. <br>
 * Digitization follows the procedure adopted in the CMS software package.
 * @param CollectionName Name of the MCParticle collection
 */
class MuonCVXDDigitiser : public Processor {
  
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
    double _tanLorentzAngle;
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
    std::vector<float> _bkgdHitsInLayer;

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
    std::vector<float> _layerLadderGap{};
    std::vector<float> _layerLadderWidth{};
 
};

#endif



