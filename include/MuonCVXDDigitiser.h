#ifndef MuonCVXDDigitiser_h
#define MuonCVXDDigitiser_h 1

// Standard
#include <string>
#include <vector>
#include <tuple>

// k4FWCore
#include <k4FWCore/Transformer.h>

// edm4hep
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/MutableSimTrackerHit.h>
#include <edm4hep/MutableTrackerHitPlane.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <edm4hep/TrackerHitSimTrackerHitLink.h>
#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>

// DD4hep
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"
#include "MyG4UniversalFluctuationForSi.h"

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

typedef std::vector<edm4hep::MutableSimTrackerHit*> MutableSimTrackerHitVec;
typedef std::vector<IonisationPoint> IonisationPointVec;
typedef std::vector<SignalPoint> SignalPointVec;

struct InternalState
{
    int    currentLayer;
    int    currentLadder;
    int    numberOfSegments;
    double currentParticleMass;
    double currentParticleMomentum;
    double currentPhi;
    double eSum;
    double segmentDepth;
    edm4hep::Vector3d  currentLocalPosition;
    edm4hep::Vector3d  currentEntryPoint;
    edm4hep::Vector3d  currentExitPoint;
    IonisationPointVec ionisationPoints;
    SignalPointVec     signalPoints;
};


/**  Digitizer for Simulated Hits in the Vertex Detector. <br>
 * Digitization follows the procedure adopted in the CMS software package. 
 * See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePixelDigitization
 * 
 * @param CollectionName name of input SimTrackerHit collection <br>
 * (default parameter value : "VXDCollection")
 * @param OutputCollectionName name of output TrackerHitsPlane collection <br>
 * (default parameter value : "VTXTrackerHits")
 * @param RelationColName name of the TrackerHitSimTrackerHitLink collection <br>
 * (default parameter value : "VTXTrackerHitRelations")
 * @param SubDetectorName name of the detector <br>
 * (default parameter value : "VertexBarrel")
 * @param TanLorentz tangent of the Lorentz angle <br>
 * (default parameter value : 0.8) <br>
 * @param TanLorentzY tangent of the Lorentz angle along Y <br>
 * (default parameter value : 0) <br>
 * @param CutOnDeltaRays cut on the energy of delta-electrons (in MeV) <br>
 * (default parameter value : 0.03) <br>
 * @param Diffusion diffusion coefficient for the nominal active layer thickness (in mm) <br>
 * (default parameter value : 0.002) <br>
 * @param PixelSizeX pixel size along direction perpendicular to beam axis (in mm) <br>
 * (default value : 0.025) <br>
 * @param PixelSizeY pixel size along beam axis (in mm) <br>
 * (default value : 0.025) <br>
 * @param ElectronsPerMeV number of electrons produced per MeV of deposited energy <br>
 * (default parameter value : 270.3) <br>
 * @param Threshold threshold on charge deposited on one pixel (in electons) <br>
 * (default parameter value : 200.0) <br>
 * @param SegmentLength segment length along track path which is used to subdivide track into segments (in mm).
 * The number of track subsegments is calculated as int(TrackLengthWithinActiveLayer/SegmentLength)+1 <br>
 * (default parameter value : 0.005) <br>
 * @param WidthOfCluster defines width in Gaussian sigmas to perform charge integration for 
 * a given pixel <br>
 * (default parameter value : 3.0) <br>
 * @param PoissonSmearing bool to switch on gaussian smearing of electrons collected on pixels <br>
 * (default parameter value : 1) <br>
 * @param ThresholdSmearSigma: sigma of Gaussian used in threshold smearing <br>
 * (default parameter value: 25) <br>
 * @param ChargeDigitize bool to switch on charge discretization <br>
 * (default parameter value: 1) <br>
 * @param ChargeDigitizeNumBits number of bits used to determine bins for charge discretization <br>
 * (default parameter value: 4) <br>
 * @param ChargeDigitizeBinning binning scheme for charge discretization
 * (default parameter value; 1) <br>
 * @param ElectronicEffects bool to switch on gaussian smearing of signal (electronic noise) <br>
 * (default parameter value : 1) <br>
 * @param ElectronicNoise electronic noise in electrons <br>
 * (default parameter value : 100) <br>
 * @param StoreFiredPixels bool to store also the fired pixels (collection names: "VTXPixels") <br>
 * (default parameter value : 0) <br>
 * @param EnergyLoss Energy loss in keV/mm <br>
 * (default parameter value : 280.0) <br>
 * @param MaxEnergyDelta max delta in energy hit (difference from the hit simulated charged and the comuputed ones) in electrons <br>
 * (default parameter value : 100) <br>
  * @param MaxTrackLength Maximum values for track path length inside the ladder (in mm)", <br>
 * (default parameter value : 10) <br> 
 * <br>
 */
class MuonCVXDDigitiser : public k4FWCore::MultiTransformer<std::tuple<edm4hep::SimTrackerHitCollection,
                                                                       edm4hep::TrackerHitPlaneCollection,
                                                                       edm4hep::TrackerHitSimTrackerHitLinkCollection,
                                                                       edm4hep::TrackerHitSimTrackerHitLinkCollection>(
                                                                 const edm4hep::SimTrackerHitCollection&)>
{  
public:  

    MuonCVXDDigitiser(const std::string& name, ISvcLocator* svcLoc);

    /** Called at the begin of the job before anything is read.
    * Use to initialize the processor, e.g. book histograms.
    */
    StatusCode initialize();

    /** Called for every event - the working horse.
    */
    std::tuple<edm4hep::SimTrackerHitCollection,
               edm4hep::TrackerHitPlaneCollection,
               edm4hep::TrackerHitSimTrackerHitLinkCollection,
               edm4hep::TrackerHitSimTrackerHitLinkCollection> operator()(
         const edm4hep::SimTrackerHitCollection& STHcol) const; 

    /** Called after data processing for clean up.
    */
    StatusCode finalize();

protected:

    bool isBarrel;
    bool isVertex;
    bool isInnerTracker;
    bool isOuterTracker;

    // processor 
    Gaudi::Property<std::string> m_subDetName{this, "SubDetectorName", std::string("VertexBarrel"), "Name of Vertex detector"};
    Gaudi::Property<double> m_tanLorentzAngleX{this, "TanLorentz", (double)0.8, "Tangent of Lorentz Angle"};
    Gaudi::Property<double> m_tanLorentzAngleY{this, "TanLorentzY", (double)0, "Tangent of Lorentz Angle along Y"};
    Gaudi::Property<double> m_cutOnDeltaRays{this, "CutOnDeltaRays", (double)0.030, "Cut on delta-ray energy (MeV)"};
    Gaudi::Property<double> m_diffusionCoefficient{this, "DiffusionCoefficient", (double)0.07, "Diffusion coefficient, sqrt(D / mu / V)."};
    Gaudi::Property<double> m_pixelSizeX{this, "PixelSizeX", (double)0.025, "Pixel Size X"};
    Gaudi::Property<double> m_pixelSizeY{this, "PixelSizeY", (double)0.025, "Pixel Size Y"};
    Gaudi::Property<double> m_electronsPerKeV{this, "ElectronsPerKeV", (double)270.3, "Electrons per keV"};
    Gaudi::Property<double> m_segmentLength{this, "SegmentLength", (double)0.005, "Segment Length in mm"};
    Gaudi::Property<double> m_threshold{this, "Threshold", (double)500., "Cell Threshold in electrons"};
    Gaudi::Property<double> m_electronicNoise{this, "ElectronicNoise", (double)80., "electronic noise in electrons"};
    Gaudi::Property<double> m_energyLoss{this, "EnergyLoss", (double)280., "Energy Loss keV/mm"};
    Gaudi::Property<double> m_deltaEne{this, "MaxEnergyDelta", (double)100., "Max delta in energy between G4 prediction and random sampling for each hit in electrons"};
    Gaudi::Property<double> m_maxTrkLen{this, "MaxTrackLength", (double)10., "Maximum values for track length (in mm)"};	
    Gaudi::Property<bool>   m_PoissonSmearing{this, "PoissonSmearing", true, "Apply Poisson smearing of electrons collected on pixels"};
    Gaudi::Property<int>    m_thresholdSmearSigma{this, "ThresholdSmearSigma", 25, "sigma of Gaussian used in threshold smearing, in electrons"};
    Gaudi::Property<bool>   m_DigitizeCharge{this, "DigitizeCharge", true, "Flag to enable Digitization of the charge collected on pixels"};
    Gaudi::Property<int>    m_ChargeDigitizeNumBits{this, "ChargeDigitizeNumBits", 4, "Number of bits used to determine bins for charge discretization"};
    Gaudi::Property<int>    m_ChargeDigitizeBinning{this, "ChargeDigitizeBinning", 1, "Binning scheme used for charge discretization"};
    Gaudi::Property<double> m_chargeMax{this, "ChargeMaximum", (double)15000., "Cell dynamic range in electrons"};
    Gaudi::Property<bool>   m_DigitizeTime{this, "DigitizeTime", true, "Flag to enable digitization of timing information."};
    Gaudi::Property<int>    m_TimeDigitizeNumBits{this, "TimeDigitizeNumBits", 10, "Number of bits used to determine bins for time discretization"};
    Gaudi::Property<double> m_timeMax{this, "TimeMaximum", (double)10., "Cell dynamic range for timing measurement [ns]"};
    Gaudi::Property<int>    m_TimeDigitizeBinning{this, "TimeDigitizeBinning", 0, "Binning scheme used for time discretization"};
    Gaudi::Property<double> m_timeSmearingSigma{this, "TimeSmearingSigma", 0.05, "Effective intrinsic time measurement resolution effects [ns]."};
    Gaudi::Property<bool>   m_electronicEffects{this, "ElectronicEffects", true, "Apply Electronic Effects"};
    Gaudi::Property<bool>   m_produceFullPattern{this, "StoreFiredPixels", false, "Store fired pixels"};

    MyG4UniversalFluctuationForSi *m_fluctuate;

    // charge discretization
    std::vector<double> m_DigitizedBins{};
    
    // geometry
    int m_numberOfLayers;
    std::vector<int>   m_laddersInLayer{};
#ifdef ZSEGMENTED
    std::vector<int>   m_sensorsPerLadder{};
#endif
    std::vector<float> m_layerRadius{};
    std::vector<float> m_layerThickness{};
    std::vector<float> m_layerHalfThickness{};
    std::vector<float> m_layerLadderLength{};
    std::vector<float> m_layerLadderHalfWidth{};
    std::vector<float> m_layerPhiOffset{};
    std::vector<float> m_layerActiveSiOffset{};
    std::vector<float> m_layerHalfPhi{};
    std::vector<float> m_layerLadderWidth{};
    // endcap specific
    std::vector<float> m_layerPetalLength{};
    std::vector<float> m_petalsInLayer{};
    std::vector<float> m_layerPetalInnerWidth{};
    std::vector<float> m_layerPetalOuterWidth{};
    const dd4hep::rec::SurfaceMap* m_map;

    /* Charge digitization helpers */
    void ProduceIonisationPoints(edm4hep::SimTrackerHit &hit, InternalState *intState) const;
    void ProduceSignalPoints(InternalState *intState) const;
    void ProduceHits(MutableSimTrackerHitVec &simTrkVec, edm4hep::SimTrackerHit &simHit, InternalState *intState) const;
    void PoissonSmearer(MutableSimTrackerHitVec &simTrkVec) const;
    void GainSmearer(MutableSimTrackerHitVec &simTrkVec) const;
    void ApplyThreshold(MutableSimTrackerHitVec &simTrkVec) const;
    void ChargeDigitizer(MutableSimTrackerHitVec &simTrkVec) const;

    /* Time digitization helpers */
    void TimeSmearer(MutableSimTrackerHitVec &simTrkVec) const;
    void TimeDigitizer(MutableSimTrackerHitVec &simTrkVec) const;

    /* Reconstruction of measurement and helpers */
    edm4hep::MutableTrackerHitPlane *ReconstructTrackerHit(MutableSimTrackerHitVec &simTrkVec, edm4hep::TrackerHitPlaneCollection *THcol, InternalState *intState) const;
    void TransformToLab(const int cellID, edm4hep::Vector3d xLoc, edm4hep::Vector3d xLab) const;
    void FindLocalPosition(edm4hep::SimTrackerHit &hit, edm4hep::Vector3d &localPosition, edm4hep::Vector3d &localDirection, InternalState *intState) const;
    void TransformXYToCellID(double x, double y, int & ix, int & iy, InternalState *intState) const;
    void TransformCellIDToXY(int ix, int iy, double & x, double & y, InternalState* instate) const;
    int GetPixelsInaRow(InternalState *intState) const;
    int GetPixelsInaColumn(InternalState *intState) const;

    StatusCode LoadGeometry();
    void PrintGeometryInfo();
    double randomTail( const double qmin, const double qmax ) const;

    /* Message Helpers */
    bool msgLevel(MSG::Level level) const{ return msgSvc()->outputLevel(name()) <= level; };

};
#endif
