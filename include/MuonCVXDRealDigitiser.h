#ifndef MuonCVXDRealDigitiser_h
#define MuonCVXDRealDigitiser_h 1

#include <string>
#include <vector>

// k4FWCore
#include <k4FWCore/Transformer.h>

// edm4hep
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/MutableTrackerHit.h>
#include <edm4hep/MutableSimTrackerHit.h>
#include <edm4hep/MutableTrackerHitPlane.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <edm4hep/TrackerHitSimTrackerHitLink.h>
#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>

// DD4hep
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

// ROOT
#include <TH1.h>

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

typedef std::vector<MutableSimTrackerHit> MutableSimTrackerHitVec;
typedef std::vector<IonisationPoint> IonisationPointVec;
typedef std::vector<SignalPoint> SignalPointVec;

/**  Digitizer for Simulated Hits in the Vertex Detector. <br>
 * Digitization follows the procedure adopted in the CMS software package. 
 * See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePixelDigitization
 * 
 * @param CollectionName name of input SimTrackerHit collection <br>
 * (default parameter value : "VXDCollection")
 * @param OutputCollectionName name of output TrackerHitsPlane collection <br>
 * (default parameter value : "VTXTrackerHits")
 * @param RelationColName name of the LCRelation <br>
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
 * @param PoissonSmearing flag to switch on gaussian smearing of electrons collected on pixels <br>
 * (default parameter value : 1) <br>
 * @param ElectronicEffects flag to switch on gaussian smearing of signal (electronic noise) <br>
 * (default parameter value : 1) <br>
 * @param ElectronicNoise electronic noise in electrons <br>
 * (default parameter value : 100) <br>
 * @param StoreFiredPixels flag to store also the fired pixels (collection names: "VTXPixels") <br>
 * (default parameter value : 0) <br>
 * @param EnergyLoss Energy loss in keV/mm <br>
 * (default parameter value : 280.0) <br>
 * @param MaxEnergyDelta max delta in energy hit (difference from the hit simulated charged and the comuputed ones) in electrons <br>
 * (default parameter value : 100) <br>
  * @param MaxTrackLength Maximum values for track path length inside the ladder (in mm)", <br>
 * (default parameter value : 10) <br> 
 * <br>
 */
class MuonCVXDRealDigitiser : public k4FWCore::MultiTransformer<std::tuple<edm4hep::SimTrackerHitCollection,
                                                                           edm4hep::TrackerHitPlaneCollection,
                                                                           edm4hep::TrackerHitSimTrackerHitCollection,
                                                                           edm4hep::TrackerHitSimTrackerHitCollection>(
                                                                     const edm4hep::SimTrackerHitCollection&)>
{  
public:
  
    MuonCVXDRealDigitiser(const std::string& name, ISvcLocator* svcLoc);

    /** Called at the begin of the job before anything is read.
    * Use to initialize the processor, e.g. book histograms.
    */
    virtual StatusCode initialize();

    /** Called for every event - the working horse.
    */
    virtual std::tuple<edm4hep::SimTrackerHitCollection,
                       edm4hep::TrackerHitPlaneCollection,
                       edm4hep::TrackerHitSimTrackerHitCollection,
                       edm4hep::TrackerHitSimTrackerHitCollection> operator(
                 const edm4hep::SimTrackerHitCollection& STHcol) const;

    /** Called after data processing for clean up.
    */
    virtual StatusCode finalize();

protected:
    void LoadGeometry() const;
    void PrintGeometryInfo() const;

    int m_barrelID = 0;

    // processor parameters
    Gaudi::Property<std::string> m_subDetName{this, "SubDetectorName", std::string("VertexBarrel"), "Name of Vertex detector"};
    Gaudi::Property<double> m_tanLorentzAngleX{this, "TanLorentz", (double)0.8, "Tangent of Lorentz Angle"};
    Gaudi::Property<double> m_tanLorentzAngleY{this, "TanLorentzY", (double)0, "Tangent of Lorentz Angle along Y"};
    Gaudi::Property<double> m_cutOnDeltaRays{this, "CutOnDeltaRays", (double)0.030, "Cut on delta-ray energy (MeV)"};
    Gaudi::Property<double> m_diffusionCoefficient{this, "DiffusionCoefficient", (double)0.07, "Diffusion coefficient, sqrt(D / mu / V)."};
    Gaudi::Property<double> m_pixelSizeX{this, "PixelSizeX", (double)0.025, "Pixel Size X"};
    Gaudi::Property<double> m_pixelSizeY{this, "PixelSizeY", (double)0.025, "Pixel Size Y"};
    Gaudi::Property<double> m_electronsPerKeV{this, "ElectronsPerKeV", (double)270.3, "Electrons per keV"};
    Gaudi::Property<double> m_segmentLength{this, "SegmentLength", (double)0.005, "Segment Length in mm"};
    Gaudi::Property<double> m_threshold{this, "Threshold", (double)200., "Cell Threshold in electrons"};
    Gaudi::Property<double> m_electronicNoise{this, "ElectronicNoise", (double)100., "electronic noise in electrons"};
    Gaudi::Property<double> m_energyLoss{this, "EnergyLoss", (double)280., "Energy Loss keV/mm"};
    Gaudi::Property<double> m_deltaEne{this, "MaxEnergyDelta", (double)100., "Max delta in energy between G4 prediction and random sampling for each hit in electrons"};
    Gaudi::Property<double> m_maxTrkLen{this, "MaxTrackLength", (double)10., "Maximum values for track length (in mm)"};

    Gaudi::Property<float> m_window_size{this, "WindowSize", (float)25, "Window size [ns]"};
    Gaudi::Property<float> m_fe_slope{this, "RD53Aslope", (float)0.1, "ADC slope for chip RD53A"};
    Gaudi::Property<int>    m_PoissonSmearing{this, "PoissonSmearing", 1, "Apply Poisson smearing of electrons collected on pixels"};
    Gaudi::Property<int>    m_electronicEffects{this, "ElectronicEffects", 1, "Apply Electronic Effects"};
    Gaudi::Property<int>    m_produceFullPattern{this, "StoreFiredPixels", 0, "Store fired pixels"};
    Gaudi::Property<int>    m_sensor_type{this, "SensorType", 1, "Sensor model to be used (0 : ChipRD53A, 1 : Trivial)"};

    // geometry
    int m_numberOfLayers;
    std::vector<int>   m_laddersInLayer{};
    std::vector<int>   m_sensorsPerLadder{};
    std::vector<float> m_layerRadius{};
    std::vector<float> m_layerThickness{};
    std::vector<float> m_layerHalfThickness{};
    std::vector<float> m_layerLadderLength{};
    std::vector<float> m_layerLadderHalfWidth{};
    std::vector<float> m_layerPhiOffset{};
    std::vector<float> m_layerActiveSiOffset{};
    std::vector<float> m_layerHalfPhi{};
    std::vector<float> m_layerLadderWidth{};
    const dd4hep::rec::SurfaceMap* m_map ;

    // Graphs
    Gaudi::Property<std::string> m_stat_filename{this, "StatisticsFilename", std::string("None"), "File name for statistics (None for disabling the feature)"};
    bool  create_stats;
    TH1F* signal_dHisto;
    TH1F* bib_dHisto;
    TH1F* signal_cSizeHisto;
    TH1F* signal_xSizeHisto;
    TH1F* signal_ySizeHisto;
    TH1F* signal_zSizeHisto;
    TH1F* signal_eDepHisto;
    TH1F* bib_cSizeHisto;
    TH1F* bib_xSizeHisto;
    TH1F* bib_ySizeHisto;
    TH1F* bib_zSizeHisto;
    TH1F* bib_eDepHisto;

};

#endif //MuonCVXDRealDigitiser_h
