#include "DetElemSlidingWindow.h"
#include <EVENT/MCParticle.h>
#include "DDRec/DetectorData.h"
#include "DD4hep/DD4hepUnits.h"
#include "marlin/VerbosityLevels.h"
#include "streamlog/streamlog.h"

#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_math.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandFlat.h"

using std::max;
using std::min;
using std::vector;
using dd4hep::rec::ISurface;
using dd4hep::rec::Vector2D;
using dd4hep::rec::Vector3D;
using CLHEP::RandGauss;
using CLHEP::RandPoisson;
using CLHEP::RandFlat;

DetElemSlidingWindow::DetElemSlidingWindow(HitTemporalIndexes& htable,
                                           PixelDigiMatrix& sensor,
                                           float tclick,
                                           float wsize,
                                           float starttime,
                                           double tanLorentzAngleX,
                                           double tanLorentzAngleY,
                                           double cutOnDeltaRays,
                                           double diffusionCoefficient,
                                           double electronsPerKeV,
                                           double segmentLength,
                                           double energyLoss,
                                           double widthOfCluster,
                                           double electronicNoise,
                                           double maxTrkLen,
                                           double maxEnergyDelta,
                                           const SurfaceMap* s_map):
    //curr_time(htable.GetMinTime() - wrad - tclick),
    curr_time(starttime),
    time_click(tclick),
    window_radius(wsize / 2),
    _htable(htable),
    _sensor(sensor),
    _tanLorentzAngleX(tanLorentzAngleX),
    _tanLorentzAngleY(tanLorentzAngleY),
    _cutOnDeltaRays(cutOnDeltaRays),
    _diffusionCoefficient(diffusionCoefficient),
    _electronsPerKeV(electronsPerKeV),
    _segmentLength(segmentLength),
    _energyLoss(energyLoss),
    _widthOfCluster(widthOfCluster),
    _electronicNoise(electronicNoise),
    _maxTrkLen(maxTrkLen),
    _deltaEne(maxEnergyDelta),
    signals(),
    surf_map(s_map),
    cell_decoder(sensor.GetCellIDFormatStr())
{
    _fluctuate = new G4UniversalFluctuation();
}

DetElemSlidingWindow::~DetElemSlidingWindow()
{
    delete(_fluctuate);
}


bool DetElemSlidingWindow::move_forward()
{
    curr_time += time_click;
    streamlog_out(DEBUG) << "Time window centered in " << curr_time << std::endl;

    streamlog_out(DEBUG) << "Hits available for " << _sensor.GetLayer() << ":" << _sensor.GetLadder()
                           << " = " << _htable.GetHitNumber(_sensor.GetLayer(), _sensor.GetLadder())
                           << std::endl;

    for (SimTrackerHit* hit = _htable.CurrentHit(_sensor.GetLayer(), _sensor.GetLadder());
         hit != nullptr && hit->getTime() - curr_time < window_radius;
         hit = _htable.CurrentHit(_sensor.GetLayer(), _sensor.GetLadder()))
    {
        StoreSignalPoints(hit);
        _htable.DisposeHit(_sensor.GetLayer(), _sensor.GetLadder());
    }

    if (!signals.empty())
    {
        streamlog_out(DEBUG) << "Signal points for " << _sensor.GetLayer() << ":" << _sensor.GetLadder()
                               << " = " << signals.size() << std::endl;

        for (TimedSignalPoint spoint = signals.front();
             curr_time - spoint.time > window_radius;
             spoint = signals.front())
        {
            signals.pop_front();
            if (signals.empty()) break;
        }
    }

    UpdatePixels();

    return _htable.GetHitNumber(_sensor.GetLayer(), _sensor.GetLadder()) > 0;
}

void DetElemSlidingWindow::UpdatePixels()
{
    _sensor.Reset();

    for (auto spoint : signals)
    {
        if (spoint.time > curr_time + window_radius) break;

        double xHFrame = _widthOfCluster * spoint.sigmaX;
        double yHFrame = _widthOfCluster * spoint.sigmaY;
        
        int ixLo, ixUp, iyLo, iyUp;
        _sensor.TransformXYToCellID(spoint.x - xHFrame, spoint.y - yHFrame, ixLo, iyLo);
        _sensor.TransformXYToCellID(spoint.x + xHFrame, spoint.y + yHFrame, ixUp, iyUp);

        for (int ix = ixLo; ix < ixUp + 1; ++ix)
        {
            if (ix < 0 || ix >= _sensor.GetSizeX()) continue;

            for (int iy = iyLo; iy < iyUp + 1; ++iy)
            {
                if (iy < 0 || iy >= _sensor.GetSizeY()) continue;

                double xCurrent, yCurrent;
                _sensor.TransformCellIDToXY(ix, iy, xCurrent, yCurrent);
                
                gsl_sf_result result;
                int status = gsl_sf_erf_Q_e((xCurrent - 0.5 * _sensor.GetPixelSizeX() - spoint.x) / spoint.sigmaX, &result);
                double LowerBound = 1 - result.val;

                status = gsl_sf_erf_Q_e((xCurrent + 0.5 * _sensor.GetPixelSizeX() - spoint.x) / spoint.sigmaX, &result);
                double UpperBound = 1 - result.val;
                double integralX = UpperBound - LowerBound;

                status = gsl_sf_erf_Q_e((yCurrent - 0.5 * _sensor.GetPixelSizeY() - spoint.y) / spoint.sigmaY, &result);
                LowerBound = 1 - result.val;

                status = gsl_sf_erf_Q_e((yCurrent + 0.5 * _sensor.GetPixelSizeY() - spoint.y) / spoint.sigmaY, &result);
                UpperBound = 1 - result.val;
                double integralY = UpperBound - LowerBound;

                float totCharge = float(spoint.charge * integralX * integralY);

                _sensor.UpdatePixel(ix, iy, { totCharge, spoint.time, PixelStatus::ok });
            }
        }
    }

    _sensor.Apply([](PixelData data) -> PixelData {
        float new_charge = 0;
        if (data.charge > 1e+03) // assume Gaussian
        {
            new_charge += float(RandGauss::shoot(data.charge, sqrt(data.charge)));
        }
        else // assume Poisson
        {
            new_charge += float(RandPoisson::shoot(data.charge));
        }
        return { new_charge, data.time, PixelStatus::ok };
    });

    if (_electronicNoise > 0)
    {
        _sensor.Apply([&, this](PixelData data) -> PixelData {
            return {
                data.charge + float(RandGauss::shoot(0., this->_electronicNoise)), 
                data.time,
                PixelStatus::ok
            };
        });
    }
}

void DetElemSlidingWindow::StoreSignalPoints(SimTrackerHit* hit)
{
    // hit and pos are in mm
    double pos[3] = {0,0,0};
    double dir[3] = {0,0,0};
    double entry[3];
    double exit[3];

    // ************************* Find local position **************************
    SurfaceMap::const_iterator sI = surf_map->find(hit->getCellID0()) ;
    const ISurface* surf = sI->second ;

    Vector3D oldPos( hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2] );
    // We need it?
    if ( ! surf->insideBounds( dd4hep::mm * oldPos ) ) {

        streamlog_out( DEBUG ) << "  hit at " << oldPos
                                << " is not on surface "
                                << *surf
                                << " distance: " << surf->distance(  dd4hep::mm * oldPos )
                                << std::endl;
        return;
    }    


    Vector2D lv = surf->globalToLocal( dd4hep::mm * oldPos  ) ;
    // Store local position in mm
    pos[0] = lv[0] / dd4hep::mm;
#ifdef ZSEGMENTED
    int segment_id = cell_decoder(hit)["sensor"];

    float s_offset = _sensor.GetSegSizeY() * _sensor.GetPixelSizeY() * (float(segment_id) + 0.5);
    s_offset -= _sensor.GetHalfLength();
    pos[1] = (lv[1] + s_offset)/ dd4hep::mm;
#else
    pos[1] = lv[1] / dd4hep::mm;
#endif

    // Add also z ccordinate
    Vector3D origin( surf->origin()[0], surf->origin()[1], surf->origin()[2]);
    pos[2] = ( dd4hep::mm * oldPos - dd4hep::cm * origin ).dot( surf->normal() ) / dd4hep::mm;

    double Momentum[3];
    for (int j = 0; j < 3; ++j) 
      if (hit->getMCParticle())
        Momentum[j] = hit->getMCParticle()->getMomentum()[j] * dd4hep::GeV;
      else
        Momentum[j] = hit->getMomentum()[j];

    // as default put electron's mass
    double particleMass = 0.510e-3 * dd4hep::GeV;
    if (hit->getMCParticle())
        particleMass = max(hit->getMCParticle()->getMass() * dd4hep::GeV, particleMass);

    double particleMomentum = sqrt(pow(Momentum[0], 2) + pow(Momentum[1], 2) + pow(Momentum[2], 2));                   
                         
    dir[0] = Momentum * surf->u();
    dir[1] = Momentum * surf->v();
    dir[2] = Momentum * surf->normal();

    // ************************************************************************

    entry[2] = -_sensor.GetHalfThickness(); 
    exit[2] = _sensor.GetHalfThickness();

    for (int i = 0; i < 2; ++i) {
        entry[i] = pos[i] + dir[i] * (entry[2] - pos[2]) / dir[2];
        exit[i]= pos[i] + dir[i] * (exit[2] - pos[2]) / dir[2];
    }

    double tanx = dir[0] / dir[2];
    double tany = dir[1] / dir[2];  
    
    // trackLength is in mm
    double trackLength = min(_maxTrkLen, _sensor.GetThickness() * sqrt(1.0 + pow(tanx, 2) + pow(tany, 2)));
  
    int _numberOfSegments = ceil(trackLength / _segmentLength );
    double dEmean = (dd4hep::keV * _energyLoss * trackLength) / ((double)_numberOfSegments);

    // TODO _segmentLength may be different from segmentLength, is it ok?
    double segmentLength = trackLength / ((double)_numberOfSegments);
    double _segmentDepth = _sensor.GetThickness() / ((double)_numberOfSegments);

    double z = -_sensor.GetHalfThickness() - 0.5 * _segmentDepth;

    double eSum = 0.0;
    vector<TimedSignalPoint> signal_buffer{_numberOfSegments};

    for (int i = 0; i < _numberOfSegments; ++i)
    {
        // ionization point
        z += _segmentDepth;
        double x = pos[0] + tanx * (z - pos[2]);
        double y = pos[1] + tany * (z - pos[2]);
        // momentum in MeV/c, mass in MeV, tmax (delta cut) in MeV, 
        // length in mm, meanLoss eloss in MeV.
        double eloss = _fluctuate->SampleFluctuations(particleMomentum * dd4hep::keV / dd4hep::MeV,
                                                      particleMass * dd4hep::keV / dd4hep::MeV,
                                                      _cutOnDeltaRays,
                                                      segmentLength,
                                                      dEmean / dd4hep::MeV) * dd4hep::MeV;

        double DistanceToPlane = _sensor.GetHalfThickness() - z;
        double xOnPlane = x + _tanLorentzAngleX * DistanceToPlane;
        double yOnPlane = y + _tanLorentzAngleY * DistanceToPlane;
        double DriftLength = DistanceToPlane * sqrt(1.0 + pow(_tanLorentzAngleX, 2) 
                                                        + pow(_tanLorentzAngleY, 2));

        double SigmaDiff = sqrt(DriftLength / _sensor.GetThickness()) * _diffusionCoefficient;

        double SigmaX = SigmaDiff * sqrt(1.0 + pow(_tanLorentzAngleX, 2));
        double SigmaY = SigmaDiff * sqrt(1.0 + pow(_tanLorentzAngleY, 2));

        // energy is in keV       
        double charge = (eloss / dd4hep::keV) * _electronsPerKeV;

        signal_buffer.push_back({
            xOnPlane,
            yOnPlane,
            SigmaX,
            SigmaY,
            charge,
            hit->getTime()
        });

        eSum += eloss;
    }

    double hEdep = hit->getEDep() / dd4hep::GeV;
    // deltaEne is a charge??
    const double thr = _deltaEne / _electronsPerKeV * dd4hep::keV;
    while (hEdep > eSum + thr) {
      // Add additional charge sampled from an 1 / n^2 distribution.
      const double       q = randomTail( thr, hEdep - eSum );
      const unsigned int h = floor(RandFlat::shoot(0.0, (double)_numberOfSegments));
      signal_buffer[h].charge += q * _electronsPerKeV / dd4hep::keV;
      eSum += q;
    }

    for(auto spoint : signal_buffer)
    {
        signals.push_back(spoint);
    }
}

//=============================================================================
// Sample charge from 1 / n^2 distribution.
//=============================================================================
double DetElemSlidingWindow::randomTail( const double qmin, const double qmax )
{
    const double offset = 1. / qmax;
    const double range  = ( 1. / qmin ) - offset;
    const double u      = offset + RandFlat::shoot() * range;
    return 1. / u;
}


