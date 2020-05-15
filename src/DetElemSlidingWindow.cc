#include "DetElemSlidingWindow.h"
#include <EVENT/MCParticle.h>
#include "DDRec/DetectorData.h"
#include "DD4hep/DD4hepUnits.h"
#include "marlin/VerbosityLevels.h"

using std::max;
using std::min;
using dd4hep::rec::ISurface;
using dd4hep::rec::Vector2D;
using dd4hep::rec::Vector3D;

DetElemSlidingWindow::DetElemSlidingWindow(HitTemporalIndexes& htable,
                                           int layer,
                                           int ladder,
                                           float tclick,
                                           float wrad,
                                           float thickness,
                                           double tanLorentzAngleX,
                                           double tanLorentzAngleY,
                                           double cutOnDeltaRays,
                                           double diffusionCoefficient,
                                           double electronsPerKeV,
                                           double segmentLength,
                                           double energyLoss,
                                           SurfaceMap* s_map):
    curr_time(htable.GetMinTime() - wrad - tclick),
    time_click(tclick),
    window_radius(wrad),
    _htable(htable),
    _layer(layer),
    _ladder(ladder),
    _layerThickness(thickness),
    _layerHalfThickness(thickness * 0.5),
    _tanLorentzAngleX(tanLorentzAngleX),
    _tanLorentzAngleY(tanLorentzAngleY),
    _cutOnDeltaRays(cutOnDeltaRays),
    _diffusionCoefficient(diffusionCoefficient),
    _electronsPerKeV(electronsPerKeV),
    _segmentLength(segmentLength),
    _energyLoss(energyLoss),
    signals(),
    surf_map(s_map)
{
    _fluctuate = new G4UniversalFluctuation();
}

DetElemSlidingWindow::~DetElemSlidingWindow()
{
    delete(_fluctuate);
}


void DetElemSlidingWindow::move_forward()
{
    curr_time += time_click;

    for(SimTrackerHit* hit = _htable.CurrentHit(_layer, _ladder);
        hit != nullptr && hit->getTime() - curr_time < window_radius;
        hit = _htable.CurrentHit(_layer, _ladder))
    {
        StoreSignalPoints(hit);
        _htable.DisposeHit(_layer, _ladder);
    }

    if (signals.empty()) return;
    for(TimedSignalPoint spoint = signals.front();
        curr_time - spoint.time > window_radius;
        spoint = signals.front())
    {
        signals.pop_front();
        if (signals.empty()) break;
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

        streamlog_out( DEBUG3 ) << "  hit at " << oldPos
                                << " is not on surface "
                                << *surf
                                << " distance: " << surf->distance(  dd4hep::mm * oldPos )
                                << std::endl;
        return;
    }    
    
    
    Vector2D lv = surf->globalToLocal( dd4hep::mm * oldPos  ) ;
    // Store local position in mm
    pos[0] = lv[0] / dd4hep::mm ;
    pos[1] = lv[1] / dd4hep::mm ;

    // Add also z ccordinate
    Vector3D origin( surf->origin()[0], surf->origin()[1], surf->origin()[2]);
    pos[2] = ( dd4hep::mm * oldPos - dd4hep::cm * origin ).dot( surf->normal() ) / dd4hep::mm;

    double Momentum[3];
    for (int j = 0; j < 3; ++j) 
      if (hit->getMCParticle())
        Momentum[j] = hit->getMCParticle()->getMomentum()[j] * dd4hep::GeV / dd4hep::keV;
      else
        Momentum[j] = hit->getMomentum()[j];

    // as default put electron's mass
    double particleMass = 0.510e-3 * dd4hep::GeV / dd4hep::keV;
    if (hit->getMCParticle())
        particleMass = max(hit->getMCParticle()->getMass() * dd4hep::GeV / dd4hep::keV, particleMass);

    double particleMomentum = sqrt(pow(Momentum[0], 2) + pow(Momentum[1], 2) + pow(Momentum[2], 2));                   
                         
    dir[0] = Momentum * surf->u();
    dir[1] = Momentum * surf->v();
    dir[2] = Momentum * surf->normal();

    // ************************************************************************

    entry[2] = -_layerHalfThickness; 
    exit[2] = _layerHalfThickness;

    for (int i = 0; i < 2; ++i) {
        entry[i] = pos[i] + dir[i] * (entry[2] - pos[2]) / dir[2];
        exit[i]= pos[i] + dir[i] * (exit[2] - pos[2]) / dir[2];
    }

    double tanx = dir[0] / dir[2];
    double tany = dir[1] / dir[2];  
    
    // trackLength is in mm
    double trackLength = min(1.0e+3, _layerThickness * sqrt(1.0 + pow(tanx, 2) + pow(tany, 2)));
  
    int _numberOfSegments = ceil(trackLength / _segmentLength );
    double dEmean = (dd4hep::keV * _energyLoss * trackLength) / ((double)_numberOfSegments);

    // TODO _segmentLength may be different from segmentLength, is it ok?
    double segmentLength = trackLength / ((double)_numberOfSegments);
    double _segmentDepth = _layerThickness / ((double)_numberOfSegments);

    double z = -_layerHalfThickness - 0.5 * _segmentDepth; 	
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

        double DistanceToPlane = _layerHalfThickness - z;
        double xOnPlane = x + _tanLorentzAngleX * DistanceToPlane;
        double yOnPlane = y + _tanLorentzAngleY * DistanceToPlane;
        double DriftLength = DistanceToPlane * sqrt(1.0 + pow(_tanLorentzAngleX, 2) 
                                                        + pow(_tanLorentzAngleY, 2));

        double SigmaDiff = sqrt(DriftLength / _layerThickness) * _diffusionCoefficient;

        double SigmaX = SigmaDiff * sqrt(1.0 + pow(_tanLorentzAngleX, 2));
        double SigmaY = SigmaDiff * sqrt(1.0 + pow(_tanLorentzAngleY, 2));

        // energy is in keV       
        double charge = (eloss / dd4hep::keV) * _electronsPerKeV;

        TimedSignalPoint spoint {
            xOnPlane,
            yOnPlane,
            SigmaX,
            SigmaY,
            charge, // electrons x KeV
            hit->getTime()
        };
        signals.push_back(spoint);
    }
}
