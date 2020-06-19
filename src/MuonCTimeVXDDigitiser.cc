#include "MuonCTimeVXDDigitiser.h"
#include "DetElemSlidingWindow.h"

MuonCTimeVXDDigitiser aMuonCTimeVXDDigitiser;

MuonCTimeVXDDigitiser::MuonCTimeVXDDigitiser() :
    MuonCVXDDigitiser()
{
    registerProcessorParameter("TimeClick",
                               "Time step",
                               _tclick,
                               (float)1);     // TODO set a default value
    registerProcessorParameter("WindowSize",
                               "Window size",
                               _window_size,
                               (float)1);     // TODO set a default value
}

void MuonCTimeVXDDigitiser::processEvent(LCEvent * evt)
{ 
    LCCollection* STHcol = nullptr;
    try
    {
        STHcol = evt->getCollection(_colName);
    }
    catch( lcio::DataNotAvailableException ex )
    {
        streamlog_out(WARNING) << _colName << " collection not available" << std::endl;
        STHcol = nullptr;
    }

    if( STHcol == nullptr ) return;

    HitTemporalIndexes t_index{STHcol};
    float start_time = t_index.GetMinTime() - _window_size / 2 - _tclick;

    // TODO pragma omp here
    for (int layer = 0; layer < _numberOfLayers; layer++)
    {
        for (int ladder = 0; ladder < _laddersInLayer[layer]; ladder++)
        {
#ifdef ZSEGMENTED
            int nun_segment = _sensorsPerLadder[layer];
#else
            int nun_segment = 1;
#endif
            PixelDigiMatrix sensor {
                layer, ladder,
                1, nun_segment,
                _layerLadderLength[layer],
                _layerLadderWidth[layer],
                _layerThickness[layer],
                _pixelSizeX, _pixelSizeY 
            };

            if (sensor.GetStatus() == MatrixStatus::pixel_number_error)
            {
                // TODO log error
                continue;
            }
            if (sensor.GetStatus() == MatrixStatus::segment_number_error)
            {
                // TODO log error
                continue;
            }

            DetElemSlidingWindow t_window {
                t_index, sensor,
                _tclick, _window_size, start_time,
                _tanLorentzAngleX, _tanLorentzAngleY,
                _cutOnDeltaRays,
                _diffusionCoefficient,
                _electronsPerKeV,
                _segmentLength,
                _energyLoss,
                _widthOfCluster,
                _electronicNoise,
                _maxTrkLen,
                _deltaEne,
                _map
            };
            
            bool goon = true;
            {
                goon = t_window.move_forward();
            }
            while(goon);
        }
    }
}
