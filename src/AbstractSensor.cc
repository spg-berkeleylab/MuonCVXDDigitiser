#include "AbstractSensor.h"
#include <cmath>

AbstractSensor::AbstractSensor( int layer,
                                int ladder,
                                int xsegmentNumber,
                                int ysegmentNumber,
                                float ladderLength,
                                float ladderWidth,
                                float thickness,
                                double pixelSizeX,
                                double pixelSizeY,
                                string enc_str,
                                int barrel_id,
                                double thr,
                                float starttime,
                                float t_step):
    _barrel_id(barrel_id),
    _layer(layer),
    _ladder(ladder),
    _thickness(thickness),
    _pixelSizeX(fabs(pixelSizeX)),
    _pixelSizeY(fabs(pixelSizeY)),
    _ladderLength(ladderLength > 0 ? ladderLength : 0),
    _ladderWidth(ladderWidth > 0 ? ladderWidth : 0),
    cellFmtStr(enc_str),
    _thr_level(thr),
    init_time(starttime),
    clock_cnt(0),
    clock_step(t_step),
    l_locate({ 0, 0 }),
    s_locate({ 0, 0 }),
    status(MatrixStatus::ok)
{
    int lwid = floor(ladderWidth * 1e4);
    int psx = floor(pixelSizeX * 1e4);
    int llen = floor(ladderLength * 1e4);
    int psy = floor(pixelSizeY * 1e4);

    l_rows = lwid / psx;
    l_columns = llen / psy;

    x_segnum = xsegmentNumber;
    y_segnum = ysegmentNumber;

    s_rows = xsegmentNumber > 0 ? l_rows / xsegmentNumber : 1;
    s_colums = ysegmentNumber > 0 ? l_columns / ysegmentNumber : 1;
    
    if (lwid % psx > 0 || llen % psy > 0)
    {
        status = MatrixStatus::pixel_number_error;
    }
    else if (l_rows % xsegmentNumber > 0 || l_columns % ysegmentNumber > 0)
    {
        status = MatrixStatus::segment_number_error;
    }

    l_locate = GridPosition(l_rows, l_columns);
    s_locate = GridPosition(x_segnum, y_segnum);
}

AbstractSensor::~AbstractSensor()
{}