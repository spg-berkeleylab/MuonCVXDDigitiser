#include "PixelDigiMatrix.h"
#include <cmath>
#include "streamlog/streamlog.h"

PixelDigiMatrix::PixelDigiMatrix(int layer,
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
                                 float fe_slope,
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
    clock_time(starttime),
    clock_step(t_step),
    q_slope(fe_slope)
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
    else
    {
        status = MatrixStatus::ok;
    }

    pixels = std::vector<PixelRawData>(l_rows * l_columns);
    pixels.assign(pixels.size(), {0, 0, false});

}

PixelDigiMatrix::~PixelDigiMatrix()
{}

void PixelDigiMatrix::Reset()
{
    pixels.assign(pixels.size(), {0, 0, false});
}

void PixelDigiMatrix::ClockSync()
{
    float delta_c = clock_step * q_slope;

    for (long unsigned int k = 0; k < pixels.size(); k++)
    {
        if (pixels[k].charge == 0) continue;

        pixels[k].charge = std::max(pixels[k].charge - delta_c, 0.f);
        pixels[k].counter += 1;
        pixels[k].thr_down = (pixels[k].charge < getThreshold()  and not pixels[k].thr_down);
    }

    clock_time += clock_step;
}

void PixelDigiMatrix::UpdatePixel(int x, int y, float chrg)
{
    if (check(x, y))
    {
        int idx = index(x, y);
        float new_chrg = pixels[idx].charge + chrg;

        pixels[idx].charge = new_chrg;
        if (new_chrg > getThreshold() and not pixels[idx].thr_down) pixels[idx].counter = 0;
    }
}

PixelData PixelDigiMatrix::GetPixel(int x, int y)
{
    if (status != MatrixStatus::ok)
    {
        return { 0, 0, PixelStatus::geometry_error };
    }
    if (check(x, y))
    {
        // TODO implement
        //return pixels[index(x, y)];
    }
    return { 0, 0, PixelStatus::out_of_bounds };
}

PixelData PixelDigiMatrix::GetPixel(int seg_x, int seg_y, int pos_x, int pos_y)
{
    return GetPixel(SensorRowToLadderRow(seg_x, pos_x), SensorColToLadderCol(seg_y, pos_y));
}

double PixelDigiMatrix::getThreshold()
{
    // TODO implement noises and threshold dispersion effects
    return _thr_level;
}
