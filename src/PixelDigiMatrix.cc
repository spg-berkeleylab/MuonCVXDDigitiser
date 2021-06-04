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
    delta_c(t_step * fe_slope),
    _active(false)
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
    _active = false;
    for (long unsigned int k = 0; k < pixels.size(); k++)
    {
        if (pixels[k].charge == 0) continue;

        if (pixels[k].active)
        {
            pixels[k].counter += 1;
        }
        else
        {
            pixels[k].counter = 0;
        }

        pixels[k].active = IsOverThreshold(pixels[k].charge);
        _active = _active || pixels[k].active;

        pixels[k].charge = std::max(pixels[k].charge - delta_c, 0.f);
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
        auto pix = pixels[index(x, y)];
        PixelData result { 0, 0, PixelStatus::off };

        if (pix.active)
        {
            result.time = clock_time - pix.counter * clock_step;
            result.status = pix.counter == 0 ? PixelStatus::start : PixelStatus::on;
        }
        else if (pix.counter > 0)
        {
            result.charge += pix.counter * delta_c;
            result.time = clock_time - pix.counter * clock_step;
            result.status = PixelStatus::ready;
        }
        return result;
    }
    return { 0, 0, PixelStatus::out_of_bounds };
}

PixelData PixelDigiMatrix::GetPixel(int seg_x, int seg_y, int pos_x, int pos_y)
{
    return GetPixel(SensorRowToLadderRow(seg_x, pos_x), SensorColToLadderCol(seg_y, pos_y));
}

bool PixelDigiMatrix::IsActive()
{
    return _active;
}

bool PixelDigiMatrix::IsOverThreshold(float charge)
{
    // TODO implement noises and threshold dispersion effects
    return charge > _thr_level;
}