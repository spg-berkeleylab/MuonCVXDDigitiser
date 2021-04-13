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
                                 float s_level,
								 int q_level):
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
    _satur_level(s_level),
	_q_level(q_level),
    max_charge(0),
    charge_valid(false)
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
    pixels = EnergyMatrix(l_rows * l_columns);
}

PixelDigiMatrix::~PixelDigiMatrix()
{}

void PixelDigiMatrix::Reset()
{
    pixels.assign(pixels.size(), {0, 0, PixelStatus::undefined});

    charge_valid = false;
}

void PixelDigiMatrix::SetTime(float time)
{
    for (long unsigned int k = 0; k < pixels.size(); k++) pixels[k].time = time;
}

void PixelDigiMatrix::UpdatePixel(int x, int y, float chrg)
{
    if (check(x, y))
    {
        // linear aggregation with threshold
        int idx = index(x, y);
        float new_chrg = pixels[idx].charge + chrg;
        pixels[idx].charge = new_chrg > _satur_level ? _satur_level : new_chrg;
    }

    charge_valid = false;
}

void PixelDigiMatrix::Apply(PixelTransformation l_expr)
{
    for (long unsigned int k = 0; k < pixels.size(); k++)
    {
        PixelData tmpd = l_expr(pixels[k]);
        tmpd.charge = std::min(tmpd.charge, _satur_level);
        pixels[k] = tmpd;
    }

    charge_valid = false;
}

PixelData PixelDigiMatrix::GetPixel(int x, int y)
{
    if (status != MatrixStatus::ok)
    {
        return { 0, 0, PixelStatus::geometry_error };
    }
    if (check(x, y))
    {
        return pixels[index(x, y)];
    }
    return { 0, 0, PixelStatus::out_of_bounds };
}

float PixelDigiMatrix::GetMaxCharge()
{
    if (charge_valid) return max_charge;

    max_charge = 0;
    for (long unsigned int k = 0; k < pixels.size(); k++)
    {
        float t_chrg = pixels[k].charge;
        if (t_chrg > max_charge) max_charge = t_chrg;
    }
    charge_valid = true;

    return max_charge;
}

PixelData PixelDigiMatrix::GetPixel(int seg_x, int seg_y, int pos_x, int pos_y)
{
    return GetPixel(SensorRowToLadderRow(seg_x, pos_x), SensorColToLadderCol(seg_y, pos_y));
}

bool PixelDigiMatrix::check(int x, int y)
{
    if (x < 0 || x >= l_rows) return false;
    if (y < 0 || y >= l_columns) return false;
    return true;
}

