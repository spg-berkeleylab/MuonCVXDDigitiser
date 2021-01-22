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
                                 int barrel_id):
    _barrel_id(barrel_id),
    _layer(layer),
    _ladder(ladder),
    _thickness(thickness),
    _pixelSizeX(abs(pixelSizeX)),
    _pixelSizeY(abs(pixelSizeY)),
    _ladderLength(ladderLength > 0 ? ladderLength : 0),
    _ladderWidth(ladderWidth > 0 ? ladderWidth : 0),
    cellFmtStr(enc_str),
    max_charge(0),
    charge_valid(false)
{
    int lwid = floor(ladderWidth * 1e4);
    int psx = floor(pixelSizeX * 1e4);
    int llen = floor(ladderLength * 1e4);
    int psy = floor(pixelSizeY * 1e4);

    x_size = lwid / psx;
    y_size = llen / psy;

    x_segnum = xsegmentNumber;
    y_segnum = ysegmentNumber;

    x_segsize = xsegmentNumber > 0 ? x_size / xsegmentNumber : 1;
    y_segsize = ysegmentNumber > 0 ? y_size / ysegmentNumber : 1;
    
    if (lwid % psx > 0 || llen % psy > 0)
    {
        status = MatrixStatus::pixel_number_error;
    }
    else if (x_size % xsegmentNumber > 0 || y_size % ysegmentNumber > 0)
    {
        status = MatrixStatus::segment_number_error;
    }
    else
    {
        status = MatrixStatus::ok;
    }
    pixels = EnergyMatrix(x_size * y_size);
}

PixelDigiMatrix::~PixelDigiMatrix()
{}

void PixelDigiMatrix::Reset()
{
    pixels.assign(pixels.size(), {0, 0, PixelStatus::undefined});

    charge_valid = false;
}

void PixelDigiMatrix::UpdatePixel(int x, int y, PixelData data)
{
    if (check(x, y))
    {
        pixels[index(x, y)] = data;
    }

    charge_valid = false;
}

void PixelDigiMatrix::Apply(PixelTransformation l_expr)
{
    for (long unsigned int k = 0; k < pixels.size(); k++)
    {
        pixels[k] = l_expr(pixels[k]);
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
    if (x < 0 || x >= x_size) return false;
    if (y < 0 || y >= y_size) return false;
    return true;
}

