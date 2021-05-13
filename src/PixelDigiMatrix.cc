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
    clock_time(0),
    q_slope(0)   // TODO read from config
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
    pixels.assign(pixels.size(), {0, 0});

}

PixelDigiMatrix::~PixelDigiMatrix()
{}

void PixelDigiMatrix::Reset()
{
    pixels.assign(pixels.size(), {0, 0});
}

void PixelDigiMatrix::ClockSync(float time)
{
    float delta_c = (time - clock_time) * q_slope;

    for (long unsigned int k = 0; k < pixels.size(); k++)
    {
        if (pixels[k].charge == 0) continue;

        pixels[k].charge = std::max(pixels[k].charge - delta_c, 0.f);
        if (pixels[k].charge < _thr_level)
        {
            // TODO store the TOT
        }
        else
        {
            pixels[k].counter += 1;
        }
    }

    clock_time = time;
}

void PixelDigiMatrix::UpdatePixel(int x, int y, float chrg)
{
    if (check(x, y))
    {
        int idx = index(x, y);
        float new_chrg = pixels[idx].charge + chrg;
        pixels[idx].charge = new_chrg;
        if (new_chrg > _thr_level) pixels[idx].counter = 0;
    }
}
/*
void PixelDigiMatrix::Apply(PixelTransformation l_expr)
{
    for (long unsigned int k = 0; k < pixels.size(); k++)
    {
        PixelData tmpd = l_expr(pixels[k]);
        tmpd.charge = std::min(tmpd.charge, _satur_level);
        pixels[k] = tmpd;
    }
}
*/
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

bool PixelDigiMatrix::check(int x, int y)
{
    if (x < 0 || x >= l_rows) return false;
    if (y < 0 || y >= l_columns) return false;
    return true;
}

