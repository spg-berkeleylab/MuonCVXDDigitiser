#include "PixelDigiMatrix.h"
#include <cmath>

PixelDigiMatrix::PixelDigiMatrix(int layer,
                                 int ladder,
                                 int xsegmentNumber,
                                 int ysegmentNumber,
                                 float ladderLength,
                                 float ladderWidth,
                                 float thickness,
                                 double pixelSizeX,
                                 double pixelSizeY):
    _layer(layer),
    _ladder(ladder),
    _thickness(thickness),
    _pixelSizeX(abs(pixelSizeX)),
    _pixelSizeY(abs(pixelSizeY)),
    _ladderLength(ladderLength > 0 ? ladderLength : 0),
    _ladderWidth(ladderWidth > 0 ? ladderWidth : 0)
{
    x_size = ceil(ladderWidth / pixelSizeX);
    y_size = ceil(ladderLength / pixelSizeY);

    x_segsize = xsegmentNumber > 0 ? x_size / xsegmentNumber : 1;
    y_segsize = ysegmentNumber > 0 ? y_size / ysegmentNumber : 1;

    if (fmod(ladderWidth, pixelSizeX) != 0 || fmod(ladderLength, pixelSizeY) != 0)
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
}

void PixelDigiMatrix::UpdatePixel(int x, int y, PixelData data)
{
    if (check(x, y))
    {
        pixels[index(x, y)] = data;
    }
}

void PixelDigiMatrix::Apply(PixelTransformation l_expr)
{
    for (long unsigned int k = 0; k < pixels.size(); k++)
    {
        pixels[k] = l_expr(pixels[k]);
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
        return pixels[index(x, y)];
    }
    return { 0, 0, PixelStatus::out_of_bounds };
}

void PixelDigiMatrix::TransformXYToCellID(double x, double y, int & ix, int & iy)
{
    double yInLadder = y + _ladderLength / 2;
    iy = int(yInLadder / _pixelSizeY);

    double xInLadder = x + _ladderWidth / 2;
    ix = int(xInLadder / _pixelSizeX);
}

void PixelDigiMatrix::TransformCellIDToXY(int ix, int iy, double & x, double & y)
{
    y = ((0.5 + double(iy)) * _pixelSizeY) - _ladderLength / 2;
    x = ((0.5 + double(ix)) * _pixelSizeX) - _ladderWidth / 2;
}

bool PixelDigiMatrix::check(int x, int y)
{
    if (x < 0 || x >= x_size) return false;
    if (y < 0 || y >= y_size) return false;
    return true;
}

