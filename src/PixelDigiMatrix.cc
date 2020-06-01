#include "PixelDigiMatrix.h"

PixelDigiMatrix::PixelDigiMatrix(int x, int y):
    x_size(x),
    y_size(y),
    pixels(x > 0 && y > 0 ? x * y : 0)
{}

PixelDigiMatrix::PixelDigiMatrix(PixelDigiMatrix&& pdm)
{}

PixelDigiMatrix::~PixelDigiMatrix()
{}

PixelDigiMatrix& PixelDigiMatrix::operator=(PixelDigiMatrix&& pdm)
{
    return *this;
}

PixelDigiMatrix::operator bool() const
{
    return x_size > 0 && y_size > 0;
}


