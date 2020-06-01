#ifndef PixelDigiMatrix_h
#define PixelDigiMatrix_h 1

#include <vector>

typedef std::vector<float> EnergyMatrix;

class DetElemSlidingWindow;

class PixelDigiMatrix
{
    friend class DetElemSlidingWindow;
public:
    PixelDigiMatrix(int x = 0, int y = 0);
    PixelDigiMatrix(PixelDigiMatrix&& pdm);
    virtual ~PixelDigiMatrix();
    PixelDigiMatrix& operator=(PixelDigiMatrix&& pdm);
    explicit operator bool() const;
private:
    int x_size;
    int y_size;
    EnergyMatrix pixels;
};

#endif //PixelDigiMatrix_h
