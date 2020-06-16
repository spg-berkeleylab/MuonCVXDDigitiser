#ifndef PixelDigiMatrix_h
#define PixelDigiMatrix_h 1

#include <vector>

struct PixelData
{
    float edep;
    float time;
};

typedef std::vector<PixelData> EnergyMatrix;

class PixelDigiMatrix
{
public:
    PixelDigiMatrix(int layer,
                    int ladder,
                    float thickness,
                    double pixelSizeX,
                    double pixelSizeY,
                    float ladderLength,
                    float ladderWidth);
    virtual ~PixelDigiMatrix();

    inline int GetLayer() { return _layer; }
    inline int GetLadder() { return _ladder; }
    inline float GetThickness() { return _thickness; }
    inline float GetHalfThickness() { return _thickness / 2; }
    inline float GetLength() { return _ladderLength; }
    inline float GetHalfLength() { return _ladderLength / 2; }
    inline float GetWidth() { return _ladderWidth; }
    inline float GetHalfWidth() { return _ladderWidth / 2; }
    inline double GetPixelSizeX() { return _pixelSizeX; }
    inline double GetPixelSizeY() { return _pixelSizeY; }
    inline int GetSizeX() { return x_size; }
    inline int GetSizeY() { return y_size; }

    void Reset();
    void UpdatePixel(int x, int y, PixelData data);
    PixelData GetPixel(int x, int y);
    void TransformXYToCellID(double x, double y, int & ix, int & iy);
    void TransformCellIDToXY(int ix, int iy, double & x, double & y);

private:
    int _layer;
    int _ladder;
    float _thickness;
    double _pixelSizeX;
    double _pixelSizeY;
    float _ladderLength;
    float _ladderWidth;
    int x_size;
    int y_size;
    EnergyMatrix pixels;
};

#endif //PixelDigiMatrix_h
