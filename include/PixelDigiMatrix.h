#ifndef PixelDigiMatrix_h
#define PixelDigiMatrix_h 1

#include <vector>
#include <functional>

enum class PixelStatus : char {
    ok,
    undefined,
    out_of_bounds,
    geometry_error
};

struct PixelData
{
    float charge;
    float time;
    PixelStatus status;
};

enum class MatrixStatus : char {
    ok,
    pixel_number_error,
    segment_number_error
};

typedef std::vector<PixelData> EnergyMatrix;
typedef std::function<PixelData(PixelData pIn)> PixelTransformation;

class PixelDigiMatrix
{
public:
    PixelDigiMatrix(int layer,
                    int ladder,
                    int xsegmentNumber,
                    int ysegmentNumber,
                    float ladderLength,
                    float ladderWidth,
                    float thickness,
                    double pixelSizeX,
                    double pixelSizeY);
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
    inline MatrixStatus GetStatus() { return status; }

    void Reset();
    void UpdatePixel(int x, int y, PixelData data);
    void Apply(PixelTransformation l_expr);
    PixelData GetPixel(int x, int y);

    void TransformXYToCellID(double x, double y, int & ix, int & iy);
    void TransformCellIDToXY(int ix, int iy, double & x, double & y);

private:
    inline int index(int x, int y) { return x * x_size + y; }
    bool check(int x, int y);

    int _layer;
    int _ladder;
    float _thickness;
    double _pixelSizeX;
    double _pixelSizeY;
    float _ladderLength;
    float _ladderWidth;
    int x_size;
    int y_size;
    int x_segsize;
    int y_segsize;
    EnergyMatrix pixels;
    MatrixStatus status;
};

#endif //PixelDigiMatrix_h
