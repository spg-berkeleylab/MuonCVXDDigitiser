#ifndef AbstractSensor_h
#define AbstractSensor_h 1

#include <string>
#include <vector>

using std::string;
using std::vector;

enum class PixelStatus : char {
    on,
    off,
    ready,
    start,
    out_of_bounds,
    geometry_error
};

enum class MatrixStatus : char {
    ok,
    pixel_number_error,
    segment_number_error
};

struct PixelData
{
    float charge;
    float time;
    PixelStatus status;
};

struct SegmentDigiHit
{
    float x;
    float y;
    float charge;
    float time;
    int cellID0;
};

typedef vector<SegmentDigiHit> SegmentDigiHitList;

struct GridCoordinate
{
    int row;
    int col;
};

static inline bool operator==(GridCoordinate a, GridCoordinate b)
{
    return a.row == b.row && a.col == b.col;
}

using LinearPosition = int;

class GridPosition
{
public:
    GridPosition(int rows, int cols) : b_size(cols) {}
    virtual ~GridPosition() {}
    LinearPosition operator()(int row, int col) { return row * b_size + col; }
    GridCoordinate operator()(LinearPosition pos)
    {
        return { pos / b_size, pos % b_size };
    }
private:
    int b_size;
};

/**
 * @class AbstractSensor
 * @brief Abstract class for segmented sensor
 *
 * Abstract class for segmented sensor.
 */
class AbstractSensor
{
public:

    AbstractSensor( int layer,
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
                    float starttime,
                    float t_step);

    ~AbstractSensor();

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
    inline int GetLadderRows() { return l_rows; }
    inline int GetLadderCols() { return l_columns; }
    inline int GetSensorRows() { return s_rows; }
    inline int GetSensorCols() { return s_colums; }
    inline int GetSegNumX() { return x_segnum; }
    inline int GetSegNumY() { return y_segnum; }

    inline int XToPixelRow(double x) { return int((x + _ladderWidth / 2) / _pixelSizeX); }
    inline int YToPixelCol(double y) { return int((y + _ladderLength / 2) / _pixelSizeY); }

    inline double PixelRowToX(int ix) { return ((0.5 + double(ix)) * _pixelSizeX) - _ladderWidth / 2; }
    inline double PixelColToY(int iy) { return ((0.5 + double(iy)) * _pixelSizeY) - _ladderLength / 2; }

    inline string GetCellIDFormatStr() { return cellFmtStr; }

    inline MatrixStatus GetStatus() { return status; }

    virtual void Reset() = 0;

    virtual void BeginClockStep() = 0;

    virtual void UpdatePixel(int x, int y, float chrg) = 0;

    virtual void EndClockStep() = 0;

    virtual PixelData GetPixel(int x, int y) = 0;

    virtual bool IsActive() = 0;

    virtual bool CheckStatus(int x, int y, PixelStatus pstat) = 0;

    virtual void buildHits(SegmentDigiHitList& output) = 0;

protected:

    inline int SensorRowToLadderRow(int seg_x, int pos_x) { return seg_x * s_rows + pos_x; }
    inline int SensorColToLadderCol(int seg_y, int pos_y) { return seg_y * s_colums + pos_y; }
    inline int LadderRowToSensorRow(int pos_x, int seg_x) { return pos_x - seg_x * s_rows; }
    inline int LadderColToSensorCol(int pos_y, int seg_y) { return pos_y - seg_y * s_colums; }

    int _barrel_id;
    int _layer;
    int _ladder;
    float _thickness;
    double _pixelSizeX;
    double _pixelSizeY;
    float _ladderLength;
    float _ladderWidth;
    int l_rows;
    int l_columns;
    int s_rows;
    int s_colums;
    int x_segnum;
    int y_segnum;
    string cellFmtStr;
    double _thr_level;
    float init_time;
    int clock_cnt;
    float clock_step;
    GridPosition l_locate;
    GridPosition s_locate;
    MatrixStatus status;
};


#endif //AbstractSensor_h