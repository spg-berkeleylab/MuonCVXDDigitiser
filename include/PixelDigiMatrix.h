#ifndef PixelDigiMatrix_h
#define PixelDigiMatrix_h 1

#include <string>
#include <vector>
#include <unordered_map>

using std::string;
using std::unordered_map;
using std::unordered_multimap;
using std::vector;

enum class PixelStatus : char {
    on,
    off,
    ready,
    start,
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
using ClockTicks = int;

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

struct LocatedPixel
{
    int row;
    int col;
    PixelData data;
};

/**
 * @class PixelDigiMatrix
 * @brief Simulation of the chip RD53A
 *
 * This class implements the basic behaviour of the chip [RD53A](https://cds.cern.ch/record/2113263).
 * It simulates a matrix of pixels. Each pixel can collect a charge and perform a linear charge
 * depletion mechanism in order to measure the charge itself. A single threshold detection, with smearing,
 * is available for each pixel.
 * The matrix of pixels corresponds to the entire ladder; the ladder is divided into a grid of sensors.
 * Each matrix is identify by a couple of ID: the layer ID and the ladder ID.
 * This class must be operated by an agent which feeds it with the charge and synchronize the actions
 * through a clock.
 */
class PixelDigiMatrix
{
public:
    /**
     * @brief The main constructor
     *
     * This constructor creates a matrix of pixels for a given ladder
     * within a layer of the vertex barrel.
     * @param layer The ID of the layer containing the pixel matrix
     * @param ladder The ID of the ladder matching this matrix of pixels
     * @param xsegmentNumber The number of sensors per ladder width
     * @param ysegmentNumber The number of sensors per ladder length
     * @param ladderLength The length of the ladder
     * @param ladderWidth The width of the ladder
     * @param thickness The thickness the ladder
     * @param pixelSizeX The width of a pixel
     * @param pixelSizeY The length of a pixel
     * @param enc_str The format string used to encode the CellID for any sensor of a ladder
     * @param barrel_id The ID of the vertex barrel inside the detector
     * @param thr The threshold for any pixel of the ladder
     * @param fe_slope The charge depletion slope of the FE
     * @param starttime The start time for the matrix evolution
     * @param t_step The clock period of the chip
     */
    PixelDigiMatrix(int layer,
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
                    float t_step);
    virtual ~PixelDigiMatrix();

    virtual void buildHits(SegmentDigiHitList& output) = 0;

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
    inline MatrixStatus GetStatus() { return status; }

    inline string GetCellIDFormatStr() { return cellFmtStr; }

    void Reset();

    void BeginClockStep();

    /**
     * @brief The charge aggregation call.
     *
     * This method must be called by the agent when a quantity of charge must be gathered
     * for a given pixel of the ladder
     * @param x The row number of the pixel
     * @param y The column number of the pixel
     * @param chrg The charge to be aggregated in the pixel
     */
    void UpdatePixel(int x, int y, float chrg);

    void EndClockStep();

    PixelData GetPixel(int x, int y);
    bool IsActive();
    bool CheckStatus(int x, int y, PixelStatus pstat);

    inline int XToPixelRow(double x) { return int((x + _ladderWidth / 2) / _pixelSizeX); }
    inline int YToPixelCol(double y) { return int((y + _ladderLength / 2) / _pixelSizeY); }

    inline double PixelRowToX(int ix) { return ((0.5 + double(ix)) * _pixelSizeX) - _ladderWidth / 2; }
    inline double PixelColToY(int iy) { return ((0.5 + double(iy)) * _pixelSizeY) - _ladderLength / 2; }

protected:

    inline int SensorRowToLadderRow(int seg_x, int pos_x) { return seg_x * s_rows + pos_x; }
    inline int SensorColToLadderCol(int seg_y, int pos_y) { return seg_y * s_colums + pos_y; }
    inline int LadderRowToSensorRow(int pos_x, int seg_x) { return pos_x - seg_x * s_rows; }
    inline int LadderColToSensorCol(int pos_y, int seg_y) { return pos_y - seg_y * s_colums; }

    PixelData GetPixel(int seg_x, int seg_y, int pos_x, int pos_y);
    bool CheckStatus(int seg_x, int seg_y, int pos_x, int pos_y, PixelStatus pstat);
    bool CheckStatusOnSensor(int seg_x, int seg_y, PixelStatus pstat);

    vector<LocatedPixel> GetPixelsFromSensor(int seg_x, int seg_y, PixelStatus pstat);

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
    float delta_c;
    GridPosition l_locate;
    GridPosition s_locate;

private:

    struct PixelRawData
    {
        float charge;
        ClockTicks t_begin;
        ClockTicks t_end;
    };

    using SensorBin = unordered_multimap<LinearPosition, LinearPosition>;

    inline bool check(int x, int y) { return (0 <= x and x < l_rows) and (0 <= y || y < l_columns); }
    PixelStatus calc_status(LinearPosition lpos);
    ClockTicks calc_end_clock(float charge);
    LinearPosition sensor_for_pixel(LinearPosition pos);

    MatrixStatus status;
    unordered_map<LinearPosition, PixelRawData> pixels;
    unordered_map<ClockTicks, SensorBin> expir_table;
    unordered_map<LinearPosition, float> charge_buffer;
    SensorBin start_table; 
};

#endif //PixelDigiMatrix_h
