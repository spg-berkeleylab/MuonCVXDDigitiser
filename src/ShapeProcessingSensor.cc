#include "ShapeProcessingSensor.h"

#include "streamlog/streamlog.h"

ShapeProcessingSensor::ShapeProcessingSensor(int layer,
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
                                             int q_level) :
    HKBaseSensor(layer,
                 ladder,
                 xsegmentNumber,
                 ysegmentNumber,
                 ladderLength,
                 ladderWidth,
                 thickness,
                 pixelSizeX,
                 pixelSizeY,
                 enc_str,
                 barrel_id,
                 thr,
                 s_level,
                 q_level)
{}

/* ***************************************************************************************
 * 
 * Radial sweep algorithm
 * TODO missing holes detection in cluster
 * 
 * ************************************************************************************ */
GridCoordinate ShapeProcessingSensor::GetNextPoint(GridCoordinate c, GridCoordinate p)
{
    int idx = (c.row - p.row) + (c.col - p.col) * 3;
    switch(idx)
    {
    case 1:
        return { c.row - 1, c.col + 1 };
    case -2:
        return { c.row,     c.col + 1 };
    case -3:
        return { c.row + 1, c.col + 1 };
    case -4:
        return { c.row + 1, c.col     };
    case -1:
        return { c.row + 1, c.col - 1 };
    case 2:
        return { c.row    , c.col - 1 };
    case 3:
        return { c.row - 1, c.col - 1 };
    case 4:
        return { c.row - 1, c.col     };
    }
    
    //unreachable statement
    return { 0, 0 };
}

vector<GridCoordinate> ShapeProcessingSensor::GetContour(ClusterOfPixel& spot)
{
    vector<GridCoordinate> result;

    // It's better to put an empty frame around the spot
    int b_rows = spot.row_max - spot.row_min + 3;
    int b_cols = spot.col_max - spot.col_min + 3;
    vector<bool> buffer { b_rows * b_cols, false };

    GridCoordinate curr_p { 0, 0 };
    GridCoordinate prev_p { 0, 0 };
    for (GridCoordinate a_coord : spot.pix)
    {
        int l_row = a_coord.row - spot.row_min + 1;
        int l_col = a_coord.col - spot.col_min + 1;
        buffer[l_row * b_cols + l_col] = true;

        // Starting points
        if (curr_p.row == 0 && a_coord.row == spot.row_min)
        {
            curr_p = { l_row, l_col };
            prev_p = { l_row - 1, l_col };
        }
    }

    // simple protection against infinite loop
    int avail_cycle = 8 * b_rows * b_cols;
    do
    {
        result.push_back(curr_p);
        avail_cycle--;
        
        if (avail_cycle == 0)
        {
            result.clear();
            break;
        }

        GridCoordinate next_p = GetNextPoint(curr_p, prev_p);
        while (!buffer[next_p.row * b_cols + next_p.col])
        {
            next_p = GetNextPoint(curr_p, next_p);
        }
        
        prev_p = curr_p;
        curr_p = next_p;
    }
    while (result.size() > 1 && result[0] == prev_p && result[1] == curr_p);
    
    //TODO last element is equal to the first

    return result;
}


/* ***************************************************************************************
 * 
 * Sub-clustering
 * 
 * ************************************************************************************ */
ClusterOfPixel ShapeProcessingSensor::processCluster(ClusterOfPixel in)
{
    if (in.pix.size() > 10)
    {
        if(streamlog::out.write<streamlog::MESSAGE>())
#pragma omp critical
        {
            streamlog::out() << "Find big cluster: " << in.row_min << ":" << in.col_min 
                             << " " << in.row_max << ":" << in.col_max << std::endl;
        }

    }

    return in;
}