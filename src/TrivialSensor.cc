#include "TrivialSensor.h"
#include "FindUnionAlgorithm.h"

TrivialSensor::TrivialSensor(int layer,
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
                            float t_step) :
    AbstractSensor(layer,
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
                   starttime,
                   t_step),
    pixels(),
    charged_pix(0),
    clock_cnt(0)
{
}

TrivialSensor::~TrivialSensor() {}

void TrivialSensor::Reset()
{
    pixels.assign(l_rows * l_columns, 0);
    charged_pix = 0;
}

void TrivialSensor::BeginClockStep()
{
    pixels.assign(l_rows * l_columns, 0);
    charged_pix = 0;
}

void TrivialSensor::UpdatePixel(int x, int y, float chrg)
{
    LinearPosition lpos = l_locate(x, y);
    if (pixels[lpos] == 0) charged_pix++;
    pixels[lpos] += chrg;
}

void TrivialSensor::EndClockStep()
{
    clock_cnt += 1;
}

PixelData TrivialSensor::GetPixel(int x, int y)
{
    PixelData result { 0, 0, PixelStatus::off };
    LinearPosition lpos = l_locate(x, y);

    if (pixels[lpos] != 0)
    {
        result.charge = pixels[lpos];
        result.time = init_time + clock_cnt * clock_step;
        result.status = PixelStatus::on;
    }

    return result;
}

bool TrivialSensor::IsActive()
{
    return charged_pix > 0;
}

bool TrivialSensor::CheckStatus(int x, int y, PixelStatus pstat)
{
    LinearPosition lpos = l_locate(x, y);
    bool result = pixels[lpos] == 0 and pstat == PixelStatus::off;
    result |= pixels[lpos] != 0 and pstat == PixelStatus::on;
    return result;
}

void TrivialSensor::buildHits(SegmentDigiHitList& output)
{
    FindUnionAlgorithm  fu_algo { s_rows, s_colums };
    BitField64 bf_encoder = getBFEncoder();

    for (int h = 0; h < GetSegNumX(); h++)
    {
        for (int k = 0; k < GetSegNumY(); k++)
        {
            //Sensor segments ordered row first
            LinearPosition sens_id = s_locate(h, k);
            bf_encoder[LCTrackerCellID::sensor()] = sens_id;

            fu_algo.init();

            for (int i = 0; i < GetSensorRows(); i++)
            {
                for (int j = 0; j < GetSensorCols(); j++)
                {
                    if (!checkStatus(h, k, i, j, PixelStatus::on))
                    {
                        fu_algo.invalidate(i, j);
                        continue;
                    }

                    bool up_is_above = checkStatus(h, k, i - 1, j, PixelStatus::on);
                    bool left_is_above = checkStatus(h, k, i, j - 1, PixelStatus::on);

                    if (up_is_above && left_is_above)
                    {
                        fu_algo.merge(i - 1, j, i, j - 1);
                        fu_algo.merge(i, j - 1, i, j);
                    }
                    else if (!up_is_above && left_is_above)
                    {
                        fu_algo.merge(i, j - 1, i, j);
                    }
                    else if (up_is_above && !left_is_above)
                    {
                        fu_algo.merge(i - 1, j, i, j);
                    }
                }
            }

            fu_algo.close();

            for (ClusterOfPixel c_item : fu_algo.get_clusters())
            {
                // Very simple implementation: geometric mean
                float x_acc = 0;
                float y_acc = 0;
                float tot_charge = 0;

                for (LinearPosition curr_pos : c_item)
                {
                    GridCoordinate gcoor = l_locate(curr_pos);
                    PixelData p_data = GetPixel(gcoor.row, gcoor.col);
                    x_acc += PixelRowToX(gcoor.row);
                    y_acc += PixelColToY(gcoor.col);
                    tot_charge += p_data.charge;
                    
                }

                SegmentDigiHit digiHit = {
                    x_acc / c_item.size(),
                    y_acc / c_item.size(),
                    tot_charge,
                    init_time + clock_cnt * clock_step,
                    bf_encoder.lowWord()
                };
                output.push_back(digiHit);
            }
        }
    }
}

