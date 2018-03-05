#include "Molecular_Model.hpp"

void Model_Segment::Measure_Radius_of_Gyration()
{
    //calculate radius of gyration
    
    //sidechain
    int num_sidechain = (int)(dp_backbone/space_sidechain);
    int len_sidechain = (nParticle-dp_backbone)/num_sidechain;
    for(int i_sidechain= 0; i_sidechain < num_sidechain; i_sidechain++)
    {
        double R_CM_side[3] =
        {
            0, 0, 0
        };
        
        for(int i = i_sidechain*len_sidechain + dp_backbone; i < (i_sidechain+1)*len_sidechain + dp_backbone; i++)
        {
            R_CM_side[0] += Segment[i].coordinate[0];
            R_CM_side[1] += Segment[i].coordinate[1];
            R_CM_side[2] += Segment[i].coordinate[2];
        }
            
        R_CM_side[0] /= len_sidechain;
        R_CM_side[1] /= len_sidechain;
        R_CM_side[2] /= len_sidechain;
            
        double Temp_Rx_side = 0, Temp_Ry_side = 0, Temp_Rz_side = 0;
            
        for(int i = i_sidechain*len_sidechain + dp_backbone; i < (i_sidechain+1)*len_sidechain + dp_backbone; i++)
        {
            double Value_side[3]
            {
                pow((Segment[i].coordinate[0] - R_CM_side[0]), 2.0),
                pow((Segment[i].coordinate[1] - R_CM_side[1]), 2.0),
                pow((Segment[i].coordinate[2] - R_CM_side[2]), 2.0)
            };
            if(Value_side[0] <= Value_side[1])
            {
                if(Value_side[1] <= Value_side[2])
                {
                    Temp_Rx_side += Value_side[0];
                    Temp_Ry_side += Value_side[1];
                    Temp_Rz_side += Value_side[2];
                }
                else if(Value_side[0] <= Value_side[2])
                {
                    Temp_Rx_side += Value_side[0];
                    Temp_Ry_side += Value_side[2];
                    Temp_Rz_side += Value_side[1];
                }
                else
                {
                    Temp_Rx_side += Value_side[2];
                    Temp_Ry_side += Value_side[0];
                    Temp_Rz_side += Value_side[1];
                }
            }
            else if(Value_side[0] <= Value_side[2])
            {
                Temp_Rx_side += Value_side[1];
                Temp_Ry_side += Value_side[0];
                Temp_Rz_side += Value_side[2];
            }
            else if(Value_side[1] <= Value_side[2])
            {
                Temp_Rx_side += Value_side[1];
                Temp_Ry_side += Value_side[2];
                Temp_Rz_side += Value_side[0];
            }
            else
            {
                Temp_Rx_side += Value_side[2];
                Temp_Ry_side += Value_side[1];
                Temp_Rz_side += Value_side[0];
            }
        }
        
        Rx_sidechain += sqrt(Temp_Rx_side/len_sidechain)/num_sidechain;
        Ry_sidechain += sqrt(Temp_Ry_side/len_sidechain)/num_sidechain;
        Rz_sidechain += sqrt(Temp_Rz_side/len_sidechain)/num_sidechain;
        Rg_sidechain += sqrt((Temp_Rx_side+Temp_Ry_side+Temp_Rz_side)/len_sidechain)/num_sidechain;
        
    }
    
    //backbone
    double R_CM_backbone[3] =
    {
        0, 0, 0
    };
    
    for(int i=0;i<dp_backbone;i++)
    {
        R_CM_backbone[0] += Segment[i].coordinate[0];
        R_CM_backbone[1] += Segment[i].coordinate[1];
        R_CM_backbone[2] += Segment[i].coordinate[2];
    }
    R_CM_backbone[0] /= dp_backbone;
    R_CM_backbone[1] /= dp_backbone;
    R_CM_backbone[2] /= dp_backbone;
    
    double Temp_Rx_backbone = 0, Temp_Ry_backbone = 0, Temp_Rz_backbone = 0;
    
    for(int i = 0;i<dp_backbone;i++)
    {
        double Value_backbone[3]
        {
            pow((Segment[i].coordinate[0] - R_CM_backbone[0]), 2.0),
            pow((Segment[i].coordinate[1] - R_CM_backbone[1]), 2.0),
            pow((Segment[i].coordinate[2] - R_CM_backbone[2]), 2.0)
        };
        if(Value_backbone[0] <= Value_backbone[1])
        {
            if(Value_backbone[1] <= Value_backbone[2])
            {
                Temp_Rx_backbone += Value_backbone[0];
                Temp_Ry_backbone += Value_backbone[1];
                Temp_Rz_backbone += Value_backbone[2];
            }
            else if(Value_backbone[0] <= Value_backbone[2])
            {
                Temp_Rx_backbone += Value_backbone[0];
                Temp_Ry_backbone += Value_backbone[2];
                Temp_Rz_backbone += Value_backbone[1];
            }
            else
            {
                Temp_Rx_backbone += Value_backbone[2];
                Temp_Ry_backbone += Value_backbone[0];
                Temp_Rz_backbone += Value_backbone[1];
            }
        }
        else if(Value_backbone[0] <= Value_backbone[2])
        {
            Temp_Rx_backbone += Value_backbone[1];
            Temp_Ry_backbone += Value_backbone[0];
            Temp_Rz_backbone += Value_backbone[2];
        }
        else if(Value_backbone[1] <= Value_backbone[2])
        {
            Temp_Rx_backbone += Value_backbone[1];
            Temp_Ry_backbone += Value_backbone[2];
            Temp_Rz_backbone += Value_backbone[0];
        }
        else
        {
            Temp_Rx_backbone += Value_backbone[2];
            Temp_Ry_backbone += Value_backbone[1];
            Temp_Rz_backbone += Value_backbone[0];
        }
    }
    
    Rx_backbone += sqrt(Temp_Rx_backbone/dp_backbone);
    Ry_backbone += sqrt(Temp_Ry_backbone/dp_backbone);
    Rz_backbone += sqrt(Temp_Rz_backbone/dp_backbone);
    Rg_backbone += sqrt((Temp_Rx_backbone+Temp_Ry_backbone+Temp_Rz_backbone)/dp_backbone);
    
    //total
    double R_CM_total[3] =
    {
        0, 0, 0
    };
    
    for(int i=0;i<nParticle;i++)
    {
        R_CM_total[0] += Segment[i].coordinate[0];
        R_CM_total[1] += Segment[i].coordinate[1];
        R_CM_total[2] += Segment[i].coordinate[2];
    }
    R_CM_total[0] /= nParticle;
    R_CM_total[1] /= nParticle;
    R_CM_total[2] /= nParticle;
    
    double Temp_Rx_total = 0, Temp_Ry_total = 0, Temp_Rz_total = 0;
    
    for(int i = 0;i<nParticle;i++)
    {
        double Value_total[3]
        {
            pow((Segment[i].coordinate[0] - R_CM_total[0]), 2.0),
            pow((Segment[i].coordinate[1] - R_CM_total[1]), 2.0),
            pow((Segment[i].coordinate[2] - R_CM_total[2]), 2.0)
        };
        if(Value_total[0] <= Value_total[1])
        {
            if(Value_total[1] <= Value_total[2])
            {
                Temp_Rx_total += Value_total[0];
                Temp_Ry_total += Value_total[1];
                Temp_Rz_total += Value_total[2];
            }
            else if(Value_total[0] <= Value_total[2])
            {
                Temp_Rx_total += Value_total[0];
                Temp_Ry_total += Value_total[2];
                Temp_Rz_total += Value_total[1];
            }
            else
            {
                Temp_Rx_total += Value_total[2];
                Temp_Ry_total += Value_total[0];
                Temp_Rz_total += Value_total[1];
            }
        }
        else if(Value_total[0] <= Value_total[2])
        {
            Temp_Rx_total += Value_total[1];
            Temp_Ry_total += Value_total[0];
            Temp_Rz_total += Value_total[2];
        }
        else if(Value_total[1] <= Value_total[2])
        {
            Temp_Rx_total += Value_total[1];
            Temp_Ry_total += Value_total[2];
            Temp_Rz_total += Value_total[0];
        }
        else
        {
            Temp_Rx_total += Value_total[2];
            Temp_Ry_total += Value_total[1];
            Temp_Rz_total += Value_total[0];
        }
    }
    
    Rx_total += sqrt(Temp_Rx_total/nParticle);
    Ry_total += sqrt(Temp_Ry_total/nParticle);
    Rz_total += sqrt(Temp_Rz_total/nParticle);
    Rg_total += sqrt((Temp_Rx_total+Temp_Ry_total+Temp_Rz_total)/nParticle);
}

