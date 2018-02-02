#include "Molecular_Model.hpp"

double Model_Segment::Measure_Radius_of_Gyration()
{
    //calculate radius of gyration
    double R_CM[3] =
    {
        0, 0, 0
    };
    
    for(int i=0;i<nParticle;i++)
    {
        R_CM[0] += Segment[i].coordinate[0];
        R_CM[1] += Segment[i].coordinate[1];
        R_CM[2] += Segment[i].coordinate[2];
    }
    R_CM[0] /= nParticle;
    R_CM[1] /= nParticle;
    R_CM[2] /= nParticle;
    
    double Rg = 0;
    
    for(int i = 0;i<nParticle;i++)
    {
        Rg +=
        pow((Segment[i].coordinate[0] - R_CM[0]), 2.0) +
        pow((Segment[i].coordinate[1] - R_CM[1]), 2.0) +
        pow((Segment[i].coordinate[2] - R_CM[2]), 2.0);
    }
    return sqrt(Rg/nParticle);
}

