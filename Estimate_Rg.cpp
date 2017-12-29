//
//  Estimate_Rg.cpp
//  Molecular Dynamics
//
//  Created by Nozomi on 2017. 12. 14..
//  Copyright © 2017년 JiHoon. All rights reserved.
//

#include "Molecular_Model.hpp"

double Model_Segment::CALCULATE_RADIUS_OF_GYRATION()
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
  
//end to end vector
  /*  double R[3] =
    {
        Segment[0].coordinate[0] - Segment[nParticle-1].coordinate[0],
        Segment[0].coordinate[1] - Segment[nParticle-1].coordinate[1],
        Segment[0].coordinate[2] - Segment[nParticle-1].coordinate[2],
    };
    return sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);*/
}
//
