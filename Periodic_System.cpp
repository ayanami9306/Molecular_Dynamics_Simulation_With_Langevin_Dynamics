//
//  Periodic_System.cpp
//  Molecular Dynamics
//
//  Created by Nozomi on 2017. 11. 14..
//  Copyright © 2017년 JiHoon. All rights reserved.
//

#include "Molecular_Model.hpp"

void Model_Segment::Periodic_Length(double * data)
{
    if(data[0] > BOUNDARY_SIZE_X/2) data[0] -= BOUNDARY_SIZE_X;
    else if(data[0] < -BOUNDARY_SIZE_X/2) data[0] += BOUNDARY_SIZE_X;
    if(data[1] > BOUNDARY_SIZE_Y/2) data[1] -= BOUNDARY_SIZE_Y;
    else if(data[1] < -BOUNDARY_SIZE_Y/2) data[1] += BOUNDARY_SIZE_Y;
    if(data[2] > BOUNDARY_SIZE_Z/2) data[2] -= BOUNDARY_SIZE_Z;
    else if(data[2] < -BOUNDARY_SIZE_Z/2) data[2] += BOUNDARY_SIZE_Z;
}

void Model_Segment::Periodic_Boundary(double * data)
{
    if(data[0] >= BOUNDARY_SIZE_X) data[0] -= BOUNDARY_SIZE_X;
    else if(data[0] < 0) data[0] += BOUNDARY_SIZE_X;
    if(data[1] >= BOUNDARY_SIZE_Y) data[1] -= BOUNDARY_SIZE_Y;
    else if(data[1] < 0) data[1] += BOUNDARY_SIZE_Y;
    if(data[2] >= BOUNDARY_SIZE_Z) data[2] -= BOUNDARY_SIZE_Z;
    else if(data[2] < 0) data[2] += BOUNDARY_SIZE_Z;
}
