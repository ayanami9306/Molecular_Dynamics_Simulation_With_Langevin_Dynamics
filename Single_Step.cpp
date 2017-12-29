//
//  Single_Step.cpp
//  Molecular Dynamics
//
//  Created by Nozomi on 2017. 11. 14..
//  Copyright © 2017년 JiHoon. All rights reserved.
//

#include "Molecular_Model.hpp"

void Model_Segment::Single_Step()
{
    time_Now += deltaT;
    Velocity_Verlet_Step();
    Compute_Forces();
}
