//
//  Molecular_Model.cpp
//  Molecular Dynamics
//

#include "Molecular_Model.hpp"

Model_Segment::Model_Segment()
{
    inv_RAND_MAX = 1.0 / RAND_MAX;
    srand((unsigned int)time(NULL));
}

Model_Segment::~Model_Segment()
{
    
}
