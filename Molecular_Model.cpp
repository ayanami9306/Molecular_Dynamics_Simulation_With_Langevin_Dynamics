#include "Molecular_Model.hpp"

Model_Segment::Model_Segment()
{
    inv_RAND_MAX = 1.0 / RAND_MAX;
    srand((unsigned int)time(NULL));
    unsigned int seed[16];
    for(int i=0;i<16;i++) seed[i] = rand();
    InitWELLRNG512a(seed);
}

Model_Segment::~Model_Segment()
{
    
}

