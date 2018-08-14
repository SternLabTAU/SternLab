/*
 FITS - Flexible Inference from Time-Series data
 (c) 2016-2018 by Tal Zinger
 tal.zinger@outlook.com
 
 ActualDataFile_Entry.cpp
 
 Stores information regarding the current generation of the population.
 Swap and operator implementations used for sorting, copying and moving objects.
 */

#include "ActualDataFile.hpp"


void ActualDataEntry::swap(ActualDataEntry& other)
{
    std::swap(pos, other.pos);
    std::swap(gen, other.gen);
    std::swap(base, other.base);
    std::swap(freq, other.freq);
    std::swap(ref, other.ref);
    std::swap(read_count, other.read_count);
}


// used for sorting only
bool ActualDataEntry::operator<( const ActualDataEntry& other ) const
{
    if (gen < other.gen) {
        return true;
    }
    
    if (gen > other.gen) {
        return false;
    }
    
    // we reached here, so gen == other.gen
    return base < other.base;
}


ActualDataEntry& ActualDataEntry::operator=(ActualDataEntry other)
{
    swap(other);
    return *this;
}


