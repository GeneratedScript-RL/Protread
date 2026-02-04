#ifndef PROTREAD_HPP
#define PROTREAD_HPP

#include <vector>

struct ContactPair
{
    int i;
    int j;
};

// Main contact calculation function
std::vector<ContactPair> calculate_map_sparse(
    const double *coords,
    int n,
    double threshold);

#endif // PROTREAD_HPP