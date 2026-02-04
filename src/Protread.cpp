#include "Protread.hpp"
#include <cmath>
#include <vector>
#include <unordered_map>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

struct SpatialKey
{
    int x, y, z;

    bool operator==(const SpatialKey &other) const
    {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct SpatialHash
{
    std::size_t operator()(const SpatialKey &k) const
    {
        std::size_t h = 73856093;
        h = h * 31 + k.x;
        h = h * 31 + k.y;
        h = h * 31 + k.z;
        return h ^ (h >> 16);
    }
};

std::vector<ContactPair> calculate_map_sparse(const double *coords, int n, double threshold)
{
    const double thresh_sq = threshold * threshold;
    const double cell_size = threshold;
    const double inv_cell_size = 1.0 / cell_size;

    std::unordered_map<SpatialKey, std::vector<int>, SpatialHash> grid;
    grid.reserve(n / 4);

    for (int i = 0; i < n; ++i)
    {
        SpatialKey key{
            static_cast<int>(std::floor(coords[i * 3 + 0] * inv_cell_size)),
            static_cast<int>(std::floor(coords[i * 3 + 1] * inv_cell_size)),
            static_cast<int>(std::floor(coords[i * 3 + 2] * inv_cell_size))};
        grid[key].push_back(i);
    }

    for (auto &entry : grid)
    {
        entry.second.shrink_to_fit();
    }

    std::vector<ContactPair> all_contacts;
    all_contacts.reserve(n * 10);

#pragma omp parallel
    {
        std::vector<ContactPair> local_contacts;
        local_contacts.reserve(1000);

#pragma omp for schedule(dynamic, 16) nowait
        for (int i = 0; i < n; ++i)
        {
            const double ix = coords[i * 3 + 0];
            const double iy = coords[i * 3 + 1];
            const double iz = coords[i * 3 + 2];

            SpatialKey center{
                static_cast<int>(std::floor(ix * inv_cell_size)),
                static_cast<int>(std::floor(iy * inv_cell_size)),
                static_cast<int>(std::floor(iz * inv_cell_size))};

            for (int dz = -1; dz <= 1; ++dz)
            {
                for (int dy = -1; dy <= 1; ++dy)
                {
                    for (int dx = -1; dx <= 1; ++dx)
                    {
                        SpatialKey neighbor{
                            center.x + dx,
                            center.y + dy,
                            center.z + dz};

                        auto it = grid.find(neighbor);
                        if (it == grid.end())
                            continue;

                        const auto &cell_points = it->second;
                        for (int j : cell_points)
                        {
                            if (j <= i)
                                continue;

                            double dx_val = ix - coords[j * 3 + 0];
                            double dy_val = iy - coords[j * 3 + 1];
                            double dz_val = iz - coords[j * 3 + 2];

                            double dist_sq = dx_val * dx_val + dy_val * dy_val + dz_val * dz_val;

                            if (dist_sq <= thresh_sq)
                            {
                                local_contacts.push_back({i, j});
                            }
                        }
                    }
                }
            }
        }

#pragma omp critical
        {
            all_contacts.insert(all_contacts.end(), local_contacts.begin(), local_contacts.end());
        }
    }

    all_contacts.shrink_to_fit();
    return all_contacts;
}