#pragma once

#include <random>

// Simple global RNG for the project.
// For now: seeded from std::random_device.
// For reproducible tests, you can add a set_rng_seed() function later.

inline std::mt19937& global_rng() {
    static thread_local std::mt19937 rng(std::random_device{}());
    return rng;
}

// Uniform(0, 1)
inline double uniform01() {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(global_rng());
}

// Normal(mean, sd^2)
inline double rnorm(double mean, double sd) {
    std::normal_distribution<double> dist(mean, sd);
    return dist(global_rng());
}

// OPTIONAL: if you want reproducible randomness in tests, uncomment this:
//
// inline void set_rng_seed(unsigned int seed) {
//     global_rng().seed(seed);
// }
