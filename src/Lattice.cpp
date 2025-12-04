#include "lattice.hpp"

// Swap f_current and f_next buffers (ping-pong)
template <isDescriptor Descriptor, std::floating_point float_type>
void Lattice<Descriptor, float_type>::swap_buffers() {
    std::swap(f_current, f_next);
}

// Variadic template for indexing (works for 2D/3D)
template <isDescriptor Descriptor, std::floating_point float_type>
template <typename... Ints>
int Lattice<Descriptor, float_type>::idx(Ints...coords) {
    static_assert(sizeof...(coords) == d, "Number of coordinates must match dimensionality");
    if constexpr (d == 2) {
        // For 2D: idx = i + j * nx
        int indices[] = {coords...};
        return indices[0] + indices[1] * nx;
    } else if constexpr (d == 3) {
        // For 3D: idx = i + j * nx + k * nx * ny
        int indices[] = {coords...};
        return indices[0] + indices[1] * nx + indices[2] * nx * ny;
    }
}

// Initialize equilibrium with u=0 and rho=1
template <isDescriptor Descriptor, std::floating_point float_type>
void Lattice<Descriptor, float_type>::initialize_equilibrium() {
    for (int cell = 0; cell < total_cells; ++cell) {
        for (int i = 0; i < q; ++i) {
            // f_eq_i = w_i * rho (when u = 0)
            f_current[i][cell] = Descriptor::w[i] * rho[cell];
            f_next[i][cell] = f_current[i][cell];
        }
    }
}

// The first type is the descriptor's weights and velocities
// THe second type is the lattice's dtributions and other datas
template class Lattice<D2Q9<float>, float>;
template class Lattice<D2Q9<double>, double>;