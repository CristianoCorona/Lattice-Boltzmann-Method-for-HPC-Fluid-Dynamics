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
    int indices[d] = {coords...};
    int index = 0;
    int stride = 1;
    for(int i = 0; i < d; ++i) {
        index += indices[i] * stride;
        stride *= sizes[i];
    }
    return index;
}

// Initialize equilibrium with u=0 and rho=1, called only in the beginning of the simulation
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
template class Lattice<D3Q19<float>, float>;
template class Lattice<D3Q19<double>, double>;