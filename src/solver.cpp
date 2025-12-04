#include "solver.hpp"

template <isDescriptor Descriptor, std::floating_point float_type = double>
void Solver::stream_collide(const int i, const float_type tau, const float_type delta_t) {
    constexpr auto scalar_prod = [](int *c_i, float_type *u_w) {
        float_type res = 0.0;
        for (int i = 0; i < Descriptor.d; ++i) {
            res += c_i[i] * u_w[i];
        }

        return res;
    }

    const int c_i[Descriptor.d] =  Descriptor.c[i];
    const float_type w_i = Descriptor.w[i];
    const std::vector<float_type> f_current = lattice.f_current[i];

    for (int y = 0; y < lattice.ny; ++y) {
        for (int x = 0; x < lattice.nx; ++x) {
            int index = y * nx + x; // lattice.idx(x, y)
            std::array<float_type, Descriptor.d> u = lattice.u[index];
            float_type rho = lattice.rho[index];
            float_type f_i = f_current[index];

            float_type f_eq = compute_eq(c_i, w_i, rho, u);

            /*
             * Here we may compute sigma (stress tensor)
             */

            float_type f_star = compute_star(f_i, f_eq, tau, delta_t);

            float_type rho_w = 0.0;
            float_type[Descriptor.d] u_w = 0.0;

            if (lattice.isAtBound({x, y}, i, rho_w, u_w)) {
                int i_opp = Descriptor.i_opp[i];
                lattice.f_next[i_opp][index] = f_star - 2.0 * w_i * rho_w * scalar_prod(c_i, u_w) * inv_cs2;
            } else {
                int next_index = (y + c_i[1]) * nx + (x + c_i[0]);
                lattice.f_next[i][next_index] = f_star;
            }
            
        }
    }
}


template <isDescriptor Descriptor, std::floating_point float_type = double>
void Solver::solve(const unsigned long n_iterations, const float_type delta_t) {
    const float_type tau = ((lattice.nu * delta_t * inv_cs2) / (lattice.delta_x * lattice.delta_x) + 0.5) * delta_t;

    // write output
    for (unsigned long n_iter = 0; n_iter < n_iterations; ++n_iter) {
        // This loop iterates over each possible direction, defined by the LatticeDescriptor, computing f_next
        for (int i = 0; i < Descriptor.q; ++i) {
            stream_collide(i, tau, delta_t);
        }
        compute_moments();
        // write output
    }
}
