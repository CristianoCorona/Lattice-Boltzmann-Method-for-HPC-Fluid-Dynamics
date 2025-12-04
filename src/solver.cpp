#include "solver.hpp"

template <isDescriptor Descriptor, std::floating_point float_type>
void Solver<Descriptor, float_type>::stream_collide(const int i, const float_type inv_tau_star) {

    constexpr auto prod = [](auto &a, auto b) {
        std::array<float_type, Descriptor::d> res();
        for (int i = 0; i < Descriptor::d; ++i) {
            res += static_cast<float_type>(a[i]) * static_cast<float_type>(b);
        }

        return res;
    }

    constexpr auto vector_sum = [](auto &a, auto &b) {
        std::array<float_type, Descriptor::d> res();
        for (int i = 0; i < Descriptor::d; ++i) {
            res[i] += static_cast<float_type>(a[i]) + static_cast<float_type>(b[i]);
        }

        return res;
    }

    const int c_i[Descriptor::d] =  Descriptor::c[i];
    const float_type w_i = Descriptor::w[i];
    const std::vector<float_type> f_current = lattice.f_current[i];



    for (int y = 0; y < lattice.ny; ++y) {
        for (int x = 0; x < lattice.nx; ++x) {
            int index = y * nx + x; // lattice.idx(x, y)
            std::array<float_type, Descriptor::d> u = lattice.u_current[index];
            float_type rho = lattice.rho_current[index];
            float_type f_i = f_current[index];

            float_type f_eq = compute_eq(c_i, w_i, rho, u);

            /*
             * Here we may compute sigma (stress tensor)
             */

            float_type f_star = compute_star(f_i, f_eq, tau, delta_t);

            float_type rho_w = 0.0;
            std::array<float_type, Descriptor::d> u_w();

            if (lattice.isAtBound({x, y}, i, rho_w, u_w)) { // adapt it according to Lattice
                int i_opp = Descriptor::i_opp[i];
                f_i = f_star - 2.0 * w_i * rho_w * scalar_prod(c_i, u_w) * inv_cs2;
                lattice.f_next[i_opp][index] = f_i;
            } else {
                int next_index = (y + c_i[1]) * nx + (x + c_i[0]);
                f_i = f_star;
                lattice.f_next[i][next_index] = f_i;
            }

            float_type &rho_next = lattice.rho_next[index];
            std::array<float_type, Descriptor::d> &u_next = lattice.u_next[index];
            if (i == 0){
                rho_next = f_i;
                u_next = prod(c_i, f_i);
            }
            else if (i == Descriptor::q - 1){
                rho_next += f_i;
                u_next = vector_sum(u_next, prod(c_i, f_i));
                u_next = prod(u_next, 1 / rho_next);

            }
            else {
                rho_next += f_i;
                u_next = vector_sum(u_next, prod(c_i, f_i));
            }
            
        }
    }
}

template <isDescriptor Descriptor, std::floating_point float_type>
void Solver<Descriptor, float_type>::solve(const unsigned long n_iterations, const float_type delta_t) {
    const float_type inv_tau_star = 1.0 / ((lattice.nu * delta_t * inv_cs2) / (lattice.delta_x * lattice.delta_x) + 0.5);

    // write output
    for (unsigned long n_iter = 0; n_iter < n_iterations; ++n_iter) {
        // This loop iterates over each possible direction, defined by the LatticeDescriptor, computing f_next
        for (int i = 0; i < Descriptor::q; ++i) {
            stream_collide(i, inv_tau_star);
        }
        // write output
        // swap
    }
}
