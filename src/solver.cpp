#include "solver.hpp"

/*
 *  stream_collide implementation, i is the direction in which it computes the step.
 */
template <isDescriptor Descriptor, std::floating_point float_type>
void Solver<Descriptor, float_type>::stream_collide(const int i, const float_type inv_tau_star) {
    /*
     *  Lambda function to compute the scalar-vector product.
     */
    constexpr auto prod = [](auto &a, auto b) {
        std::array<float_type, Descriptor::d> res{};
        for (int i = 0; i < Descriptor::d; ++i) {
            res[i] += static_cast<float_type>(a[i]) * static_cast<float_type>(b);
        }

        return res;
    };

    /*
     *  Lambda function to compute the vector sum.
     */
    constexpr auto vector_sum = [](auto &a, auto &b) {
        std::array<float_type, Descriptor::d> res{};
        for (int i = 0; i < Descriptor::d; ++i) {
            res[i] += static_cast<float_type>(a[i]) + static_cast<float_type>(b[i]);
        }

        return res;
    };

    const int c_i[Descriptor::d] =  Descriptor::c[i];
    const float_type w_i = Descriptor::w[i];
    const std::vector<float_type> &f_current = lattice.f[i];

    /*
     *  Iterates over the entire lattice computing f_eq, f_star and f_i for the
     *  given direction i and summing the contribution to the next density and
     *  velocity of the cells.
     */
    #pragma omp parallel for schedule(static) collapse(2)
    for (int y = 0; y < lattice.sizes[1]; ++y) {
        for (int x = 0; x < lattice.sizes[0]; ++x) {
            int index = lattice.idx(x, y);
            std::array<float_type, Descriptor::d> u = lattice.u[index];
            float_type rho = lattice.rho[index];
            float_type f_i = f_current[index];

            float_type f_eq = compute_eq(c_i, w_i, rho, u);

            // Here we may compute sigma (stress tensor)

            float_type f_star = compute_star(f_i, f_eq, inv_tau_star);

            float_type rho_w = 0.0;
            std::array<float_type, Descriptor::d> u_w{};

            /*
             *  Checks if the current cell is a boundary one, if so applies
             *  the bounce back method with the boundary parameters rho_w and u_w,
             *  otherwise compute the propagation with the push method.
             */
            if (lattice.isAtBound({x, y}, i, rho_w, u_w)) { // adapt it according to Lattice
                int i_opp = Descriptor::i_opp[i];
                f_i = f_star - 2.0 * w_i * rho_w * scalar_prod(c_i, u_w) * inv_cs2;
                lattice.f_next[i_opp][index] = f_i;
            } else {
                int next_index = (y + c_i[1]) * nx + (x + c_i[0]);
                f_i = f_star;
                lattice.f_next[i][next_index] = f_i;
            }

            /*
             *  Compute the contribute that f_i gives to rho_next and u_next.
             */
            if (i == 0) {
                rho_next = f_i;
                u_next = prod(c_i, f_i);
            } else {
                rho_next += f_i;
                u_next = vector_sum(u_next, prod(c_i, f_i));
                if (i == Descriptor::q - 1) {
                    u_next = prod(u_next, 1.0 / rho_next);
                }
            }
        }
    }
}

/*
 *  solve implementation, n_iterations is the length of the simulation
 *  and delta_t is the time interval length of each step.
 */
template <isDescriptor Descriptor, std::floating_point float_type>
void Solver<Descriptor, float_type>::solve(
        const unsigned long n_iterations,
        const float_type delta_t) {
    /*
     *  Since tau and delta_t are always used as a ratio between them
     *  and it is constant over the entire simulation, we compute it
     *  with respect of tau's relation with problem's parameters.
     */
    const float_type inv_tau_star = 1.0 / (
                (lattice.nu * delta_t * inv_cs2) / 
                (lattice.dx * lattice.dx) + 
                0.5
            );

    std::array<std::vector<float_type>, Descriptor::q> f_next;
    std::array<std::vector<float_type>, Descriptor::d> u_next;
    std::vector<float_type> rho_next;

    // write output
    for (unsigned long n_iter = 0; n_iter < n_iterations; ++n_iter) {
        /*
         *  This loop iterates over each possible direction, defined by the LatticeDescriptor,
         *  computing f_next by calling stream_collide.
         */
        for (int i = 0; i < Descriptor::q; ++i) {
            stream_collide(i, inv_tau_star);
        }
        // write output
        // swap
    }
}
