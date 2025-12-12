#include "LBM/solver.hpp"

/*
 *  stream_collide implementation, i is the direction in which it computes the step.
 */
template <isDescriptor Descriptor, std::floating_point float_type>
void Solver<Descriptor, float_type>::stream_collide(
            const int i,
            const float_type inv_tau_star,
            std::array<std::vector<float_type>, Descriptor::q> &f_next,
            std::array<std::vector<float_type>, Descriptor::d> &u_next,
            std::vector<float_type> &rho_next) {

    const std::array<int, Descriptor::d> c_i = Descriptor::c[i];
    const float_type w_i = Descriptor::w[i];
    const std::vector<float_type> &f_current = lattice.f[i];

    /*
     *  Iterates over the entire lattice computing f_eq, f_star and f_i for the
     *  given direction i and summing the contribution to the next density and
     *  velocity of the cells.
     */
    #pragma omp parallel for schedule(static)
    for(int index = 0; index < lattice.total_cells; ++index) {
        std::array<float_type, Descriptor::d> u;
        for (int d = 0; d < Descriptor::d; ++d) {
            u[d] = lattice.u[d][index];
        }
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
        if (lattice.isAtBound(index, i, rho_w, u_w)) { // adapt it according to Lattice
            int i_opp = Descriptor::i_opp[i];
            f_i = f_star - 2.0 * w_i * rho_w * scalar_prod(c_i, u_w) * inv_cs2;
            f_next[i_opp][index] = f_i;
        } else {
            int next_index = lattice.get_next_index(index, c_i);
            f_i = f_star;
            f_next[i][next_index] = f_i;
        }

        /*
         *  Compute the contribute that f_i gives to rho_next and u_next.
         */
        if (i == 0) {
            rho_next[index] = f_i;
            for (int d = 0; d < Descriptor::d; ++d) {
                u_next[d][index] = static_cast<float_type>(c_i[d]) * f_i;
            }
        } else {
            rho_next[index] += f_i;
            for (int d = 0; d < Descriptor::d; ++d) {
                u_next[d][index] += static_cast<float_type>(c_i[d]) * f_i;
            }
            if (i == Descriptor::q - 1) {
                float_type inv_rho = 1.0 / rho_next[index];
                for (int d = 0; d < Descriptor::d; ++d) {
                    u_next[d][index] *= inv_rho;
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

    for (int i = 0; i < Descriptor::q; ++i) {
        f_next[i].resize(lattice.total_cells, 0.0);
    }
    for (int i = 0; i < Descriptor::d; ++i) {
        u_next[i].resize(lattice.total_cells, 0.0);
    }
    rho_next.resize(lattice.total_cells, 0.0);

    lattice.write_vtk(lattice.rho, "rho", 0);
    lattice.write_vtk(lattice.u, "u", 0);
    for (unsigned long n_iter = 0; n_iter < n_iterations; ++n_iter) {
        /*
         *  This loop iterates over each possible direction, defined by the LatticeDescriptor,
         *  computing f_next by calling stream_collide.
         */
        for (int i = 0; i < Descriptor::q; ++i) {
            stream_collide(i, inv_tau_star, f_next, u_next, rho_next);
        }

        lattice.swap_buffers(f_next, rho_next, u_next);
        if ((n_iter + 1) % 20 == 0) {
            lattice.write_vtk(lattice.rho, "rho", n_iter + 1);
            lattice.write_vtk(lattice.u, "u", n_iter + 1);
        }
    }
    lattice.write_vtk(lattice.rho, "rho", n_iterations);
    lattice.write_vtk(lattice.u, "u", n_iterations);
}

template class Solver<D2Q9<float>, float>;
template class Solver<D2Q9<double>, double>;
