#include "LBM/solver.hpp"

/*
 *  stream_collide implementation, i is the direction in which it computes the step.
 */
template <isDescriptor Descriptor, std::floating_point float_type>
void Solver<Descriptor, float_type>::stream_collide(
            const int i,
            const float_type inv_tau_star,
            std::array<std::vector<float_type>, Descriptor::q> &f_next) {

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
        if (lattice.is_at_bound(index, i, rho_w, u_w)) { // adapt it according to Lattice
            int i_opp = Descriptor::opposite[i];
            float_type f_bounced = f_star - 2.0 * w_i * rho_w * scalar_prod(c_i, u_w) * inv_cs2;
            f_next[i_opp][index] = f_bounced;
        } else {
            int next_index = lattice.get_next_index(index, i);
            f_next[i][next_index] = f_star;
        }
    }
}

/*
 *  Function to compute rho_next and u_next.
 */
template <isDescriptor Descriptor, std::floating_point float_type>
void Solver<Descriptor, float_type>::update_moments(
            const std::array<std::vector<float_type>, Descriptor::q> &f_next,
            std::vector<float_type> &rho_next,
            std::array<std::vector<float_type>, Descriptor::d> &u_next) {
    
    std::array<const float_type*, Descriptor::q> f_ptrs;
    for (int i = 0; i < Descriptor::q; ++i) {
        f_ptrs[i] = f_next[i].data();
    }

    #pragma omp parallel for schedule(static)
    for (int index = 0; index < lattice.total_cells; ++index) {
        float_type local_rho = 0.0;
        std::array<float_type, Descriptor::d> local_u{};
        
        for (int i = 0; i < Descriptor::q; ++i) {
            float_type val = f_ptrs[i][index];
            local_rho += val;
            for (int d = 0; d < Descriptor::d; ++d) {
                local_u[d] += static_cast<float_type>(Descriptor::c[i][d]) * val;
            }
        }

        rho_next[index] = local_rho;
        if (local_rho > 1e-9) {
            float_type inv_rho = 1.0 / local_rho;
            for (int d = 0; d < Descriptor::d; ++d) {
                u_next[d][index] = local_u[d] * inv_rho;
            }
        } else {
            for (int d = 0; d < Descriptor::d; ++d) {
                u_next[d][index] = 0.0;
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
        const float_type delta_t,
        int output_interval) {
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

    lattice.initialize_equilibrium();

    lattice.write_vtk(0);
    for (unsigned long n_iter = 0; n_iter < n_iterations; ++n_iter) {
        /*
         *  This loop iterates over each possible direction, defined by the LatticeDescriptor,
         *  computing f_next by calling stream_collide.
         */
        for (int i = 0; i < Descriptor::q; ++i) {
            stream_collide(i, inv_tau_star, f_next);
        }

        update_moments(f_next, rho_next, u_next);

        lattice.swap_buffers(f_next, rho_next, u_next);
        if ((n_iter + 1) % output_interval == 0) {
            lattice.write_vtk(n_iter + 1);
        }
    }
    lattice.write_vtk(n_iterations);
}

template class Solver<D2Q9<float>, float>;
template class Solver<D2Q9<double>, double>;
template class Solver<D3Q19<float>, float>;
template class Solver<D3Q19<double>, double>;
