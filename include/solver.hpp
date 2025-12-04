// Solver for LBM simulation
// Handles:
// - Stream and collide operations (merged to avoid race conditions)
// - Boundary conditions
// - Moment computation
// - Output writing
//
// Notes:
// - feq (equilibrium) computed as register variable
// - fi* (collision and propagation are merged)
// - Each cell pulls from its neighbors (no push) to avoid race conditions

// assumption: Lattice has a function to which we pass coordinates, direction, rho_w, u_w and it gives us a boolean (to indicate if it's a 
// boundary or not) and, if true, also rho_w, u_w of the wall/lid through parameters
// so boundary_values is a bool: bool boundary_values(const std::vector<int> &coords, int i, float_type &rho_b, float_type &u_b);

#include "lattice.hpp"

#ifndef SOLVER_HPP
#define SOLVER_HPP

template <isDescriptor Descriptor, std::floating_point float_type = double>
class Solver {

private:

    Lattice<Descriptor, float_type> &lattice;
    
    // Since cs^2 is constant and it's a problem property it's more efficent to store it as a constexpr inside the Solver
    // to make the compiler able to inster its value at compile time inside the code preventing the read of the same value
    // from Lattice at each iteration
    constexpr float_type inv_cs2 = 3.0; // 1.0 / (1.0 / 3.0)
    constexpr float_type inv_2cs2 = 1.5; // 2.0 * 1.0 / (1.0 / 3.0)
    constexpr float_type inv_2cs4 = 4.5; // 2.0 * 1.0 / (1.0 / 3.0)^2

    constexpr auto scalar_prod = [](auto &a, auto &b) {
        float_type res = 0.0;
        for (int i = 0; i < Descriptor::d; ++i) {
            res += static_cast<float_type>(a[i]) * static_cast<float_type>(b[i]);
        }

        return res;
    }

    // helper function to compute equilibrium in a single direction
    inline float_type compute_eq(const std::array<int, Descriptor::d> &c_i, const float_type w_i, const float_type rho, const std::array<float_type, Descriptor::d> &u){
        return w_i * rho * (1.0 + inv_cs2 * scalar_prod(c_i, u) + inv_2cs4 * scalar_prod(c_i, u) * scalar_prod(c_i, u) - inv_2cs2 * scalar_prod(u, u));
    }; 
    
    // helper function to compute collision in a single direction
    inline float_type compute_star(const float_type f_i, const float_type feq_i, const float_type inv_tau_star){
        return f_i * (1 - inv_tau_star) + feq_i * inv_tau_star;
    };
    
    // function that computes f_next in a single direction
    void stream_collide(const int i, const float_type inv_tau_star);
    
public:

    // constructor
    Solver(Lattice<Descriptor, float_type> &lattice) : lattice(lattice); 

    // wraps the entire simulation over a user-defined number of iterations (n_iterations) where each step reprents a user-defined time interval
    // of lenght delta_t
    void solve(const unsigned long n_iterations, const float_type delta_t);
    
};

#endif
