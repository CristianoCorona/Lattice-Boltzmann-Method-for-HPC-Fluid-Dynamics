#include "lattice.hpp"
#include <array>

#ifndef SOLVER_HPP
#define SOLVER_HPP

// Solver for LBM simulation
// Handles:
// - Stream and collide operations (merged to avoid race conditions)
// - Boundary conditions
// - Moment computation
template <isDescriptor Descriptor, std::floating_point float_type = double>
class Solver {

private:

    /*
     *  Represents the lattice where the problem is defined.
     */
    Lattice<Descriptor, float_type> &lattice;
    
    /*  Since cs^2 is constant and it's a problem property it's more efficient to store
     *  it as a constexpr inside the Solver to make the compiler able to insert its value
     *  at compile time inside the code, preventing the read of the same value from
     *  Lattice at each iteration.
    */
    static constexpr float_type inv_cs2 = 3.0; // 1.0 / (1.0 / 3.0)
    static constexpr float_type inv_2cs2 = 1.5; // 2.0 * 1.0 / (1.0 / 3.0)
    static constexpr float_type inv_2cs4 = 4.5; // 2.0 * 1.0 / (1.0 / 3.0)^2

    /*
     *  Lambda function used to compute the scalar product between two arrays,
     *  it returns the value as a float_type.
     */
    constexpr auto scalar_prod = [](auto &a, auto &b) {
        float_type res = 0.0;
        for (int i = 0; i < Descriptor::d; ++i) {
            res += static_cast<float_type>(a[i]) * static_cast<float_type>(b[i]);
        }

        return res;
    }

    /*
     *  Helper function to compute equilibrium in a single direction.
     */
    inline float_type compute_eq(
            const std::array<int,
            Descriptor::d> &c_i,
            const float_type w_i,
            const float_type rho,
            const std::array<float_type,
            Descriptor::d> &u) {
        return w_i * rho * (
                1.0 + inv_cs2 * scalar_prod(c_i, u) + 
                inv_2cs4 * scalar_prod(c_i, u) * scalar_prod(c_i, u) - 
                inv_2cs2 * scalar_prod(u, u)
                );
    }; 
    
    /*
     *  Helper function to compute collision in a single direction.
     */
    inline float_type compute_star(
            const float_type f_i,
            const float_type feq_i,
            const float_type inv_tau_star) {
        return f_i * (1 - inv_tau_star) + feq_i * inv_tau_star;
    };
    
    /*
     *  This function encapsulates the computation of equilibrium (via compute_eq),
     *  collision (via compute_star), propagation in both inner and boundary cells
     *  for every direction. The function also computes density and velocity of each
     *  cell of the lattice.
     */
    void stream_collide(const int i, const float_type inv_tau_star);
    
public:

    /*
     *  Constructor takes the lattice as argument.
     */
    Solver(Lattice<Descriptor, float_type> &lattice_) : lattice(lattice_); 

    /*
     *  Performs the simulation over a user-defined number of iterations (n_iterations)
     *  where each step reprents a user-defined time interval of lenght delta_t.
     */
    void solve(const unsigned long n_iterations, const float_type delta_t);
    
};

#endif
