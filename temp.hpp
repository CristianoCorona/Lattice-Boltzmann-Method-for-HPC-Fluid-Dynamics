// no Node class (we should not iterate over objects)
// merge collision and propagation: stream_collide_save (page 565)
// (((each cell pulls from its neighbors (no push) to avoid race conditions (page 647))))
// ghost/halo cells => extra rows/columns
// (((pipeline idea: compute tile-boundary cells, send boundary data with MPI to adjacent tiles
//         -> compute inner cells while receiving boundary data with MPI -> compute boundary cells)))

// to store (Lattice):
// - f_current and f_next (ping-pong to avoid race conditions while computing collision and propagation)
// - u 
// - rho

// not to store (Solver):
// - feq (register variable)
// - fi* (collision and propagation are merged)

// ?:
// kinematic viscosity
// viscous stress tensor

#include <concepts>
#include <vector>
#include <array>

#ifndef LBM_HPP
#define LBM_HPP

template <int D, int Q>
struct LatticeDescriptor {

    static constexpr int d = D;
    static constexpr int q = Q;

};


template <typename T>
concept isDescriptor = requires{

    T::d;
    T::q;

} && std::derived_from<T, LatticeDescriptor<T::d, T::q>>;


struct D2Q9 : public LatticeDescriptor<2, 9> {

    static constexpr int c[9][2] = {
        { 0, 0},
        { 1, 0}, { 0, 1}, {-1, 0}, { 0,-1},
        { 1, 1}, {-1, 1}, {-1,-1}, { 1,-1}
    };
    static constexpr double w[9] = {
        4.0 / 9.0,
        1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
    };
    static constexpr double cs2 = 1.0 / 3.0;

};


template <isDescriptor Descriptor>
class Lattice {

public:

    Lattice(int nx, int ny) : u(Descriptor.d), f_current(Descriptor.q), f_next(Descriptor.q), u_w(Descriptor.d * 2, std::vector<double>(Descriptor.d)); // constructor

    std::array<std::vector<double>, Descriptor::q> f_current;
    std::array<std::vector<double>, Descriptor::q> f_next;
    std::array<std::vector<double>, Descriptor::d> u;
    std::vector<double> rho;
    std::vector<double> sigma;
    std::array<std::array<std::vector<double>, Descriptor::d>, 2> u_w;

    // Do we need this if we have a single function for the solve?
    void swap_buffers();

    // functions for boundaries
    void boundary_values(const std::vector<int> &coords, double &rho_b, std::vector<double> &u_w);

    // from matrix indeces return vector indicex
    inline int idx(const std::vector<int> &coords);
};


template <isDescriptor Descriptor>
class Solver {

private:

    double tau;
    Lattice<Descriptor> &lattice;

    // helper function to compute equilibrium
    inline double feq(int i, double rho, double ux, double uy); 

public:

    // constructor
    Solver(double kinematic_viscosity, Lattice<Descriptor> &lattice) : lattice(lattice); 

    void compute_moments(Lattice<Descriptor>& grid);
    void stream_collide_inner(Lattice<Descriptor>& grid);
    // (((between these there could be MPI communication)))
    void stream_collide_boundaries(Lattice<Descriptor>& grid);
    void apply_bcs(Lattice<Descriptor>& grid, double u_lid);

    // function for the output, maybe change second parameter, 
    // use flag to set what output you want
    void write_output(Lattice<Descriptor>& grid, char* out_buffer, char flag);

};

#endif