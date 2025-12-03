// no Node class (we should not iterate over objects)
// merge collision and propagation: stream_collide_save (page 565)
// (((each cell pulls from its neighbors (no push) to avoid race conditions (page 647))))
// ghost/halo cells => extra rows/columns
// (((pipeline idea: send boundary data with MPI -> compute inner cells -> receive boundary data with MPI -> compute boundary cells)))

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

template<std::floating_point float_type>
struct D2Q9 : public LatticeDescriptor<2, 9> {

    static constexpr int c[9][2] = {
        { 0, 0},
        { 1, 0}, { 0, 1}, {-1, 0}, { 0,-1},
        { 1, 1}, {-1, 1}, {-1,-1}, { 1,-1}
    };
    static constexpr float_type w[9] = {
        4.0 / 9.0,
        1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
    };
    static constexpr float_type cs2 = 1.0 / 3.0;

};

template <isDescriptor Descriptor, std::floating_point float_type = double> //customizable precision, default double
class Lattice {

public:

    static constexpr int d = Descriptor::d;
    static constexpr int q = Descriptor::q;
    int nx; 
    int ny; 
    const int total_cells = nx * ny;

    Lattice(int nx, int ny, float_type lid_velocity) : nx(nx), ny(ny), u(d), f_current(q), f_next(q), u_lid(lid_velocity){
        //TODO: ghost cells handling -> +2 padding?
        for(int i = 0; i < q; ++i) {
            f_current[i].resize(total_cells, 0.0); 
            f_next[i].resize(total_cells, 0.0);
        }
        for(int k = 0; k < d; ++k) {
            u[k].resize(total_cells, 0.0);
        }
        rho.resize(total_cells, 1.0); // Init rho a 1
        mask.resize(total_cells, FLUID);
        
        setup_CellType();
    }

    std::array<std::vector<float_type>, q> f_current;
    std::array<std::vector<float_type>, q> f_next;
    std::vector<float_type> rho;
    std::array<std::vector<float_type>, d> u;
    float_type u_lid;

    enum CellType : uint8_t {FLUID=0, WALL=1, LID=2};
    std::vector<CellType> mask;

    void swap_buffers();

    // functions for boundaries-> NOT needed with the mask approach
    //void boundary_values(const std::vector<int> &coords, float_type &rho_b);

    // functions for indices
    template <typename... Ints>
    [[nodiscard]] inline int idx(Ints...coords); //variadic template for indexing (general for 2D/3D)

    //Initialize equilibrium with u=0 and rho=1
    void initialize_equilibrium();

private:
    // function to setup the mask; called once in the constructor
    void setup_CellType();
    
};

template <isDescriptor Descriptor>
class Solver {

private:

    double tau;
    Lattice<Descriptor> &lattice;

    inline double feq(int i, double rho, double ux, double uy); // helper function to compute equilibrium

public:

    Solver(double kinematic_viscosity, Lattice<Descriptor> &lattice) : lattice(lattice); // constructor

    void compute_moments(Lattice<Descriptor>& grid);
    void stream_collide_inner(Lattice<Descriptor>& grid);
    // (((between these there will be MPI communication)))
    void stream_collide_boundaries(Lattice<Descriptor>& grid);
    void apply_bcs(Lattice<Descriptor>& grid, double u_lid);

    // function for the output

};

#endif