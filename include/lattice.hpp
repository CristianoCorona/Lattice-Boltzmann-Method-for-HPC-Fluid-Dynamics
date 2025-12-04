// Lattice data structures for LBM simulation
// - f_current and f_next (ping-pong to avoid race conditions while computing collision and propagation)
// - u (velocity field)
// - rho (density field)

#include <concepts>
#include <vector>
#include <array>

#ifndef LATTICE_HPP
#define LATTICE_HPP

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

    static constexpr std::array<std::array<int, 2>, 9> c = {{
        {{ 0, 0}},
        {{ 1, 0}}, {{ 0, 1}}, {{-1, 0}}, {{ 0,-1}},
        {{ 1, 1}}, {{-1, 1}}, {{-1,-1}}, {{ 1,-1}}
    }};
    static constexpr std::array<float_type, 9> w = {
        4.0 / 9.0,
        1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
    };
    static constexpr std::array<int, 9> opposite = {0, 3, 4, 1, 2, 7, 8, 5, 6};

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

        for(int i = 0; i < q; ++i) {
            f_current[i].resize(total_cells, 0.0); 
            f_next[i].resize(total_cells, 0.0);
        }

        for(int k = 0; k < d; ++k) {
            u[k].resize(total_cells, 0.0);
        }
        rho.resize(total_cells, 1.0); // Init rho a 1
    }

    std::array<std::vector<float_type>, q> f_current;
    std::array<std::vector<float_type>, q> f_next;
    std::vector<float_type> rho;
    std::array<std::vector<float_type>, d> u;
    float_type u_lid;

    // To swap f_current and f_next
    void swap_buffers();

    // functions for boundaries
    // #################################################################################
    // we should instead compute the boundares first, to avoid the if condition on every cell at every iteration
    //void boundary_values(const std::vector<int> &coords, float_type &rho_b);

    // functions for indices
    template <typename... Ints>
    [[nodiscard]] inline int idx(Ints...coords); //variadic template for indexing (general for 2D/3D)

    //Initialize equilibrium with u=0 and rho=1
    void initialize_equilibrium();
    
};

#endif
