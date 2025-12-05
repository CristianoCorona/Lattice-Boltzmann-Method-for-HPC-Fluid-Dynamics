// Lattice data structures for LBM simulation
// - f_current and f_next (ping-pong to avoid race conditions while computing collision and propagation)
// - u (velocity field)
// - rho (density field)

#include <concepts>
#include <vector>
#include <array>
#include <string>

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


template<std::floating_point float_type>
struct D3Q19 : public LatticeDescriptor<3, 19> {

    static constexpr std::array<std::array<int, 3>, 19> c = {{
        {{ 0, 0, 0}},
        {{ 1, 0, 0}}, {{ -1, 0, 0}}, {{ 0, 1, 0}}, {{ 0, -1, 0}}, {{ 0, 0, 1}}, {{ 0, 0, -1}},
        {{ 1, 1, 0}}, {{ 1, -1, 0}}, {{ -1, 1, 0}}, {{ -1, -1, 0}},
        {{ 1, 0, 1}}, {{ 1, 0, -1}}, {{ -1, 0, 1}}, {{ -1, 0, -1}},
        {{ 0, 1, 1}}, {{ 0, 1, -1}}, {{ 0, -1, 1}}, {{ 0, -1, -1}}
    }};
    static constexpr std::array<float_type, 19> w = {
        1.0 / 3.0,
        1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
    };
    static constexpr std::array<int, 19> opposite = {0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15};
};


template <isDescriptor Descriptor, std::floating_point float_type = double> //customizable precision, default double
class Lattice {

public:

    static constexpr int d = Descriptor::d;
    static constexpr int q = Descriptor::q;

    // sizes in each dimention
    std::array<int, d> sizes;
    int total_cells;

    // pre-computed strides for indexing
    std::array<int, d> strides;

    Lattice(std::array<int, d> dimensions, float_type lid_velocity, float_type dx = 1.0, std::string output_file_) : sizes(dimensions), f_current(q), f_next(q), u_lid(lid_velocity), dx(dx), output_file(output_file_){

        total_cells = 1;
        for(int s : dimensions) total_cells *= s;

        for(int i = 0; i < q; ++i) {
            f_current[i].resize(total_cells, 0.0); 
            f_next[i].resize(total_cells, 0.0);
        }

        u.resize(total_cells);
        for(int cell = 0; cell < total_cells; ++cell) {
            for(int j = 0; j < d; ++j) {
                u[cell][j] = 0.0;
            }
        }
        rho.resize(total_cells, 1.0);

        // strides computation
        int current_stride = 1;
        for (int i = 0; i < d; ++i) {
            strides[i] = current_stride;
            current_stride *= sizes[i];
        }
    }

    std::array<std::vector<float_type>, q> f_current;
    std::array<std::vector<float_type>, q> f_next;
    std::vector<float_type> rho;
    float_type u_lid, nu, delta_t, dx;
    const std::string output_file;

    // u has 2 or 3 components depending on dim
    template<int dim>
    struct u_struct {

        std::array<float_type, dim> data;
        float_type& operator[](int i) { return data[i]; }
        const float_type& operator[](int i) const { return data[i]; }
    };
    std::vector<u_struct<d>> u;

    // To swap f_current and f_next
    void swap_buffers();

    // functions for boundaries
    //void boundary_values(const std::vector<int> &coords, float_type &rho_b);

    // functions for indices
    //variadic template for indexing (general for 2D/3D)
    template <typename... Ints>
    [[nodiscard]] inline int idx(Ints...coords) const {
        static_assert(sizeof...(coords) == d, "Number of coordinates must match lattice dimension");
        // unroll stack array
        const int c[d] = {static_cast<int>(coords)...};
        // y * nx + x  (2D)
        if constexpr (d == 2) {
            return c[1] * strides[1] + c[0];
        } else if constexpr (d == 3) {
            // z * (nx * ny) + y * nx + x  (3D)
            return c[2] * strides[2] + c[1] * strides[1] + c[0];
        }
    } 

    //Initialize equilibrium with u=0 and rho=1
    void initialize_equilibrium();
    
    // Write a scalar field vector to VTK file for visualization
    void write_vtk(const std::vector<float_type>& data, const std::string& field_name);
    
};

#endif
