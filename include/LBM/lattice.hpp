// Lattice data structures for LBM simulation
// - u (velocity field)
// - rho (density field)
#include <concepts>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include "descriptor.hpp"
#include "boundary.hpp"

#ifndef LATTICE_HPP
#define LATTICE_HPP

// Lattice for LBM simulation
// Handles:
// - Distribution functions (f_current and f_next)
// - Data storage for velocity field (u) and density (rho) 
// - Indexing and boundary conditions
// - VTK output for visualization
template <isDescriptor Descriptor, std::floating_point float_type = double> //customizable precision, default double
class Lattice {

public:

    using DirEnum = typename Direction<Descriptor::d>::Value;

    /*
     *   Lattice dimensions and properties
     */ 
    static constexpr int d = Descriptor::d;
    static constexpr int q = Descriptor::q;
    std::array<int, d> sizes;
    int total_cells;

    /*
     *   Pre-computed strides for indexing
     */
    std::array<int, d> strides;


    /*
     *   Neighbor offsets for each direction c_i -> used in get_next_index
     */
    std::array<int,q> neighbor_offsets;

    Lattice(std::array<int, d> dimensions, 
            float_type lid_velocity, 
            float_type nu_,
            DirEnum moving_wall_,
            std::string output_file_ = "output.vtk"
            ) 
            : sizes(dimensions), 
            u_lid(lid_velocity), 
            nu(nu_), dx(dx), 
            rho_init(rho_init_),
            rho_wall(rho_wall_),
            rho_lid(rho_lid_),
            output_file(output_file_),
            float_type dx = 1.0, 
            float_type rho_init_=1.0,
            float_type rho_wall_=1.0,
            float_type rho_lid_=1.0,
            {
        // total cells computation
        total_cells = 1;
        for(int s : dimensions) total_cells *= s;

        // resize vectors
        for(int i = 0; i < q; ++i) {
            f[i].resize(total_cells, 0.0); 
        }

        for(int i = 0; i < d; ++i) {
            u[i].resize(total_cells, 0.0);
        }
        rho.resize(total_cells, 0.0);
        rho.assign(total_cells, rho_init);

        // strides computation
        int current_stride = 1;
        for (int i = 0; i < d; ++i) {
            strides[i] = current_stride;
            current_stride *= sizes[i];
        } //strides -> (x -> 1, y -> nx, z -> nx*ny)

        // neighbor offsets computation
        for (int i = 0; i < q; ++i) {
            int offset = 0;
            if constexpr (d == 2) {
                // D2Q9 offsets -> (c_0 -> 0, c_1 -> 1, c_2 -> nx, 
                             //c_3 -> -1, c_4 -> -nx, c_5 -> nx+1, 
                             //c_6 -> nx-1, c_7 -> -nx-1, c_8 -> -nx+1)
                offset = Descriptor::c[i][0] * strides[0] + Descriptor::c[i][1] * strides[1];
            } else if constexpr (d == 3) {
                // D3Q19 offsets -> (c_0 -> 0, c_1 -> 1, c_2 -> -1, c_3 -> nx, c_4 -> -nx, 
                             //c_5 -> nx*ny, c_6 -> -nx*ny, c_7 -> nx+1, c_8 -> nx-1,
                             //c_9 -> -nx+1, c_10 -> -nx-1, c_11 -> nx*ny+1,
                             //c_12 -> nx*ny-1, c_13 -> -nx*ny+1, c_14 -> -nx*ny-1,
                             //c_15 -> ny+1, c_16 -> ny-1, c_17 -> -ny+1, c_18 -> -ny-1)
                offset = Descriptor::c[i][0] * strides[0] + 
                         Descriptor::c[i][1] * strides[1] + 
                         Descriptor::c[i][2] * strides[2];
            }
            neighbor_offsets[i] = offset;
        }

        // initialize walls boundary
        walls_boundary = WallsBoundary<Descriptor, float_type>(sizes, u_lid, moving_wall_);
    }

    /*
     *   Distribution functions (f_current and f_next)
     *   Density (rho) and velocity (u) fields
     *   Physical parameters and output file
     */
    std::array<std::vector<float_type>, q> f;
    float_type rho_init, rho_wall, rho_lid, u_lid, nu, dx;
    std::vector<float_type> rho;
    std::array<std::vector<float_type>, d> u;
    const std::string output_file;

    /*
     * Boundary handling
     */
    WallsBoundary <Descriptor, float_type> walls_boundary;


    /*
     *   Swap the distribution function buffers (f_current and f_next)
     */
    void swap_buffers(std::array<std::vector<float_type>, q> f,
                      std::vector<float_type> rho, 
                      std::array<std::vector<float_type>, d> u);

    /*
     *   This function checks if the given node (passed through linearized index for generality)
     *   is on a boundary in the given direction.
     *   If so, it sets the boundary velocity value (u_b) and
     *   the boundary density value (rho_b) accordingly.
     *   Returns true if on boundary, false otherwise
     */
    bool isAtBound(int index, int direction, float_type &rho_b, float_type &u_b);

    /*
     *   Function to get the next index in the lattice given the current index and direction c_i
     */
    [[nodiscard]] inline
    int get_next_index(int current_index, int direction) const {
        return current_index + neighbor_offsets[direction];
    }

    /*
     *   Optimized indexing function
     *   Variadic template for indexing (general for 2D/3D)
     */
    template <typename... Ints>
    [[nodiscard]] inline 
    int idx(Ints...coords) const {
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

    /*
     *   Initialize equilibrium with u=0 and rho=1
     */
    void initialize_equilibrium();
    
    /*
     *   Write a scalar field vector to VTK file for visualization
     */
    void write_vtk(const std::vector<float_type>& data, const std::string& field_name, int iter);
    
    /*
     *   Write a vector field to VTK file for visualization (overload)
     */
    void write_vtk(const std::array<std::vector<float_type>, Descriptor::d>& data, const std::string& field_name, int iter);
    
};

#endif
