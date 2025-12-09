#include <vector>
#include <array>
#include "descriptor.hpp"

#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

template<isDescriptor Descriptor, std::floating_point float_type>
class WallsBoundary {
    public: 
        
        WallsBoundary(int Nx, int Ny, int Nz, float_type wall_speed, Direction<Descriptor::d> moving_wall) : 
            speed(wall_speed),
            index_moving_wall(static_cast<unsigned char>(moving_wall))
            requires (Descriptor::d == 3)
        {
            walls = {
                0, Nx - 1,    // LEFT, RIGHT
                0, Ny - 1,    // BOTTOM, TOP
                0, Nz - 1     // BACK, FRONT
            };
        }

        WallsBoundary(int Nx, int Ny, float_type wall_speed, Direction<Descriptor::d> moving_wall) : 
            speed(wall_speed),
            index_moving_wall(static_cast<unsigned char>(moving_wall))
            requires (Descriptor::d == 2)
        {
            walls = {
                0, Nx - 1,    // LEFT, RIGHT
                0, Ny - 1     // BOTTOM, TOP
            };
        }

        std::array<int, 2 * (Descriptor::d) > walls;

        float_type wall_speed;

        unsigned char index_moving_wall;

        float_type get_speed_of_wall(int cell_index, int direction);

        bool is_at_bound(int cell_index, int direction);

};

#endif // BOUNDARY_HPP
