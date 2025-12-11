#include <vector>
#include <array>
#include "descriptor.hpp"

#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

template<isDescriptor Descriptor, std::floating_point float_type>
class WallsBoundary {
    public: 

        static constexpr int d = Descriptor::d;
        static constexpr int q = Descriptor::q;
        
        WallsBoundary(int Nx, int Ny, int Nz, float_type wall_speed, Direction<Descriptor::d> moving_wall) : 
            speed(wall_speed),
            requires (d == 3)
        {
            if(moving_wall != Direction<d>::NODIR) {
                is_moving_wall = true;
                index_moving_wall = static_cast<unsigned char>(moving_wall - 1);
            } else {
                is_moving_wall = false;
            }
            walls = {
                0, Nx - 1,    // LEFT, RIGHT
                0, Ny - 1,    // BOTTOM, TOP
                0, Nz - 1     // BACK, FRONT
            };
        }

        WallsBoundary(int Nx, int Ny, float_type wall_speed, Direction<Descriptor::d> moving_wall) : 
            speed(wall_speed),
            requires (d == 2)
        {
            if(moving_wall != Direction<d>::NODIR) {
                is_moving_wall = true;
                index_moving_wall = static_cast<unsigned char>(moving_wall - 1);
            } else {
                is_moving_wall = false;
            }
            walls = {
                0, Nx - 1,    // LEFT, RIGHT
                0, Ny - 1     // BOTTOM, TOP
            };
        }

        float_type get_speed_of_wall(Direction<Descriptor::d> wall);

        Direction<Descriptor::d> is_at_bound(int cell_index);

        bool will_get_bounced_back(Direction<Descriptor::d> wall, int direction);


    protected:
    
        static const std::array<int, 2 * (d) > walls;

        const float_type wall_speed;

        const bool is_moving_wall;

        const unsigned char index_moving_wall;

};

#endif // BOUNDARY_HPP
