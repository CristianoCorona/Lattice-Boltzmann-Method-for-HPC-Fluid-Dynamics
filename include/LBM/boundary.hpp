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
        
        WallsBoundary(std::array<int, Descriptor::d> dim, float_type wall_speed, Direction<Descriptor::d> moving_wall) : 
            speed(wall_speed),
        {
            if(moving_wall != Direction<d>::NODIR) {
                is_moving_wall = true;
                index_moving_wall = static_cast<unsigned char>(moving_wall - 1);
            } else {
                is_moving_wall = false;
            }
            for(int i = 0; i < d; ++i) {
                walls[2*i] = 0;
                walls[2*i + 1] = dim[i] - 1;
            }
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
