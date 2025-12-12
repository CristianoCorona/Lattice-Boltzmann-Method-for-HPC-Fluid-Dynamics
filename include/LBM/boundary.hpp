#include <vector>
#include <array>
#include "descriptor.hpp"

#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

template<isDescriptor Descriptor, std::floating_point float_type>
class WallsBoundary {
    public: 

        using DirEnum = typename Direction<Descriptor::d>::Value;

        static constexpr int d = Descriptor::d;
        static constexpr int q = Descriptor::q;
        
        WallsBoundary(std::array<int, Descriptor::d> dim, std::array<float_type, d> wall_speed, DirEnum moving_wall) : 
            wall_speed(wall_speed)
        {
            if(!(moving_wall == Direction<d>::NODIR)) {
                is_moving_wall = true;
                index_moving_wall = static_cast<unsigned char>(moving_wall) - 1;
            } else {
                is_moving_wall = false;
            }
            for(int i = 0; i < d; ++i) {
                walls[2*i] = 0;
                walls[2*i + 1] = dim[i] - 1;
            }
        }
        

        std::array<float_type, Descriptor::d> get_speed_of_wall(DirEnum wall);

        DirEnum is_at_bound(int cell_index);

        bool will_get_bounced_back(DirEnum wall, int direction);

        bool isMovingWall(DirEnum wall);

    protected:
    
        std::array<int, 2 * (d) > walls;

        const std::array<float_type, d> wall_speed;

        bool is_moving_wall;

        unsigned char index_moving_wall;

};

#endif // BOUNDARY_HPP
