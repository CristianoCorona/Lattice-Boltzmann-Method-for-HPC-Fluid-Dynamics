#include "LBM/boundary.hpp"

template<isDescriptor Descriptor, std::floating_point float_type>
std::array<float_type, Descriptor::d> WallsBoundary<Descriptor, float_type>::get_speed_of_wall(DirEnum wall) {
    if(is_moving_wall && wall == static_cast<DirEnum>(index_moving_wall + 1)) {
        return wall_speed;
    }
    return std::array<float_type, Descriptor::d>{};
}


template<isDescriptor Descriptor, std::floating_point float_type>
typename WallsBoundary<Descriptor, float_type>::DirEnum WallsBoundary<Descriptor, float_type>::is_at_bound(int cell_index, int direction) {
    int coord[d];
    for(int i = 0; i < d; ++i) {
        coord[i] = cell_index % (walls[2*i + 1] + 1);
        cell_index /= (walls[2*i + 1] + 1);
    }
    for(int i = 0; i < walls.size(); i++){
        if(coord[i/2] == walls[i]){
            DirEnum wall = static_cast<DirEnum>(i + 1);
            if (will_get_bounced_back(wall, direction)) {
                return wall;
            }
        }
    }
    return Direction<d>::NODIR;
}


template<isDescriptor Descriptor, std::floating_point float_type>
bool WallsBoundary<Descriptor, float_type>::will_get_bounced_back(DirEnum wall, int direction)
{
    if(wall != Direction<d>::NODIR) {
        unsigned char w = static_cast<unsigned char>(wall) - 1;
        if(Descriptor::c[direction][w / 2] == (w % 2 == 1) * 2 - 1) {
            return true;
        }
        else {
            return false;
        }
    }
    return false;
}



template<isDescriptor Descriptor, std::floating_point float_type>
bool WallsBoundary<Descriptor, float_type>::isMovingWall(DirEnum wall) {
    if(is_moving_wall && (wall == static_cast<DirEnum>(index_moving_wall + 1))) {
        return true;
    }
    return false;
}

template class WallsBoundary<D2Q9<float>, float>;
template class WallsBoundary<D2Q9<double>, double>;
template class WallsBoundary<D3Q19<float>, float>;
template class WallsBoundary<D3Q19<double>, double>;