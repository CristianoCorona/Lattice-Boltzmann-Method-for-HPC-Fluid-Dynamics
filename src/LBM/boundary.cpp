#include "LBM/boundary.hpp"

template<isDescriptor Descriptor, std::floating_point float_type>
float_type WallsBoundary<Descriptor, float_type>::get_speed_of_wall(Direction<Descriptor::d> wall) {
    if(is_moving_wall && wall == static_cast<Direction<d>>(index_moving_wall + 1)) {
        return wall_speed;
    }
    return 0.0;
}


template<isDescriptor Descriptor, std::floating_point float_type>
Direction<Descriptor::d> WallsBoundary<Descriptor, float_type>::is_at_bound(int cell_index) {
    int coord[d];
    for(int i = 0; i < d; ++i) {
        coord[i] = cell_index % (walls[2*i + 1] + 1);
        cell_index /= (walls[2*i + 1] + 1);
    }
    for(int i = 0; i < walls.size(); i++){
        if(coord[i/2] == walls[i] || coord[i/2] == walls[i]){
            return static_cast<Direction<d>>(i);
        }
    }
    return Direction<d>::NODIR;
}


template<isDescriptor Descriptor, std::floating_point float_type>
bool WallsBoundary<Descriptor, float_type>::will_get_bounced_back(Direction<Descriptor::d> wall, int direction)
{
    if(wall != Direction<d>::NODIR) {
        if(Descriptor::c[direction][static_cast<int>(wall) / 2] == (static_cast<int>(wall) % 2 == 1) * 2 - 1) {
            return true;
        }
        else {
            return false;
        }
    }
    return false;
}

template class WallsBoundary<D2Q9<float>, float>;
template class WallsBoundary<D2Q9<double>, double>;
template class WallsBoundary<D3Q19<float>, float>;
template class WallsBoundary<D3Q19<double>, double>;
