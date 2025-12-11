#include "lattice.hpp"
#include <utility>

// Swap f, rho and u
template <isDescriptor Descriptor, std::floating_point float_type>
void Lattice<Descriptor, float_type>::swap_buffers(std::array<std::vector<float_type>, q> f,
                                                    std::vector<float_type> rho, 
                                                    std::array<std::vector<float_type>, d> u) {
    std::swap(this->f, f);
    std::swap(this->rho, rho);
    std::swap(this->u, u);
}

template <isDescriptor Descriptor, std::floating_point float_type>
bool Lattice<Descriptor, float_type>::isAtBound(int index, int direction, float_type &rho_b, float_type &u_b) {
    // Implement boundary checking logic here
    // For example, check if index corresponds to a wall or lid
    // Set rho_b and u_b accordingly
    return false; // Placeholder return value
}

// Initialize equilibrium with u=0 and rho=rho_init, called only in the beginning of the simulation
template <isDescriptor Descriptor, std::floating_point float_type>
void Lattice<Descriptor, float_type>::initialize_equilibrium() {
    for (int cell = 0; cell < total_cells; ++cell) {
        for (int i = 0; i < q; ++i) {
            // f_eq_i = w_i * rho (when u = 0)
            f[i][cell] = Descriptor::w[i] * rho_init;
        }
    }
}

// Write a scalar field vector to VTK file for visualization
template <isDescriptor Descriptor, std::floating_point float_type>
void Lattice<Descriptor, float_type>::write_vtk(const std::vector<float_type>& data, const std::string& field_name, int iter) {
    const std::string filename = output_file + "-" + std::to_string(iter) + ".vtk";

    std::ofstream file(filename);
    if (!file.is_open()) return;

    file << "# vtk DataFile Version 2.0\n";
    file << "LBM Scalar Field\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";

    if (Descriptor::d == 2) {
        file << "DIMENSIONS " << sizes[0] << " " << sizes[1] << " 1\n";
        file << "ORIGIN 0 0 0\n";
        file << "SPACING " << dx << " " << dx << " " << dx << "\n";
    } else {
        file << "DIMENSIONS " << sizes[0] << " " << sizes[1] << " " << sizes[2] << "\n";
        file << "ORIGIN 0 0 0\n";
        file << "SPACING " << dx << " " << dx << " " << dx << "\n";
    }

    file << "POINT_DATA " << total_cells << "\n";
    file << "SCALARS " << field_name << " float\n";
    file << "LOOKUP_TABLE default\n";

    for (auto val : data) {
        file << val << "\n";
    }

    file.close();
}

// Write a vector field to VTK file for visualization
template <isDescriptor Descriptor, std::floating_point float_type>
void Lattice<Descriptor, float_type>::write_vtk(const std::array<std::vector<float_type>, Descriptor::d>& data, const std::string& field_name, int iter) {
    const std::string filename = output_file + "-" + std::to_string(iter) + ".vtk";

    std::ofstream file(filename);
    if (!file.is_open()) return;

    file << "# vtk DataFile Version 2.0\n";
    file << "LBM Vector Field\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";

    if (Descriptor::d == 2) {
        file << "DIMENSIONS " << sizes[0] << " " << sizes[1] << " 1\n";
        file << "ORIGIN 0 0 0\n";
        file << "SPACING " << dx << " " << dx << " " << dx << "\n";
    } else {
        file << "DIMENSIONS " << sizes[0] << " " << sizes[1] << " " << sizes[2] << "\n";
        file << "ORIGIN 0 0 0\n";
        file << "SPACING " << dx << " " << dx << " " << dx << "\n";
    }

    file << "POINT_DATA " << total_cells << "\n";
    file << "VECTORS " << field_name << " float\n";

    for (int cell = 0; cell < total_cells; ++cell) {
        for (int k = 0; k < Descriptor::d; ++k) {
            file << data[k][cell] << " ";
        }
        if (Descriptor::d == 2) {
            file << "0.0";
        }
        file << "\n";
    }

    file.close();
}

// The first type is the descriptor's weights and velocities
// THe second type is the lattice's dtributions and other datas
template class Lattice<D2Q9<float>, float>;
template class Lattice<D2Q9<double>, double>;
template class Lattice<D3Q19<float>, float>;
template class Lattice<D3Q19<double>, double>;

// function to compute the stess tensor of the lattice -- add in the Solver class
// sigma a,b = -(1-delta_t/2tau)sum(c_ia*c_ib*f_neq_i), with f_neq_i being f_i - f_eq