#include <iostream>
#include <vector>

/*Simple MD code for simple monoatomic liquids
 *Energy and Force calculated using Lennard-Jones
 */

struct particle{
    std::vector<float> coordinates;
    std::vector<float> velocities;
    std::vector<float> forces;
};

std::vector<particle> lattice_pos(int, float);
void print_particle_info(std::vector<particle>);

int main(){
    std::vector<particle> particles = lattice_pos(5, 1);
    print_particle_info(particles);
    return 0;
}

//What kind of arcane trickery is involved in generating a 3D cube of points
//End me
std::vector<particle> lattice_pos(int npart, float grid_spacing){
    std::vector<particle> particles;
    int counter = 0;
    for (float x = 0; x < npart*grid_spacing; x+=grid_spacing){
        for (float y = 0; y < npart*grid_spacing; y+=grid_spacing){
            for (float z = 0; z < npart*grid_spacing; z +=grid_spacing){
                particles.push_back({{x, y, z}, {0, 0, 0}, {0, 0, 0}});
            }
        }
    }
    std::cout << "Lattice Gen Routine Done\n";
    return particles;
}

void print_particle_info(std::vector<particle> particles){
    int counter = 0;
    for (auto const i : particles){
        std::cout << "Particle " << counter << '\n'; 
        std::cout << "Coordinates\n";
        for (auto const coord : i.coordinates){
            std::cout << coord << '\t';
        }
        std::cout << "\nVelocities\n";
        for (auto const vels : i.velocities){
            std::cout << vels << '\t';
        }
        std::cout << "\nForces\n";
        for (auto const forces : i.forces){
            std::cout << forces << '\t';
        }
        std::cout << std::endl;
        counter++;
    }
}
void initialise(int nparts, std::vector<float>* coordinates, std::vector<float>* velocities){

}