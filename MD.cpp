#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <numeric>

/*Simple MD code for simple monoatomic liquids
 *Energy and Force calculated using Lennard-Jones
 */

struct particle{
    std::vector<float> coordinates;
    std::vector<float> velocities;
    std::vector<float> forces;
};

std::vector<particle> lattice_pos(int, float);
void initialise(std::vector<particle>&, float, float);
void print_particle_info(std::vector<particle>);
float ranf();

int main(){
    std::vector<particle> particles = lattice_pos(5, 1);
    initialise(particles, 273, 0.1);
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

float ranf(){
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<float> distribution(-0.5, 0.5);
    return distribution(generator);
}

void initialise(std::vector<particle> &particles, float temp, float dt){
    float sumv = 0;
    float sumv2 = 0;
    float npart = particles.size();
    for (particle &i : particles){
        i.velocities = {ranf(), ranf(), ranf()};
        sumv += std::accumulate(i.velocities.begin(), i.velocities.end(), 0.0);
        sumv2 += std::pow(std::accumulate(i.velocities.begin(), i.velocities.end(), 0.0), 2);
    }
    sumv /= npart;
    sumv2 /= npart;
    float fs = std::sqrt(3*temp/sumv2);
    for (particle &i : particles){
        i.velocities = {(i.velocities[0]-sumv)*fs, (i.velocities[1]-sumv)*fs, (i.velocities[2]-sumv)*fs};
        i.coordinates = {(i.coordinates[0] - i.velocities[0])*dt, (i.coordinates[1] - i.velocities[1])*dt, (i.coordinates[2] - i.velocities[2])*dt};
    }
}

