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
void force(std::vector<particle>&, float&);
void print_particle_info(std::vector<particle>);
float ranf();

int main(){
    std::vector<particle> particles = lattice_pos(5, 1);
    float energy = 0;
    initialise(particles, 273, 0.1);
    force(particles, energy);
    print_particle_info(particles);
    std::cout << "Energy = " << energy << std::endl;
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

void force(std::vector<particle> &particles, float &energy){
    int npart = particles.size();
    std::vector<float> r;
    std::vector<float> box = {10, 10, 10};
    float r2;
    float rc = 4.5;
    float ecut = 4 * ((1/std::pow(rc, 12)) - (1/std::pow(rc, 6)));
    for (int i = 0; i < npart-1; i++){
        for (int j = i+1; j < npart; j++){
            r[0] = particles[i].coordinates[0] - particles[j].coordinates[0];
            r[1] = particles[i].coordinates[1] - particles[j].coordinates[1];
            r[2] = particles[i].coordinates[2] - particles[j].coordinates[2];
            r[0] -= box[0]*std::nearbyint(r[0]/box[0]);
            r[1] -= box[1]*std::nearbyint(r[1]/box[1]);
            r[2] -= box[2]*std::nearbyint(r[2]/box[2]);
            r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
            float rc2 = rc*rc;
            if (r2 < rc2){
                float r2i = 1/r2;
                float r6i = std::pow(r2i, 3);
                float ff = 48*r2i*r6i*(r6i-0.5);
                particles[i].forces[0] += ff*r[0];
                particles[i].forces[1] += ff*r[1];
                particles[i].forces[2] += ff*r[2];
                particles[j].forces[0] -= ff*r[0];
                particles[j].forces[1] -= ff*r[1];
                particles[j].forces[2] -= ff*r[2];
                energy += 4*r6i*(r6i-1) - ecut;
            }
        }
    }
}

