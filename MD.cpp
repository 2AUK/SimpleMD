#include <iostream>
#include <vector>
#include <tuple>

std::vector<float> lattice_pos(int, float);

int main(){
    std::vector<float> newc = lattice_pos(5, 0.5);
    for (auto it = newc.begin(); it != newc.end(); it++){
        float x = *std::next(it, 0); float y = *std::next(it, 1); float z = *std::next(it, 2);
        std::cout << x << '\t' << y << '\t' << z << '\n';
        std::advance(it, 2);
    }
    return 0;
}

//2D arrays [[0, 0, 0], [0, 0, 0.5], ..., [5, 5, 5]] for however many particles there are.
//Essentially building a lattice of particles. This would be easiest by building a cube 
//using whatever number evenly divides in to the number of particles
//probably the stupidest thing i've ever done
std::vector<float> lattice_pos(int npart, float grid_spacing){
    std::vector<float> coordinates;
    for (float x = 0; x <= npart*grid_spacing; x+= grid_spacing){
        for (float y = 0; y <= npart*grid_spacing; y+= grid_spacing){
            for (float z = 0; z <= npart*grid_spacing; z += grid_spacing){
                coordinates.insert(coordinates.end(), {x, y, z});
            }
        }
    }
    std::cout << "Lattice Gen Routine Done\n";
    return coordinates;
}