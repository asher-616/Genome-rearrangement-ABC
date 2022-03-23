#include "simulator.h"
#include "genomes.h"


int main() {

    std::vector<int> chromosomes = {123, 131, 102, 110, 117, 97, 88, 126, 106};
    Simulator sim("/home/elyawy/Development/Msc/projects/GR_pipeline/RAxML_tree.tree");
    sim.initSim(chromosomes, 1.03279829509959927, 50, 0.9657311105390779, 0.165151629627951, 0.8412344377699854, 0.18006701554605503, 0, 0, 1.01);
    sim.setSeed(1);
    sim.simulateBasedOnTree();

    return 0;
}