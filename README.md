# Genome Rearrangement Simulator

A simple and fast simulator written in c++ with python bindings, Allowing advanced research in evolutionary biology.

## ⚡️ Quick start

- To compile the simulator c++ python package do the following:
  1. Clone this repo to a linux machine.
  2. Create a python 3.6+ environment.
  3. cd into the GenomeRearrangement_package folder
  4.  > ```pip install .```
   
- Example use of the simulator:
  
```python
import GenomeRearrangement as GR

tree_file = "path_to_newick_tree_file"
simulator = GR.Sim(tree_file)

sim_params = {
    "Chromosomes": [200,300,500],
    "AParam": 1.1,
    "maxBlockSize": 50,
    "InvertionRate": 0.5,
    "TranslocateRatio": 0.7,
    "FusionRate": 0.1,
    "FissionRate": 0.02,
    "DuplicationRate": 0.01,
    "LossRate":  0.001,
    "randRootAparam":  1.1
}
sim_params = list(sim_params.values())

simulator.init_sim(*sim_params)
simulator.set_seed(1)

simulated_genomes = simulator.run_sim()
```

- Generating summary statistics:
  
```python

genomes_object = GR.genomes(simulated_genome, simulator.get_tree())
sumarry_statistics = genomes_object.get_sum_stats()
```

