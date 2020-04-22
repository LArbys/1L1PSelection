#!/bin/bash

container=/cluster/tufts/wongjiradlab/larbys/larbys-containers/singularity_dldependencies_pytorch1.3.sing

module load singularity
singularity shell --nv $container

