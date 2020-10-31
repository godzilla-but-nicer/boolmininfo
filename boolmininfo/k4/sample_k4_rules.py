import numpy as np

# get the samples
sample_size = 1000
num_combos = 2**(2**4)
sample = np.random.randint(num_combos, size=sample_size)

# write to file
np.savetxt(snakemake.output[0], sample)
