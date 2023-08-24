from brainsmash.mapgen.base import Base
from brainsmash.mapgen.eval import base_fit
import numpy as np

# The parameters were selected based on extensive trial and error

#base_fit(
#    x="C:/Users/Stuart/Documents/ThalamicGradients/data/preprocessed/PC1_thal.txt",
#    D="C:/Users/Stuart/Documents/ThalamicGradients/data/preprocessed/SeedDists.txt",
#    nsurr=100,
#    nh=25,
#    deltas=np.arange(0.1, 1, 0.1),
#    pv=75,
#    kernel='gaussian',
#    resample=False)

PCs = ["PC1","PC2","PC3"]

for PC in PCs:
	base = Base(x="./data/processed/"+PC+"_thal.txt", D="./data/preprocessed/SeedDists.txt", resample=False, kernel='gaussian', nh=25, deltas=np.arange(0.1, 1, 0.1), pv=75)
	surrogates = base(n=1000)
	np.savetxt("./data/processed/thal_surrogates_"+PC+".csv", surrogates, delimiter=",")
	print(PC)