
from cc3d import CompuCellSetup
        

from NoisyChemotaxisSteppables import NoisyChemotaxisSteppable


CompuCellSetup.register_steppable(steppable=NoisyChemotaxisSteppable(frequency=1))


CompuCellSetup.run()

