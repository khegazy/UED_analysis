import numpy as np
import sys
sys.path.append('../../../plots/scripts')
from plotClass import plotCLASS
from plotParams import pltParams


params = pltParams()
plc = plotCLASS()

#plc.print2d("../data/centering_smoothed[1024,1024].dat",
#   "../test_smoothed")
plc.print2d("../data/centering_selectedRange[1024,1024].dat",
    "../test_selectedRange")

