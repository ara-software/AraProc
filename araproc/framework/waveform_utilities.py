import numpy as np
import ROOT

def tgraph_to_arrays(the_graph):
    xvals = np.asarray(the_graph.GetX())
    yvals = np.asarray(the_graph.GetY())
    return xvals, yvals
