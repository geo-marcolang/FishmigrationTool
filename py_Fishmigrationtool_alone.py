#!/usr/bin/env python
# ------------------------------------------------------------------------------
__author__ = 'Marco Lang'
__contact__ = 'marco.lang1@students.unibe.ch'
__date__ = '05 april 2020'
__version__ = '1.1'
__status__ = "initial release"

"""
Name:           py_Fishmigrationtool_alone.py
Compatibility:  Python 3.x
Description:    ...

Requires:       tkinter, numpy, pandas, networkx, matplotlib

AUTHOR:         Marco Lang
ORGANIZATION:   University Bern
"""
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
from tkinter.ttk import Combobox
import os
import sys
import math
import networkx as nx
import matplotlib.pyplot as plt
import itertools
import numpy as np
from numpy import random
import re
import os
import sys
from scipy.spatial import distance
import shutil
import copy
import pandas as pd
import csv
import re
from itertools import islice
from matplotlib.pyplot import cm
import matplotlib
matplotlib.use('Agg')
# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------
# ---------------------------------____Alone____------------------------------
def createGraph(edgfilepath, twodmpath, ndsdepthsol, ndsabsvelocitysol, ndswsesol):
    """
    # Function definition:
        This function reads the .edg edge file and the .2dm mesh file and creates first an graph with nodes and edges.
        Then the attributes from the sol files (depth, abs_velocity and wse) are added to the nodes
    # Input (1) needed:
        full filename and filepath of .edg file
    # Input (2) needed:
        full filename and filepath of .2dm file
    # Input (3) needed:
        full filename and filepath of nds_depth.sol file
    # Input (4) needed:
        full filename and filepath of nds_abs_velocity.sol file
    # Input (5) needed:
        full filename and filepath of nds_wse.sol file
    # Output created:
        returns the networkx graph and an array of mesh nodes with: first column is: NODEID, second, third, fourth column is: x,y,z
    """
    # Create an empty graph with no nodes and no edges
    G = nx.Graph()

    # Load the .edg data (Input (1))
    with open(edgfilepath, 'r') as edg:
        lines = edg.readlines()  # type: List[str]

    # The data, which we need to calculate are only from the sixth line. That's why we skip the lines in front of it
    edgdata = lines[5:]

    # The lines are separated by a whitespace, so you can split them into rows on this whitespace
    columns = [line.split() for line in edgdata]

    # Because the columns are of data type integer, they must be converted into a float
    Edge = [int(column[0]) for column in columns]
    Node1 = [int(column[1]) for column in columns]

    # Make a copy of Node1 while expanding the list of Node2 values further down
    Elem_L = [int(column[3]) for column in columns]
    Elem_R = [column[4] for column in columns]
    Node1copy = Node1
    Node2 = [int(column[2]) for column in columns]

    # In the Elem_R series, the NULL values are replaced with -9999.
    Elem_R = [rep.replace('NULL', '-9999') for rep in Elem_R]

    # Now the list of Node2 is extended by the copied list Node1copy.
    Node1copy.extend(Node2)

    # Convert strings in List Elem_R to integer
    Elem_R = map(int, Elem_R)

    # Add leading zeros to my data
    [str(item).zfill(5) for item in Edge]

    # To create edges, the two node rows must be combined in a single tuple. => Edges=(Node1_line1,Node2_line1),(Node1_line2,Node2_line2),(Node1_line3,Node2_line3)
    tuple_N = zip(Node1, Node2)

    # Setting Node1 back to the original list
    Node1 = [int(column[1]) for column in columns]

    # Creating Edges and check the graph info
    G.add_edges_from(tuple_N, weight=float(1.00))

    # Use .2dm file for the position of the nodes (Input (2))
    mesh = open(twodmpath, "r")
    countnodes = 0
    for line in mesh:
        checknd = line.strip().split()
        if checknd[0] == 'ND' and len(checknd) == 5:
            countnodes = countnodes + 1
    nodesarray = np.zeros((countnodes, 4), dtype=float)
    mesh.close()
    mesh = open(twodmpath, "r")
    i = 0
    for line2 in mesh:
        checknd = line2.strip().split()
        if checknd[0] == 'ND' and len(checknd) == 5:
            nodesarray[i, 0] = checknd[1]
            nodesarray[i, 1] = checknd[2]
            nodesarray[i, 2] = checknd[3]
            nodesarray[i, 3] = checknd[4]
            i = i + 1
    mesh.close()

    # Extract nodes, x- and y-position from attribute table from the .2dm file
    Attribute_nodes = nodesarray[:, 0]
    Attribute_xPos = nodesarray[:, 1]
    Attribute_yPos = nodesarray[:, 2]
    Attribute_pos = zip(Attribute_xPos, Attribute_yPos)
    Attribute_pos_ID = zip(Attribute_pos, Attribute_nodes)

    # Adding nodes to Graph
    G.add_nodes_from(Attribute_nodes)

    # Creat an empty dictionary
    emd = {}

    # Fill dictionary with nodes as keys and xy-position as values
    emd = dict(zip(Attribute_nodes, Attribute_pos))

    for n, p in emd.items():
        G.nodes[n]['pos'] = p

    # show graph
    # nx.draw(G, pos=emd, node_size=10)
    # plt.savefig
    # plt.close("all")

    # Make the attributes of the function available for other functions
    createGraph.Attribute_nodes = Attribute_nodes
    createGraph.Attribute_xPos = Attribute_xPos
    createGraph.Attribute_yPos = Attribute_yPos
    createGraph.Attribute_pos = Attribute_pos
    createGraph.nodesarray = nodesarray

    del (mesh)
    del (countnodes)
    del (line)
    del (checknd)
    del (i)

    # Load the nds_depth.sol data (Input (3))
    meshd = open(twodmpath, "r")
    ndsdepth = open(ndsdepthsol, "r")
    # Count nodes in mesh
    countnodesd = 0
    nodeslistd = []
    for line in meshd:
        checkndd = line.strip().split()
        if checkndd[0] == 'ND':
            nodeslistd.append(checkndd[1])
            countnodesd += 1
    meshd.close()
    # Count timesteps in nds_depth file
    counttimestepsd = 0
    timestepslistd = []
    for line in ndsdepth:
        checkndd = line.strip().split()
        if checkndd[0] == 'TS' and len(checkndd) == 3:
            if checkndd[2] == '0x0.000000000000p+0' or '0x0p+0':
                timestepslistd.append(0.0)
            else:
                timestepslistd.append(float(checkndd[2]))
            counttimestepsd += 1
    ndsdepth.close()
    rowsd = len(nodeslistd)
    columsd = len(timestepslistd) + 2  # one column added first for NODEID and one column at the end for max(flowdepth)
    # Create the array with zeros
    ndsdeptharray = np.zeros((rowsd, columsd), dtype=float)
    i = 0

    for item in nodeslistd:
        ndsdeptharray[i, 0] = float(item)
        i += 1

    # Loop though nds_depth.sol file and add flowdepths of each timestep to ndsdeptharray
    ndsdepth = open(ndsdepthsol, "r")
    linecounterd = 0
    columncounterd = 0
    for line in ndsdepth:
        checkndd = line.strip().split()
        if len(checkndd) < 1:
            continue
        if checkndd[0] in ["MESH2d", "MESH2D", "MESHNAME", "SCALAR", "ND", "DATASET", "OBJTYPE", "RT_JULIAN", "BEGSCL",
                           "NC", "NAME", "TIMEUNITS", "seconds"]:
            continue
        if checkndd[0] == 'TS' and len(checkndd) > 1:
            columncounterd += 1
            linecounterd = 0
        elif len(checkndd) == 1 and checkndd[0] != "ENDDS":
            ndsdeptharray[linecounterd, columncounterd] = float(checkndd[0])
            linecounterd += 1
    ndsdepth.close()
    # Loop through the rows in the last column of the array and write the max of flowdepth values
    i = 0
    while i < linecounterd:
        ndsdeptharray[i, -1] = max(ndsdeptharray[i, 1:-2])
        i += 1

    ndsdeptharraytranspose = ndsdeptharray.transpose()

    # Create outer dictionary with TS? - TS? as keys
    Stepsd = np.arange(start=0, stop=counttimestepsd, step=1)
    Timestepsd = []
    for i in Stepsd:
        Timestepsd.append("TS" + str(i))
        i = i + 1

    # Maxflowdepth as timestep at the end of the list
    Timestepsd.append("MAX")
    Timestepsd = np.insert(Timestepsd, 0, "NODE")

    TS_tupled = tuple(Timestepsd)

    Depth_TS_out = {}
    Depth_TS_out = dict.fromkeys(Timestepsd)

    # Make a tuple of the values
    depthvalues_tuple = tuple(ndsdeptharraytranspose)

    # Creating dictionary with timesteps as keys and depthvalues as values because networkx needs one
    Depth_TS = {}
    Depth_TS = dict(zip(TS_tupled, depthvalues_tuple))

    # Adding Depth_TS as attributes to the graph
    for timestepsd in Depth_TS:
        nx.set_node_attributes(G, Depth_TS[timestepsd], 'Water_Depth')

    # nx get nodes attributes
    ndsdepthattributes = nx.get_node_attributes(G, 'Water_Depth')

    # Load the nds_depth.sol data (Input (4))
    meshv = open(twodmpath, "r")
    ndsabsvelocity = open(ndsabsvelocitysol, "r")
    # Count nodes in mesh
    countnodesv = 0
    nodeslistv = []
    for line in meshv:
        checkndv = line.strip().split()
        if checkndv[0] == 'ND':
            nodeslistv.append(checkndv[1])
            countnodesv += 1
    meshv.close()
    # Count timesteps in nds_velocity file
    counttimestepsv = 0
    timestepslistv = []
    for line in ndsabsvelocity:
        checkndv = line.strip().split()
        if checkndv[0] == 'TS' and len(checkndv) == 3:
            if checkndv[2] == '0x0.000000000000p+0' or '0x0p+0':
                timestepslistv.append(0.0)
            else:
                timestepslistv.append(float(checkndv[2]))
            counttimestepsv += 1
    ndsabsvelocity.close()
    rowsv = len(nodeslistv)
    columsv = len(
        timestepslistv) + 2  # one column added first for NODEID and one column at the end for max(flow velocity)
    # Create the array with zeros
    ndsabsvelocityarray = np.zeros((rowsv, columsv), dtype=float)
    i = 0

    for item in nodeslistv:
        ndsabsvelocityarray[i, 0] = float(item)
        i += 1

    # Loop though nds_velocity.sol file and add flow velocity of each timestep to ndsabsvelocityarray
    ndsabsvelocity = open(ndsabsvelocitysol, "r")
    linecounterv = 0
    columncounterv = 0
    for line in ndsabsvelocity:
        checkndv = line.strip().split()
        if len(checkndv) < 1:
            continue
        if checkndv[0] in ["MESH2d", "MESH2D", "MESHNAME", "SCALAR", "ND", "DATASET", "OBJTYPE", "RT_JULIAN", "BEGSCL",
                           "NC", "NAME", "TIMEUNITS", "seconds"]:
            continue
        if checkndv[0] == 'TS' and len(checkndv) > 1:
            columncounterv += 1
            linecounterv = 0
        elif len(checkndv) == 1 and checkndv[0] != "ENDDS":
            ndsabsvelocityarray[linecounterv, columncounterv] = float(checkndv[0])
            linecounterv += 1
    ndsabsvelocity.close()
    # Loop through the rows in the last column of the array and write the max of flow velocity values
    i = 0
    while i < linecounterv:
        ndsabsvelocityarray[i, -1] = max(ndsabsvelocityarray[i, 1:-2])
        i += 1

    ndsabsvelocityarraytranspose = ndsabsvelocityarray.transpose()

    # Create outer dictionary with TS? - TS? as keys
    Stepsv = np.arange(start=0, stop=counttimestepsv, step=1)
    Timestepsv = []
    for i in Stepsv:
        Timestepsv.append("TS" + str(i))
        i = i + 1

    # Maxflowvelocity as timestep at the end of the list
    Timestepsv.append("MAX")
    Timestepsv = np.insert(Timestepsv, 0, "NODE")

    TS_tuplev = tuple(Timestepsv)

    Velocity_TS_out = {}
    Velocity_TS_out = dict.fromkeys(Timestepsv)

    # Make a tuple of the values
    velocityvalues_tuple = tuple(ndsabsvelocityarraytranspose)

    # Creating dictionary with timesteps as keys and velocityvalues as values because networkx needs one
    Velocity_TS = {}
    Velocity_TS = dict(zip(TS_tuplev, velocityvalues_tuple))

    # Adding Velocity_TS as attributes to the graph
    for timestepsv in Velocity_TS:
        nx.set_node_attributes(G, Velocity_TS[timestepsv], 'Flow_Velocity')

    # nx get nodes attributes
    ndsvelocityattributes = nx.get_node_attributes(G, 'Flow_Velocity')

    # Load the nds_wse.sol data (Input (5))
    meshw = open(twodmpath, "r")
    ndswse = open(ndswsesol, "r")
    # Count nodes in mesh
    countnodesw = 0
    nodeslistw = []
    for line in meshw:
        checkndw = line.strip().split()
        if checkndw[0] == 'ND':
            nodeslistw.append(checkndw[1])
            countnodesw += 1
    meshw.close()
    # Count timesteps in nds_wse file
    counttimestepsw = 0
    timestepslistw = []
    for line in ndswse:
        checkndw = line.strip().split()
        if checkndw[0] == 'TS' and len(checkndw) == 3:
            if checkndw[2] == '0x0.000000000000p+0' or '0x0p+0':
                timestepslistw.append(0.0)
            else:
                timestepslistw.append(float(checkndw[2]))
            counttimestepsw += 1
    ndswse.close()
    rowsw = len(nodeslistw)
    columsw = len(timestepslistw) + 2  # one column added first for NODEID and one column at the end for max(wse)
    # Create the array with zeros
    ndswsearray = np.zeros((rowsw, columsw), dtype=float)
    i = 0

    for item in nodeslistw:
        ndswsearray[i, 0] = float(item)
        i += 1

    # Loop though nds_wse.sol file and add wse of each timestep to ndswsearray
    ndswse = open(ndswsesol, "r")
    linecounterw = 0
    columncounterw = 0
    for line in ndswse:
        checkndw = line.strip().split()
        if len(checkndw) < 1:
            continue
        if checkndw[0] in ["MESH2d", "MESH2D", "MESHNAME", "SCALAR", "ND", "DATASET", "OBJTYPE", "RT_JULIAN", "BEGSCL",
                           "NC", "NAME", "TIMEUNITS", "seconds"]:
            continue
        if checkndw[0] == 'TS' and len(checkndw) > 1:
            columncounterw += 1
            linecounterw = 0
        elif len(checkndw) == 1 and checkndw[0] != "ENDDS":
            ndswsearray[linecounterw, columncounterw] = float(checkndw[0])
            linecounterw += 1
    ndswse.close()
    # Loop through the rows in the last column of the array and write the max of wse values
    i = 0
    while i < linecounterw:
        ndswsearray[i, -1] = max(ndswsearray[i, 1:-2])
        i += 1

    ndswsearraytranspose = ndswsearray.transpose()

    # Create outer dictionary with TS? - TS? as keys
    Stepsw = np.arange(start=0, stop=counttimestepsw, step=1)
    Timestepsw = []
    for i in Stepsw:
        Timestepsw.append("TS" + str(i))
        i = i + 1

    # Maxwse as timestep at the end of the list
    Timestepsw.append("MAX")
    Timestepsw = np.insert(Timestepsw, 0, "NODE")

    TS_tuplew = tuple(Timestepsw)

    Wse_TS_out = {}
    Wse_TS_out = dict.fromkeys(Timestepsw)

    # Make a tuple of the values
    wsevalues_tuple = tuple(ndswsearraytranspose)

    # Creating dictionary with timesteps as keys and wsevalues as values because networkx needs one
    Wse_TS = {}
    Wse_TS = dict(zip(TS_tuplew, wsevalues_tuple))

    ## Adding Wse_TS as attributes to the graph
    for timestepsw in Wse_TS:
        nx.set_node_attributes(G, Wse_TS[timestepsw], 'WSE')

    # nx get nodes attributes
    ndswseattributes = nx.get_node_attributes(G, 'WSE')

    # Make the attributes of the function available for other functions
    createGraph.G = G
    createGraph.emd = emd
    createGraph.Velocity_TS = Velocity_TS
    createGraph.Depth_TS = Depth_TS
    createGraph.Wse_TS = Wse_TS

    print("Das Netzwerk wurde erfolgreich erstellt und  Wassertiefe, Fliessgeschwindigkeit und Wasserspiegellage wurden als Attribute angehängt")
    #Infos about the graph
    print("Folgende Informationen zum gesamten Netzwerk")
    print(nx.info(G))
# END - creating Graph with flow velocity and water depth attributes
# ------------------------------------------------------------------------------
def massconservation(inflowdat,outflowdat, prozentabweichung, plotformat, outfiledirectory):
    """
    # Function definition:
        The Inflow_th.dat and Outflow_th.dat files are used to check the mass conservation of the model.
        After a certain time, the outflow compensates the continuous inflow.
    # Input (1) needed:
        full filename and filepath of inflow.dat file
    # Input (2) needed:
        full filename and filepath of outflow.dat file
    # Input (3) needed:
        format of the output graphic
    # Input (4) needed:
        filepath for the outputdirectory
    # Output created:
        returns the networkx graph and an array of mesh nodes with: first column is: NODEID, second, third, fourth column is: x,y,z
    """
    #load inputdata Input (1)
    inflow = inflowdat
    # load inputdata Input (2)
    outflow = outflowdat

    inflowr = np.loadtxt(inflow)
    outflowr = np.loadtxt(outflow)

    #Steady inflow hydrograph and outflow hydrograph
    x = inflowr[:, 0]
    y1 = inflowr[:, 1]
    y2 = outflowr[:, 1]
    fig, ax = plt.subplots()

    sumtimesteps = len(outflowr[:, 1])
    timeoftimestepts = inflowr[:, 0]

    timesteps = []
    bigger = []
    smaller = []
    index = len(outflowr[:, 1]) - 1
    index1 = len(outflowr[:, 1]) - 1

    i = 0
    while i <= index:
        timesteps.append(i)
        i = i + 1
    #timesteps.insert(0,0)
    prozent = prozentabweichung
    prozentsatz = (100-prozent)/100

    i = 1
    while i <= index1:
        if outflowr[:, 1][i] >= prozentsatz * (inflowr[:, 1][i]): #It will reflect the time steps of the outflow which deviate by a maximum of 1 percent from the inflow.
            bigger.append("TS" + str(i))
        else:
            smaller.append("TS" + str(i))
        i = i + 1
    print('Es gibt ' + str(sumtimesteps) + ' Timesteps für diese Simulation')
    if bigger == 0:
        print('Kein Timestep ist in einem stationären Zustand zwischen Einfluss (inflow) und Abfluss (outflow)')
    else:
        print('Die folgenden Timesteps sind in einem stationären Zustand zwischen Einfluss (inflow) und Abfluss (outflow): ' + str(bigger))

    plt.title("Kontinuierliche Zufluss- und Abflussganglinie",fontsize=10)
    ax.set_xlabel('Zeit [s]',fontsize=6)
    ax.set_ylabel('Abfluss [$m^3/s$]',fontsize=6)
    ax.plot(x, y1, label="Q inflow")
    ax.plot(x, y2, label="Q outflow")
    ax.legend(loc='best', frameon=False,fontsize=6)
    plt.grid(True)
    #plt.text(0.95, 0.05, 'Die folgenden Timesteps sind in einem stationären Zustand zwischen Einfluss (inflow) und Abfluss (outflow):\n' + str(bigger), horizontalalignment='right',verticalalignment='bottom', fontsize= 4,transform=ax.transAxes)
    ax3 = ax.twiny()
    ax3.set_xlabel('Timesteps', color="black",fontsize=6)
    ax3.plot(timesteps, y1, label='Q inflow')
    ax3.tick_params(axis='x')

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()
    plt.savefig(outfiledirectory + '/' + 'Massenkonservierung des Modells - inflow Hydrograph und outflow Hydrograph' + '.'+ plotformat)
    plt.close("all")
# END - mass conservation of inflow and outflow hydrograph
# ------------------------------------------------------------------------------
def sourcenodes(timestep):
    """
    # Function definition:
        This function takes the previously built network (createGraph) and looks at a neighborhood analysis which of the entered nodes has the greatest water depth at the entered timestep.
        Between source and targetnode the thalweg is searched in a following function.
    # Input (1) needed:
        just a timestep where we have steady state conditions.
    # Input (2) needed:
        nodes that can be used as sourcenodes
    # Output created:
        returns the node in the network that has the greatest water depth based on a neighborhood analysis of the nodes entered.
    """
    # Load the networkx graph from the previous function
    F = nx.Graph(createGraph.G)
    # List of the sourcenodes
    nodesremovelist = startsourcetargetnodes.nodesremoveint
    sourcenodeslist = startsourcetargetnodes.nodessint

    F.remove_nodes_from(nodesremovelist)
    print("Folgende Nodes wurden aus dem Netzwerk gelöscht" + str(nodesremovelist))
    #H = nx.Graph(createGraph.G)
    #Attribute_nodes = createGraph.Attribute_nodes
    #depth = dict(zip(Attribute_nodes, Depth_TS[("TS" + str(timestep))]))
    #nx.set_node_attributes(H, depth, 'Water_Depth')
    Depth_TS = createGraph.Depth_TS
    nodesarray = createGraph.nodesarray

    #Adding waterdepth to nodesarray
    depthttimestep = Depth_TS[("TS" + str(timestep))]
    depthttimesteptranspose = depthttimestep.transpose()
    depthtoarray = np.insert(nodesarray, 4, values=depthttimesteptranspose, axis=1)

    sourcenodeslist.sort()
    neighbourmeandepth = []

    #Analyzes which of the entered nodes has the greatest water depth, based on the nearest neighbor.
    for l in sourcenodeslist:
        try:
                idx=[]
                idx.append([n for n in F.neighbors(l)])
                neigbournodesarray = depthtoarray[np.in1d(depthtoarray[:, 0], idx)]
                sumdepth = np.mean(neigbournodesarray[:,4])
                neighbourmeandepth.append(sumdepth)
        except nx.NetworkXError:
            print('Der Node ' + str(l) + ' ist nicht im Graph')
            continue
        except nx.NodeNotFound:
            print('Sourcenode ' + str(l) + ' ist nicht im Graph')
            continue
        except KeyError:
            print('Node ' + str(l) + ' ist nicht verfuegbar')
            continue

    meanlist = neighbourmeandepth
    meanlistarray = np.asarray(meanlist)
    meanlistarraytranspose = meanlistarray.transpose()
    idx1 = sourcenodeslist
    sourcenodesarray = depthtoarray[np.in1d(depthtoarray[:, 0], idx1)]
    sourcenodesarrayneighbour = np.insert(sourcenodesarray, 5, values=meanlistarraytranspose, axis=1)
    maxdepthindex = np.where(sourcenodesarrayneighbour[:, 5] == np.amax(sourcenodesarrayneighbour[:, 5]))
    #print(sourcenodesarrayneighbour)
    print("Starte den Talweg mit dem Node " + str(int(sourcenodesarrayneighbour[:, 0][maxdepthindex])))
    sourcenodes.source = sourcenodesarrayneighbour[:, 0][maxdepthindex]
# END - finding the deepest nodes in the inflow string
# ------------------------------------------------------------------------------
def targetnodes(timestep):
    """
    # Function definition:
        This function takes the previously built network (createGraph) and looks at a neighborhood analysis which of the entered nodes has the greatest water depth at the entered timestep.
        Between source and targetnode the thalweg is searched in a following function.
    # Input (1) needed:
        just a timestep where we have steady state conditions.
    # Input (2) needed:
        nodes that can be used as targetnodes
    # Output created:
        returns the node in the network that has the greatest water depth based on a neighborhood analysis of the nodes entered.
    """
    # Load the networkx graph from the previous function
    F = nx.Graph(createGraph.G)
    # List of the sourcenodes
    targetnodeslist = startsourcetargetnodes.nodestint
    nodesremovelist = startsourcetargetnodes.nodesremoveint

    F.remove_nodes_from(nodesremovelist)
    #H = nx.Graph(createGraph.G)
    #Attribute_nodes = createGraph.Attribute_nodes
    #depth = dict(zip(Attribute_nodes, Depth_TS[("TS" + str(timestep))]))
    #nx.set_node_attributes(H, depth, 'Water_Depth')
    Depth_TS = createGraph.Depth_TS
    nodesarray = createGraph.nodesarray

    #Adding waterdepth to nodesarray
    depthttimestep = Depth_TS[("TS" + str(timestep))]
    depthttimesteptranspose = depthttimestep.transpose()
    depthtoarray = np.insert(nodesarray, 4, values=depthttimesteptranspose, axis=1)

    targetnodeslist.sort()
    neighbourmeandepth = []

    #Analyzes which of the entered nodes has the greatest water depth, based on the nearest neighbor.
    for l in targetnodeslist:
        try:
            idx=[]
            idx.append([n for n in F.neighbors(l)])
            neigbournodesarray = depthtoarray[np.in1d(depthtoarray[:, 0], idx)]
            sumdepth = np.sum(neigbournodesarray[:,4])
            neighbourmeandepth.append(sumdepth)
        except nx.NetworkXError:
            print('Der Node ' + str(l) + ' ist nicht im Graph')
            continue
        except nx.NodeNotFound:
            print('Targetnode ' + str(l) + ' ist nicht im Graph')
            continue
        except KeyError:
            print('Node ' + str(l) + ' ist nicht verfuegbar')
            continue

    meanlist = neighbourmeandepth
    meanlistarray = np.asarray(meanlist)
    meanlistarraytranspose = meanlistarray.transpose()
    idx1 = targetnodeslist
    targetnodesarray = depthtoarray[np.in1d(depthtoarray[:, 0], idx1)]
    targetnodesarrayneighbour = np.insert(targetnodesarray, 5, values=meanlistarraytranspose, axis=1)
    maxdepthindex = np.where(targetnodesarrayneighbour[:, 5] == np.amax(targetnodesarrayneighbour[:, 5]))
    #print(targetnodesarrayneighbour)
    print("Vollende den Talweg mit dem Node " + str(int(targetnodesarrayneighbour[:, 0][maxdepthindex])))
    targetnodes.target = targetnodesarrayneighbour[:, 0][maxdepthindex]
# END - finding the deepest nodes in the outflow string
# ------------------------------------------------------------------------------
def thalweg(outfiledirectory, inflowdat, outflowdat, mindestwassertiefe_generell, mindestwassertiefe_untiefen, mindestwassertiefe_einzelfall, maximallpassuntiefen, maximallpasseinzelfall, maximalgeschwindigkeit, maximalabsturzhoehe, sourcenode, targetnode, plotformat):
    """
    # Function definition:
        In this function the thalweg between two nodes is calculated. This can be done for several time steps.
        After a successful search of the thalweg, different graphics are plotted along the deepest channel.
    # Input (1) needed:
        filepath for the outputdirectory
    # Input (2) needed:
        full filename and filepath of inflow.dat file
    # Input (3) needed:
        full filename and filepath of outflow.dat file
    # Input (4) needed:
        general minimum water depth for fish
    # Input (5) needed:
        minimum water depth for shallows for fish
    # Input (6) needed:
        maximum velocity for fish
    # Input (7) needed:
        maximum fall height for fish, e.g. at threshold buoats
    # Input (8) needed:
        sourcenode
    # Input (9) needed:
        general minimum water depth for fish
    # Input (10) needed:
        targetnode
    # Input (11) needed:
        memory format for graphics e.g. pdf or png
    # Input (12) needed:
        timesteps
    # Output created:
        returns graphs and a csv file with the cumulative distribution of water depths, flow velocities and fall heights along the thalweg.
        As well as the water depths and flow velocities in the longitudinal profile of this deepest channel.
        In addition, the frequency and length of the undershoots of the minimum water depth.
    """
    # creat list with the desired timesteps
    #global H
    thalwegtimestepslist = startthalweg.timestepsint
    thalwegtimestepslistnew = np.array([], dtype=int)
    bedingung = startthalweg.bedingung
    talweg = startthalweg.talweg
    dischargelist = []
    waterdepthlist = []
    flowvelocitylist = []
    inflowdischargelist = []
    laengelist = []
    waterdepthprozentdislist = []
    waterdepthcumlist = []
    flowvelocityprozentdislist = []
    flowvelocitycumlist = []
    deltazprozentdislist = []
    deltazcumlist = []
    lengthkmglist = []
    lengthkmulist = []
    dischargelengthkmglist = []
    dischargelengthkmulist = []

    allpaths = []

    # Starts the loop to output the thalweg as csv for each desired timestep
    for i in thalwegtimestepslist:
        try:
            print("Berechnungen für den Timestep " + str(i))
            H = nx.Graph(createGraph.G)
            emd = createGraph.emd
            Attribute_nodes = createGraph.Attribute_nodes
            Depth_TS = createGraph.Depth_TS
            Velocity_TS = createGraph.Velocity_TS
            Wse_TS = createGraph.Wse_TS
            nodesarray = createGraph.nodesarray
            depth = dict(zip(Attribute_nodes, Depth_TS[("TS" + str(i))]))
            nx.set_node_attributes(H, depth, 'Water_Depth')

            # Change the weight of the edges
            nx.get_node_attributes(H, 'Water_Depth')
            H.edges(data=True)

            # add velocity to nodes
            velocity = dict(zip(Attribute_nodes, Velocity_TS[("TS" + str(i))]))
            nx.set_node_attributes(H, velocity, 'Flow_Velocity')

            # add wse to nodes
            wse = dict(zip(Attribute_nodes, Wse_TS[("TS" + str(i))]))
            nx.set_node_attributes(H, wse, 'WSE')

            if bedingung == 'mitBedingungen' and talweg == "tiefsterTalweg":
                print("Tiefster Talweg wird mit Bedingungen berechnet")
                # Delte edges with maximalgeschwindigkeit
                for u, v, f in H.edges(data=True):
                    if maximalgeschwindigkeit != 0:
                        if abs(((H.nodes[u]['Flow_Velocity']) + (
                        H.nodes[v]['Flow_Velocity'])) / 2) < maximalgeschwindigkeit:
                            f['weightv'] = abs(((H.nodes[u]['Flow_Velocity']) + (H.nodes[v]['Flow_Velocity'])) / 2)
                        else:
                            f['weightv'] = 10000000000000.0
                    elif maximalgeschwindigkeit == 0:
                        f['weightv'] = (((H.nodes[u]['Flow_Velocity']) + (H.nodes[v]['Flow_Velocity'])) / 2)

                # Delete edges with a weight of 1000000000000 otherwise the path will take the path with the fewest nodes at the edge because the distances are larger and it has fewer nodes.
                selected_edgesv = [(u, v) for u, v, f in H.edges(data=True) if f['weightv'] == 10000000000000.0]
                H.remove_edges_from(selected_edgesv)

                # Delte edges with maximalabsturzhoehe
                for u, v, g in H.edges(data=True):
                    if maximalabsturzhoehe != 0:
                        if abs(((H.nodes[u]['WSE']) - (H.nodes[v]['WSE'])) / math.sqrt(((((H.nodes[u]['pos'])[0]) - ((H.nodes[v]['pos'])[0])) ** 2) + (((((H.nodes[u]['pos'])[1]) - ((H.nodes[v]['pos'])[1]))) ** 2))) < maximalabsturzhoehe:
                            g['weighta'] = abs(((H.nodes[u]['WSE']) - (H.nodes[v]['WSE'])) / math.sqrt(((((H.nodes[u]['pos'])[0]) - ((H.nodes[v]['pos'])[0])) ** 2) + (((((H.nodes[u]['pos'])[1]) - ((H.nodes[v]['pos'])[1]))) ** 2)))
                        else:
                            g['weighta'] = 10000000000000.0
                    elif maximalabsturzhoehe == 0:
                        g['weighta'] = abs(((H.nodes[u]['WSE']) - (H.nodes[v]['WSE'])) / math.sqrt(((((H.nodes[u]['pos'])[0]) - ((H.nodes[v]['pos'])[0])) ** 2) + (((((H.nodes[u]['pos'])[1]) - ((H.nodes[v]['pos'])[1]))) ** 2)))
                # Delete edges with a weight of 1000000000000 otherwise the path will take the path with the fewest nodes at the edge because the distances are larger and it has fewer nodes.
                selected_edgesw = [(u, v) for u, v, g in H.edges(data=True) if g['weighta'] == 10000000000000.0]
                H.remove_edges_from(selected_edgesw)

                # 1/depth = x and then x/depth = y and again y/depth and so on to give the edges the right weighting.
                for u, v, d in H.edges(data=True):
                    # d['weight'] = ((G.nodes[u]['Water_Depth'])+(G.nodes[v]['Water_Depth']))/2
                    if mindestwassertiefe_einzelfall != 0:
                        if (H.nodes[u]['Water_Depth']) > 0.005 and (H.nodes[v]['Water_Depth']) > 0.005 and ((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth']) / 2) >= mindestwassertiefe_einzelfall:
                            zwg1 = 1 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg2 = zwg1 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg3 = zwg2 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg4 = zwg3 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg5 = zwg4 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg6 = zwg5 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg7 = zwg6 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg8 = zwg7 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg9 = zwg8 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg10 = zwg9 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg11 = zwg10 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg12 = zwg11 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg13 = zwg12 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg14 = zwg13 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg15 = zwg14 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg16 = zwg15 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            d['weight'] = zwg16 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        else:
                            d['weight'] = 10000000000000.0
                    elif mindestwassertiefe_einzelfall == 0:
                        if (H.nodes[u]['Water_Depth']) > 0.005 and (H.nodes[v]['Water_Depth']) > 0.005 and ((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth']) / 2) >= 0.005:
                            zwg1 = 1 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg2 = zwg1 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg3 = zwg2 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg4 = zwg3 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg5 = zwg4 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg6 = zwg5 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg7 = zwg6 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg8 = zwg7 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg9 = zwg8 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg10 = zwg9 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg11 = zwg10 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg12 = zwg11 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg13 = zwg12 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg14 = zwg13 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg15 = zwg14 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg16 = zwg15 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            d['weight'] = zwg16 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        else:
                            d['weight'] = 10000000000000.0

                # Delete edges with a weight of 1000000000000 otherwise the path will take the path with the fewest nodes at the edge because the distances are larger and it has fewer nodes.
                selected_edges = [(u, v) for u, v, d in H.edges(data=True) if d['weight'] == 10000000000000.0]
                H.remove_edges_from(selected_edges)
                #H.edges(data=True)

            elif bedingung == 'ohneBedingungen' and talweg == "tiefsterTalweg":
                print("Tiefster Talweg wird ohne Bedingungen berechnet")
                # 1/depth = x and then x/depth = y and again y/depth and so on to give the edges the right weighting.
                for u, v, d in H.edges(data=True):
                    # d['weight'] = ((G.nodes[u]['Water_Depth'])+(G.nodes[v]['Water_Depth']))/2
                    if mindestwassertiefe_einzelfall != 0:
                        if (H.nodes[u]['Water_Depth']) > 0.005 and (H.nodes[v]['Water_Depth']) > 0.005 and ((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth']) / 2) >= 0.005:
                            zwg1 = 1 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg2 = zwg1 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg3 = zwg2 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg4 = zwg3 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg5 = zwg4 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg6 = zwg5 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg7 = zwg6 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg8 = zwg7 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg9 = zwg8 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg10 = zwg9 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg11 = zwg10 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg12 = zwg11 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg13 = zwg12 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg14 = zwg13 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg15 = zwg14 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg16 = zwg15 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            d['weight'] = zwg16 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        else:
                            d['weight'] = 10000000000000.0
                    elif mindestwassertiefe_einzelfall == 0:
                        if (H.nodes[u]['Water_Depth']) > 0.005 and (H.nodes[v]['Water_Depth']) > 0.005 and ((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth']) / 2) >= 0.005:
                            zwg1 = 1 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg2 = zwg1 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg3 = zwg2 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg4 = zwg3 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg5 = zwg4 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg6 = zwg5 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg7 = zwg6 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg8 = zwg7 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg9 = zwg8 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg10 = zwg9 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg11 = zwg10 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg12 = zwg11 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg13 = zwg12 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg14 = zwg13 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg15 = zwg14 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            zwg16 = zwg15 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                            d['weight'] = zwg16 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        else:
                            d['weight'] = 10000000000000.0
                # Delete edges with a weight of 1000000000000 otherwise the path will take the path with the fewest nodes at the edge because the distances are larger and it has fewer nodes.
                selected_edges = [(u, v) for u, v, e in H.edges(data=True) if e['weight'] == 10000000000000.0]
                H.remove_edges_from(selected_edges)
                #H.edges(data=True)
            elif bedingung == 'mitBedingungen' and talweg == "kuerzesterTalweg":
                print("Kürzester Talweg wird mit Bedingungen berechnet")
                # Delte edges with maximalgeschwindigkeit
                for u, v, f in H.edges(data=True):
                    if maximalgeschwindigkeit != 0:
                        if abs(((H.nodes[u]['Flow_Velocity']) + (
                        H.nodes[v]['Flow_Velocity'])) / 2) < maximalgeschwindigkeit:
                            f['weightv'] = abs(((H.nodes[u]['Flow_Velocity']) + (H.nodes[v]['Flow_Velocity'])) / 2)
                        else:
                            f['weightv'] = 10000000000000.0
                    elif maximalgeschwindigkeit == 0:
                        f['weightv'] = (((H.nodes[u]['Flow_Velocity']) + (H.nodes[v]['Flow_Velocity'])) / 2)

                # Delete edges with a weight of 1000000000000 otherwise the path will take the path with the fewest nodes at the edge because the distances are larger and it has fewer nodes.
                selected_edgesv = [(u, v) for u, v, f in H.edges(data=True) if f['weightv'] == 10000000000000.0]
                H.remove_edges_from(selected_edgesv)

                # Delte edges with maximalabsturzhoehe
                for u, v, g in H.edges(data=True):
                    if maximalabsturzhoehe != 0:
                        if abs(((H.nodes[u]['WSE']) - (H.nodes[v]['WSE'])) / math.sqrt(((((H.nodes[u]['pos'])[0]) - ((H.nodes[v]['pos'])[0])) ** 2) + (((((H.nodes[u]['pos'])[1]) - ((H.nodes[v]['pos'])[1]))) ** 2))) < maximalabsturzhoehe:
                            g['weighta'] = abs(((H.nodes[u]['WSE']) - (H.nodes[v]['WSE'])) / math.sqrt(((((H.nodes[u]['pos'])[0]) - ((H.nodes[v]['pos'])[0])) ** 2) + (((((H.nodes[u]['pos'])[1]) - ((H.nodes[v]['pos'])[1]))) ** 2)))
                        else:
                            g['weighta'] = 10000000000000.0
                    elif maximalabsturzhoehe == 0:
                        g['weighta'] = abs(((H.nodes[u]['WSE']) - (H.nodes[v]['WSE'])) / math.sqrt(((((H.nodes[u]['pos'])[0]) - ((H.nodes[v]['pos'])[0])) ** 2) + (((((H.nodes[u]['pos'])[1]) - ((H.nodes[v]['pos'])[1]))) ** 2)))
                # Delete edges with a weight of 1000000000000 otherwise the path will take the path with the fewest nodes at the edge because the distances are larger and it has fewer nodes.
                selected_edgesw = [(u, v) for u, v, g in H.edges(data=True) if g['weighta'] == 10000000000000.0]
                H.remove_edges_from(selected_edgesw)

                # 1/depth = x and then x/depth = y and again y/depth and so on to give the edges the right weighting.
                for u, v, d in H.edges(data=True):
                    if mindestwassertiefe_einzelfall != 0:
                        # d['weight'] = ((G.nodes[u]['Water_Depth'])+(G.nodes[v]['Water_Depth']))/2
                        if (H.nodes[u]['Water_Depth']) > 0.005 and (H.nodes[v]['Water_Depth']) > 0.005 and ((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth']) / 2) >= mindestwassertiefe_einzelfall:
                            d['weight'] = ((1/((H.nodes[u]['Water_Depth']))) + ((1/(H.nodes[v]['Water_Depth'])))) / 2
                        else:
                            d['weight'] = 10000000000000.0
                    elif mindestwassertiefe_einzelfall == 0:
                        # d['weight'] = ((G.nodes[u]['Water_Depth'])+(G.nodes[v]['Water_Depth']))/2
                        if (H.nodes[u]['Water_Depth']) > 0.005 and (H.nodes[v]['Water_Depth']) > 0.005 and ((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth']) / 2) >= 0.005:
                            d['weight'] = ((1/((H.nodes[u]['Water_Depth']))) + ((1/(H.nodes[v]['Water_Depth'])))) / 2
                        else:
                            d['weight'] = 10000000000000.0
                # Delete edges with a weight of 1000000000000 otherwise the path will take the path with the fewest nodes at the edge because the distances are larger and it has fewer nodes.
                selected_edges = [(u, v) for u, v, d in H.edges(data=True) if d['weight'] == 10000000000000.0]
                H.remove_edges_from(selected_edges)
                #H.edges(data=True)

            elif bedingung == 'ohneBedingungen' and talweg == "kuerzesterTalweg":
                print("Kürzester Talweg wird ohne Bedingungen berechnet")
                # 1/depth = x and then x/depth = y and again y/depth and so on to give the edges the right weighting.
                for u, v, d in H.edges(data=True):
                    if mindestwassertiefe_einzelfall != 0:
                        # d['weight'] = ((G.nodes[u]['Water_Depth'])+(G.nodes[v]['Water_Depth']))/2
                        if (H.nodes[u]['Water_Depth']) > 0.005 and (H.nodes[v]['Water_Depth']) > 0.005 and ((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth']) / 2) >= mindestwassertiefe_einzelfall:
                            d['weight'] = ((1/((H.nodes[u]['Water_Depth']))) + ((1/(H.nodes[v]['Water_Depth'])))) / 2
                        else:
                            d['weight'] = 10000000000000.0
                    elif mindestwassertiefe_einzelfall == 0:
                        # d['weight'] = ((G.nodes[u]['Water_Depth'])+(G.nodes[v]['Water_Depth']))/2
                        if (H.nodes[u]['Water_Depth']) > 0.005 and (H.nodes[v]['Water_Depth']) > 0.005 and ((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth']) / 2) >= 0.005:
                            d['weight'] = ((1/((H.nodes[u]['Water_Depth']))) + ((1/(H.nodes[v]['Water_Depth'])))) / 2
                        else:
                            d['weight'] = 10000000000000.0
                # Delete edges with a weight of 1000000000000 otherwise the path will take the path with the fewest nodes at the edge because the distances are larger and it has fewer nodes.
                selected_edges = [(u, v) for u, v, e in H.edges(data=True) if e['weight'] == 10000000000000.0]
                H.remove_edges_from(selected_edges)
                #H.edges(data=True)

            # nx.info(H,pos=emd)
            nodesremovelist = startthalweg.nodesremoveint
            H.remove_nodes_from(nodesremovelist)
            print("Folgende Nodes wurden aus dem Netzwerk gelöscht" + str(nodesremovelist))
            print("Folgende Informationen zum angepassten Netzwerk")
            print(nx.info(H))
            # Find shortest path
            allpaths1 = list(islice(nx.shortest_simple_paths(H, sourcenode, targetnode, weight="weight"), 1))
            allpaths.append(allpaths1)
            allpaths2 = np.asarray(allpaths1)
            thalwegtimestepslistnew = np.append(thalwegtimestepslistnew, i)

            # Adding waterdepth to nodesarray
            depthttimestep = Depth_TS[("TS" + str(i))]
            depthttimesteptranspose = depthttimestep.transpose()
            depthtoarray = np.insert(nodesarray, 4, values=depthttimesteptranspose, axis=1)

            # Adding flow velocity to nodesarray
            velocitytimestep = Velocity_TS[("TS" + str(i))]
            velocitytimesteptranspose = velocitytimestep.transpose()
            velocitytoarray = np.insert(depthtoarray, 5, values=velocitytimesteptranspose, axis=1)

            #Calculate the Variability of Water Depth and Flow Velocity
            indexd = np.argwhere(velocitytoarray[:,4] == 0.0)
            indexv = np.argwhere(velocitytoarray[:,5] == 0.0)
            velocitytoarraynewd = np.delete(velocitytoarray[:,4], indexd)
            velocitytoarraynewv = np.delete(velocitytoarray[:, 5], indexv)

            stdabd = np.std(velocitytoarraynewd)
            meand = np.mean(velocitytoarraynewd)
            vield = (1+(stdabd/meand))**2

            print("Die Vielfalt der Fliesstiefe beträgt " + str(vield) + " mit einem Mittelwert von " + str(meand) + " und einer Standardabweichung von " + str(stdabd))

            stdabv = np.std(velocitytoarraynewv)
            meanv = np.mean(velocitytoarraynewv)
            vielv = (1+(stdabv/meanv))**2
            print("Die Vielfalt der Fliessgeschwindigkeit beträgt " + str(vielv) + " mit einem Mittelwert von " + str(meanv) + " und einer Standardabweichung von " + str(stdabv))

            hmid = vield*vielv
            print("Der HMID beim TS " + str(i) + " beträgt " + str(hmid))

            # Adding wse to nodesarray
            wsetimestep = Wse_TS[("TS" + str(i))]
            wsetimesteptranspose = wsetimestep.transpose()
            wsetoarray = np.insert(velocitytoarray, 6, values=wsetimesteptranspose, axis=1)

            # Adding inflow to nodesarray
            inflow = inflowdat
            inflow1 = np.loadtxt(inflow)
            dischargein = inflow1[:, 1]
            dischargein1 = np.full((1, len(nodesarray)), dischargein[i])

            dischargeintoarray = np.insert(wsetoarray, 7, values=dischargein1, axis=1)

            # Adding outflow to nodesarray
            outflow = outflowdat
            outflow1 = np.loadtxt(outflow)
            discharge = outflow1[:, 1]
            discharge1 = np.full((1, len(nodesarray)), discharge[i])

            dischargetoarray = np.insert(dischargeintoarray, 8, values=discharge1, axis=1)

            # Leave only selected nodes in array test
            idx = allpaths1
            thalwegarray = dischargetoarray[np.in1d(dischargetoarray[:, 0], idx)]

            # Adding the drawing position
            test = thalwegarray.astype(int)
            test1 = test[:, 0]
            test2 = allpaths1[0]
            test3 = np.asarray(test2)

            position = np.array([])
            for element in test1:
                result = np.where(test3 == element)
                position = np.append(position, result)

            positiontranspose = position.transpose()

            thalwegarraypos = np.insert(thalwegarray, 9, values=positiontranspose, axis=1)

            # Sort the array with the thalweg nodes
            thalwegarraysort = np.sort(thalwegarraypos.view('f8,f8,f8,f8,f8,f8,f8,f8,f8,f8'), order=['f9'],
                                       axis=0).view(np.float)

            # Adding laenge and distance between nodes to nodesarray
            laenge = np.array([0.0])
            distancexy = np.array([0.0])

            index = len(thalwegarraysort) - 1
            e = 0
            while e < index:
                distance12 = math.sqrt(((thalwegarraysort[:, 1][e] - thalwegarraysort[:, 1][e + 1]) ** 2) + (
                            (thalwegarraysort[:, 2][e] - thalwegarraysort[:, 2][e + 1]) ** 2))
                distanceall = laenge[e] + distance12
                laenge = np.append(laenge, distanceall)
                distancexy = np.append(distancexy, distance12)

                e = e + 1

            thalwegarraylaenge = np.insert(thalwegarraysort, 10, values=laenge, axis=1)
            thalwegarraydis = np.insert(thalwegarraylaenge, 11, values=distancexy, axis=1)

            # Adding deltaz to nodesarray
            deltaz = np.array([0.0])
            index1 = len(thalwegarraysort) - 1
            o = 0
            while o < index1:
                distancedeltaz = math.sqrt(((thalwegarraysort[:, 1][o] - thalwegarraysort[:, 1][o + 1]) ** 2) + (
                            (thalwegarraysort[:, 2][o] - thalwegarraysort[:, 2][o + 1]) ** 2))
                deltazwse = abs(thalwegarraysort[:, 6][o] - thalwegarraysort[:, 6][o + 1]) / distancedeltaz
                deltaz = np.append(deltaz, deltazwse)
                o = o + 1

            thalwegarraydeltaz = np.insert(thalwegarraydis, 12, values=deltaz, axis=1)

            #Adding hmid to nodesarray
            hmid1 = np.full((1, len(thalwegarraydeltaz)), hmid)

            hmidtoarray = np.insert(thalwegarraydeltaz, 13, values=hmid1, axis=1)

            # Output to csv test
            float_formatter = lambda x: "%.8f" % x
            np.set_printoptions(formatter={'float_kind': float_formatter})
            allpaths3 = allpaths2.transpose()
            matrix = pd.DataFrame(data=hmidtoarray,
                                  columns=["Node", "X", "Y", "Z", "Wassertiefe", "Fliessgeschwindigkeit",
                                           "Wasserspiegel", "Inflow", "Outflow", "Position", "Gesamtdistanz",
                                           "Einzeldistanz", "Hindernishoehe", "HMID"])
            matrix.to_csv(outfiledirectory + '\\' + 'Talweg-TS' + str(i) + "-Abfluss" + str(
                round(thalwegarraydeltaz[:, 8][0], 2)) + '.csv', index=True, header=True, sep=',')

            # Sort the array with the water depth for the cummulative distribution of water depth
            thalwegarraydepthsort = np.sort(thalwegarraydeltaz.view('f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8'),
                                            order=['f4'], axis=0).view(np.float)

            # Adding the procentual distance of a waterdepth to the nodesarray
            prozentdepth = np.array([0.0])
            index2 = len(thalwegarraysort) - 1
            j = 0
            while j < index2:
                prozentdepthdis = ((thalwegarraydepthsort[:, 11][j]) / thalwegarraydeltaz[:, 10][index2]) * 100
                prozentdepthlaenge = prozentdepth[j] + prozentdepthdis
                prozentdepth = np.append(prozentdepth, prozentdepthlaenge)
                j = j + 1

            thalwegarrayprozentdepth = np.insert(thalwegarraydepthsort, 13, values=prozentdepth, axis=1)

            # Sort the array with the flow velocity depth for the cummulative distribution of flow velocity
            thalwegarrayvelocitysort = np.sort(thalwegarraydeltaz.view('f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8'),
                                               order=['f5'], axis=0).view(np.float)
            thalwegarrayvelocitysortreverse = thalwegarrayvelocitysort[::-1]

            # Adding the procentual distance of a flow velocity to the nodesarray
            prozentvelocity = np.array([0.0])
            index3 = len(thalwegarraysort) - 1
            r = 0
            while r < index3:
                prozentvelocitydis = ((thalwegarraydepthsort[:, 11][r]) / thalwegarraydeltaz[:, 10][index3]) * 100
                prozentvelocitylaenge = prozentvelocity[r] + prozentvelocitydis
                prozentvelocity = np.append(prozentvelocity, prozentvelocitylaenge)
                r = r + 1

            thalwegarrayprozentvelocity = np.insert(thalwegarrayvelocitysortreverse, 13, values=prozentvelocity, axis=1)

            # Sort the array with the deltaz for the cummulative distribution of deltaz
            thalwegarraydeltazsort = np.sort(thalwegarraydeltaz.view('f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8'),
                                             order=['f12'], axis=0).view(np.float)

            # Adding the procentual distance of a waterdepth to the nodesarray
            prozentdeltaz = np.array([0.0])
            index4 = len(thalwegarraysort) - 1
            s = 0
            while s < index4:
                prozentdeltazdis = ((thalwegarraydeltazsort[:, 11][s]) / thalwegarraydeltaz[:, 10][index4]) * 100
                prozentdeltazlaenge = prozentdeltaz[s] + prozentdeltazdis
                prozentdeltaz = np.append(prozentdeltaz, prozentdeltazlaenge)
                s = s + 1

            thalwegarrayprozentdeltaz = np.insert(thalwegarraydeltazsort, 13, values=prozentdeltaz, axis=1)

            # Frequency and length of the undershoots of the general and minimum water depths for shallows along the deepest channel
            deptharray = np.array([])
            deptharray1 = thalwegarraydeltaz[:, 4]
            index5 = len(thalwegarraydeltaz[:, 4]) - 1
            p = 0
            while p <= index5:
                deptharray = np.append(deptharray, deptharray1[p])
                p = p + 1

            lengthkmu = []
            lengthkmg = []
            meankmg = []
            meankmu = []

            # Different conditions to make different statements depending on the input
            if mindestwassertiefe_untiefen == 0 and mindestwassertiefe_generell == 0:
                lengthkmu = thalwegarraydeltaz[:, 11]
                meankmu = deptharray

                frequencykmu = len(lengthkmu)
            elif mindestwassertiefe_untiefen > 0.0 and mindestwassertiefe_generell == 0:
                m = np.where(deptharray < mindestwassertiefe_untiefen)  # mindestwassertiefe_untiefen
                m = m[0]
                kmu = m
                if kmu.size == 0:
                    print("no Points under the Mindestwassertiefe Untiefen")
                else:
                    m = np.where(deptharray < mindestwassertiefe_untiefen)  # mindestwassertiefe_untiefen
                    m = m[0]
                    kmu = m
                    kmui = [(x + 1) for x, y in zip(m, m[1:]) if y - x != 1]
                    kmuiconsecutive = []
                    prev = 0
                    for index in kmui:
                        new_list = [x for x in kmu[prev:] if x < index]
                        kmuiconsecutive.append(new_list)
                        prev += len(new_list)
                    kmuiconsecutive.append([x for x in kmu[prev:]])

                    lengthkmu = []
                    for listu in kmuiconsecutive:
                        lengthkmu.append(sum(thalwegarraydeltaz[:, 11][listu]))

                    meankmu = []
                    for listu in kmuiconsecutive:
                        meankmu.append(
                            (sum(thalwegarraydeltaz[:, 4][listu])) / float(len(thalwegarraydeltaz[:, 4][listu])))
            elif mindestwassertiefe_untiefen == 0 and mindestwassertiefe_generell > 0.0:
                n = np.where(deptharray < mindestwassertiefe_generell)  # mindestwassertiefe_generell
                n = n[0]
                kmg = n
                if kmg.size == 0:
                    print("no Points under the Mindestwassertiefe generell")
                else:
                    n = np.where(deptharray < mindestwassertiefe_generell)  # mindestwassertiefe_generell
                    n = n[0]
                    kmg = n
                    kmgi = [(x + 1) for x, y in zip(n, n[1:]) if y - x != 1]
                    kmgiconsecutive = []
                    prev = 0
                    for index in kmgi:
                        new_list = [x for x in kmg[prev:] if x < index]
                        kmgiconsecutive.append(new_list)
                        prev += len(new_list)
                    kmgiconsecutive.append([x for x in kmg[prev:]])

                    lengthkmg = []
                    for listg in kmgiconsecutive:
                        lengthkmg.append(sum(thalwegarraydeltaz[:, 11][listg]))

                    meankmg = []
                    for listg in kmgiconsecutive:
                        meankmg.append(
                            (sum(thalwegarraydeltaz[:, 4][listg])) / float(len(thalwegarraydeltaz[:, 4][listg])))
            elif mindestwassertiefe_generell > 0.0 and mindestwassertiefe_untiefen > 0.0 and mindestwassertiefe_generell != mindestwassertiefe_untiefen:
                m = np.where(deptharray < mindestwassertiefe_untiefen)  # mindestwassertiefe_untiefen
                m = m[0]
                kmu = m
                n = np.where(deptharray < mindestwassertiefe_generell)  # mindestwassertiefe_generell
                n = n[0]
                kmg = n
                if kmu.size == 0 and kmg.size > 0:
                    n = np.where(deptharray < mindestwassertiefe_generell)  # mindestwassertiefe_generell
                    n = n[0]
                    kmg = n
                    kmgi = [(x + 1) for x, y in zip(n, n[1:]) if y - x != 1]
                    kmgiconsecutive = []
                    prev = 0
                    for index in kmgi:
                        new_list = [x for x in kmg[prev:] if x < index]
                        kmgiconsecutive.append(new_list)
                        prev += len(new_list)
                    kmgiconsecutive.append([x for x in kmg[prev:]])

                    lengthkmg = []
                    for listg in kmgiconsecutive:
                        lengthkmg.append(sum(thalwegarraydeltaz[:, 11][listg]))

                    meankmg = []
                    for listg in kmgiconsecutive:
                        meankmg.append(
                            (sum(thalwegarraydeltaz[:, 4][listg])) / float(len(thalwegarraydeltaz[:, 4][listg])))
                elif kmg.size == 0 and kmu.size > 0:
                    m = np.where(deptharray < mindestwassertiefe_untiefen)  # mindestwassertiefe_untiefen
                    m = m[0]
                    kmu = m
                    kmui = [(x + 1) for x, y in zip(m, m[1:]) if y - x != 1]
                    kmuiconsecutive = []
                    prev = 0
                    for index in kmui:
                        new_list = [x for x in kmu[prev:] if x < index]
                        kmuiconsecutive.append(new_list)
                        prev += len(new_list)
                    kmuiconsecutive.append([x for x in kmu[prev:]])

                    lengthkmu = []
                    for listu in kmuiconsecutive:
                        lengthkmu.append(sum(thalwegarraydeltaz[:, 11][listu]))

                    meankmu = []
                    for listu in kmuiconsecutive:
                        meankmu.append(
                            (sum(thalwegarraydeltaz[:, 4][listu])) / float(len(thalwegarraydeltaz[:, 4][listu])))
                elif kmg.size > 0 and kmu.size > 0:
                    if mindestwassertiefe_generell > mindestwassertiefe_untiefen:
                        m = np.where(deptharray < mindestwassertiefe_untiefen)  # mindestwassertiefe_untiefen
                        m = m[0]
                        kmu = m
                        kmui = [(x + 1) for x, y in zip(m, m[1:]) if y - x != 1]
                        kmuiconsecutive = []
                        prev = 0
                        for index in kmui:
                            new_list = [x for x in kmu[prev:] if x < index]
                            kmuiconsecutive.append(new_list)
                            prev += len(new_list)
                        kmuiconsecutive.append([x for x in kmu[prev:]])

                        lengthkmu = []
                        for listu in kmuiconsecutive:
                            lengthkmu.append(sum(thalwegarraydeltaz[:, 11][listu]))

                        meankmu = []
                        for listu in kmuiconsecutive:
                            meankmu.append(
                                (sum(thalwegarraydeltaz[:, 4][listu])) / float(len(thalwegarraydeltaz[:, 4][listu])))
                        n = np.where(deptharray < mindestwassertiefe_generell)  # mindestwassertiefe_generell
                        n = n[0]
                        for indi in m:
                            delindex = np.where(n == indi)
                            n = np.delete(n, delindex[0][0])
                        kmg = n
                        kmgi = [(x + 1) for x, y in zip(n, n[1:]) if y - x != 1]
                        kmgiconsecutive = []
                        prev = 0
                        for index in kmgi:
                            new_list = [x for x in kmg[prev:] if x < index]
                            kmgiconsecutive.append(new_list)
                            prev += len(new_list)
                        kmgiconsecutive.append([x for x in kmg[prev:]])

                        lengthkmg = []
                        for listg in kmgiconsecutive:
                            lengthkmg.append(sum(thalwegarraydeltaz[:, 11][listg]))

                        meankmg = []
                        for listg in kmgiconsecutive:
                            meankmg.append(
                                (sum(thalwegarraydeltaz[:, 4][listg])) / float(len(thalwegarraydeltaz[:, 4][listg])))
                    else:
                        n = np.where(deptharray < mindestwassertiefe_generell)  # mindestwassertiefe_generell
                        n = n[0]
                        kmg = n
                        kmgi = [(x + 1) for x, y in zip(n, n[1:]) if y - x != 1]
                        kmgiconsecutive = []
                        prev = 0
                        for index in kmgi:
                            new_list = [x for x in kmg[prev:] if x < index]
                            kmgiconsecutive.append(new_list)
                            prev += len(new_list)
                        kmgiconsecutive.append([x for x in kmg[prev:]])

                        lengthkmg = []
                        for listg in kmgiconsecutive:
                            lengthkmg.append(sum(thalwegarraydeltaz[:, 11][listg]))

                        meankmg = []
                        for listg in kmgiconsecutive:
                            meankmg.append(
                                (sum(thalwegarraydeltaz[:, 4][listg])) / float(len(thalwegarraydeltaz[:, 4][listg])))

                        m = np.where(deptharray < mindestwassertiefe_untiefen)  # mindestwassertiefe_untiefen
                        m = m[0]
                        for indi in n:
                            delindex = np.where(m == indi)
                            m = np.delete(m, delindex[0][0])
                        kmu = m
                        kmui = [(x + 1) for x, y in zip(m, m[1:]) if y - x != 1]
                        kmuiconsecutive = []
                        prev = 0
                        for index in kmui:
                            new_list = [x for x in kmu[prev:] if x < index]
                            kmuiconsecutive.append(new_list)
                            prev += len(new_list)
                        kmuiconsecutive.append([x for x in kmu[prev:]])

                        lengthkmu = []
                        for listu in kmuiconsecutive:
                            lengthkmu.append(sum(thalwegarraydeltaz[:, 11][listu]))

                        meankmu = []
                        for listu in kmuiconsecutive:
                            meankmu.append(
                                (sum(thalwegarraydeltaz[:, 4][listu])) / float(len(thalwegarraydeltaz[:, 4][listu])))

            elif mindestwassertiefe_generell > 0.0 and mindestwassertiefe_untiefen > 0.0 and mindestwassertiefe_generell == mindestwassertiefe_untiefen:
                m = np.where(deptharray < mindestwassertiefe_untiefen)  # mindestwassertiefe_untiefen
                m = m[0]
                kmu = m
                n = np.where(deptharray < mindestwassertiefe_generell)  # mindestwassertiefe_generell
                n = n[0]
                kmg = n
                if kmu.size == 0 and kmg.size > 0:
                    n = np.where(deptharray < mindestwassertiefe_generell)  # mindestwassertiefe_generell
                    n = n[0]
                    kmg = n
                    kmgi = [(x + 1) for x, y in zip(n, n[1:]) if y - x != 1]
                    kmgiconsecutive = []
                    prev = 0
                    for index in kmgi:
                        new_list = [x for x in kmg[prev:] if x < index]
                        kmgiconsecutive.append(new_list)
                        prev += len(new_list)
                    kmgiconsecutive.append([x for x in kmg[prev:]])

                    lengthkmg = []
                    for listg in kmgiconsecutive:
                        lengthkmg.append(sum(thalwegarraydeltaz[:, 11][listg]))

                    meankmg = []
                    for listg in kmgiconsecutive:
                        meankmg.append(
                            (sum(thalwegarraydeltaz[:, 4][listg])) / float(len(thalwegarraydeltaz[:, 4][listg])))
                elif kmg.size == 0 and kmu.size > 0:
                    m = np.where(deptharray < mindestwassertiefe_untiefen)  # mindestwassertiefe_untiefen
                    m = m[0]
                    kmu = m
                    kmui = [(x + 1) for x, y in zip(m, m[1:]) if y - x != 1]
                    kmuiconsecutive = []
                    prev = 0
                    for index in kmui:
                        new_list = [x for x in kmu[prev:] if x < index]
                        kmuiconsecutive.append(new_list)
                        prev += len(new_list)
                    kmuiconsecutive.append([x for x in kmu[prev:]])

                    lengthkmu = []
                    for listu in kmuiconsecutive:
                        lengthkmu.append(sum(thalwegarraydeltaz[:, 11][listu]))

                    meankmu = []
                    for listu in kmuiconsecutive:
                        meankmu.append(
                            (sum(thalwegarraydeltaz[:, 4][listu])) / float(len(thalwegarraydeltaz[:, 4][listu])))
                elif kmg.size > 0 and kmu.size > 0:
                    m = np.where(deptharray < mindestwassertiefe_untiefen)  # mindestwassertiefe_untiefen
                    m = m[0]
                    kmu = m
                    kmui = [(x + 1) for x, y in zip(m, m[1:]) if y - x != 1]
                    kmuiconsecutive = []
                    prev = 0
                    for index in kmui:
                        new_list = [x for x in kmu[prev:] if x < index]
                        kmuiconsecutive.append(new_list)
                        prev += len(new_list)
                    kmuiconsecutive.append([x for x in kmu[prev:]])

                    lengthkmu = []
                    for listu in kmuiconsecutive:
                        lengthkmu.append(sum(thalwegarraydeltaz[:, 11][listu]))

                    meankmu = []
                    for listu in kmuiconsecutive:
                        meankmu.append(
                            (sum(thalwegarraydeltaz[:, 4][listu])) / float(len(thalwegarraydeltaz[:, 4][listu])))
                    n = np.where(deptharray < mindestwassertiefe_generell)  # mindestwassertiefe_generell
                    n = n[0]
                    kmg = n
                    kmgi = [(x + 1) for x, y in zip(n, n[1:]) if y - x != 1]
                    kmgiconsecutive = []
                    prev = 0
                    for index in kmgi:
                        new_list = [x for x in kmg[prev:] if x < index]
                        kmgiconsecutive.append(new_list)
                        prev += len(new_list)
                    kmgiconsecutive.append([x for x in kmg[prev:]])

                    lengthkmg = []
                    for listg in kmgiconsecutive:
                        lengthkmg.append(sum(thalwegarraydeltaz[:, 11][listg]))

                    meankmg = []
                    for listg in kmgiconsecutive:
                        meankmg.append(
                            (sum(thalwegarraydeltaz[:, 4][listg])) / float(len(thalwegarraydeltaz[:, 4][listg])))

            # Discharge mindestwassertiefe_generell und mindestwassertiefe_untiefen
            lengthkmgl = thalwegarraydeltaz[:, 8][0] * 1
            lengthkmul = thalwegarraydeltaz[:, 8][0] * 1

            ###Plots###
            #WSE, riverbed and fall height
            dis = thalwegarraydeltaz[:,10]
            z = thalwegarraydeltaz[:,3]
            wse = thalwegarraydeltaz[:,6]
            deltaz = thalwegarraydeltaz[:,12]

            plt.figure()
            color = 'tab:grey'
            color1 = 'tab:blue'
            plt.subplot(211)
            plt.plot(dis, z, label='Sohle', color=color, linewidth=0.8)
            plt.plot(dis, wse, label='Wasserspiegel', color=color1, linewidth=0.8)
            plt.xlabel('Distanz (m)', color="black",fontsize=6)
            plt.ylabel('Höhe (m.ü.M)', color="black",fontsize=6)
            plt.title('Wasserspiegel und Sohlenlage entlang der tiefsten Rinne' + "\n TS: " + str(i) + " und Abfluss: " + str(round(thalwegarraydeltaz[:, 8][0],1)) + r" $m^3/s$", fontsize=10)

            plt.legend(loc='best', frameon=False,fontsize=6)
            plt.tick_params(axis='y', labelcolor="black")

            color = 'tab:brown'
            plt.subplot(212)
            plt.bar(dis, deltaz, label='Hindernishöhe', color=color)
            plt.title('Hindernishöhen entlang der tiefsten Rinne', fontsize=10)
            plt.xlabel('Distanz (m)', color="black", fontsize=6)
            plt.ylabel('Höhe (m)', color="black",fontsize=6)  # we already handled the x-label with ax1
            if maximalabsturzhoehe == 0:
                plt.legend(loc='best', frameon=False, fontsize=6)
                plt.tick_params(axis='y', labelcolor="black")

                # fig.tight_layout()  # otherwise the right y-label is slightly clipped
                plt.grid(linestyle='-', linewidth=0.5)
                plt.tight_layout()
                plt.show()
                plt.savefig(
                    outfiledirectory + '\\' + 'Wasserspiegel-Sohlenlage und Hindernishöhen entlang der tiefsten Rinne' + " TS" + str(
                        i) + '.' + plotformat)
            else:
                plt.axhline(y=maximalabsturzhoehe, color='r', linestyle='-', label='maximale Sprunghöhe', linewidth=0.8)
                plt.legend(loc='best', frameon=False,fontsize=6)
                plt.tick_params(axis='y', labelcolor="black")

                #fig.tight_layout()  # otherwise the right y-label is slightly clipped
                plt.grid(linestyle='-', linewidth=0.5)
                plt.tight_layout()
                plt.show()
                plt.savefig(outfiledirectory + '\\' + 'Wasserspiegel-Sohlenlage und Hindernishöhen entlang der tiefsten Rinne' + " TS" + str(i) + '.' + plotformat)

            #Water Depth
            dis = thalwegarraydeltaz[:,10]
            depth = thalwegarraydeltaz[:,4]

            fig, ax1 = plt.subplots()

            color = 'tab:blue'
            plt.title('Wassertiefen entlang der tiefsten Rinne' + "\n TS: " + str(i) + " und Abfluss: " + str(round(thalwegarraydeltaz[:, 8][0],1)) + r" $m^3/s$", fontsize=10)
            ax1.set_xlabel('Distanz (m)', color="black", fontsize=6)
            ax1.set_ylabel('Wassertiefe (m)', color="black", fontsize=6)
            ax1.plot(dis, depth, label='Wassertiefe', color=color, linewidth=0.8)
            if mindestwassertiefe_generell == 0 and mindestwassertiefe_untiefen == 0:
                ax1.legend(loc='best', frameon=False, fontsize=6)
                ax1.tick_params(axis='y', labelcolor="black")

                fig.tight_layout()  # otherwise the right y-label is slightly clipped
                ax1.set_axisbelow(True)
                ax1.grid(linestyle='-', linewidth=0.5)
                plt.show()
                plt.savefig(
                    outfiledirectory + '\\' + 'Wassertiefen entlang der tiefsten Rinne' + " TS" + str(i) + '.' + plotformat)
            elif mindestwassertiefe_generell > 0 and mindestwassertiefe_untiefen == 0:
                ax1.axhline(y=mindestwassertiefe_generell, color='r', linestyle='-', label='Generelle Mindestwassertiefe', linewidth=0.8)
                ax1.legend(loc='best', frameon=False,fontsize=6)
                ax1.tick_params(axis='y', labelcolor="black")

                fig.tight_layout()  # otherwise the right y-label is slightly clipped
                ax1.set_axisbelow(True)
                ax1.grid(linestyle='-', linewidth=0.5)
                plt.show()
                plt.savefig(outfiledirectory + '\\' + 'Wassertiefen entlang der tiefsten Rinne' + " TS" + str(i) + '.' + plotformat)
            elif mindestwassertiefe_generell == 0 and mindestwassertiefe_untiefen > 0:
                ax1.axhline(y=mindestwassertiefe_untiefen, color = 'r', linestyle = 'dashed', label = 'Mindestwassertiefe für Untiefen', linewidth=0.8)
                ax1.legend(loc='best', frameon=False,fontsize=6)
                ax1.tick_params(axis='y', labelcolor="black")

                fig.tight_layout()  # otherwise the right y-label is slightly clipped
                ax1.set_axisbelow(True)
                ax1.grid(linestyle='-', linewidth=0.5)
                plt.show()
                plt.savefig(outfiledirectory + '\\' + 'Wassertiefen entlang der tiefsten Rinne' + " TS" + str(i) + '.' + plotformat)
            else:
                ax1.axhline(y=mindestwassertiefe_generell, color='r', linestyle='-', label='Generelle Mindestwassertiefe', linewidth=0.8)
                ax1.axhline(y=mindestwassertiefe_untiefen, color = 'r', linestyle = 'dashed', label = 'Mindestwassertiefe für Untiefen', linewidth=0.8)
                ax1.legend(loc='best', frameon=False,fontsize=6)
                ax1.tick_params(axis='y', labelcolor="black")

                fig.tight_layout()  # otherwise the right y-label is slightly clipped
                ax1.set_axisbelow(True)
                ax1.grid(linestyle='-', linewidth=0.5)
                plt.show()
                plt.savefig(outfiledirectory + '\\' + 'Wassertiefen entlang der tiefsten Rinne' + " TS" + str(i) + '.' + plotformat)

            #Flow Velocity
            dis = thalwegarraydeltaz[:, 10]
            velocity = thalwegarraydeltaz[:, 5]

            fig, ax1 = plt.subplots()

            color = 'tab:cyan'
            plt.title('Fliessgeschwindigkeiten entlang der tiefsten Rinne'+ "\n TS: " + str(i) + " und Abfluss: " + str(round(thalwegarraydeltaz[:, 8][0],1)) + r" $m^3/s$",fontsize = 10)
            ax1.set_xlabel('Distanz (m)', color="black", fontsize=6)
            ax1.set_ylabel('Geschwindigkeit (m/s)', color="black", fontsize=6)
            ax1.plot(dis, velocity, label='Geschwindigkeit', color=color, linewidth=0.8)
            if maximalgeschwindigkeit == 0:
                ax1.legend(loc='best', frameon=False, fontsize=6)
                ax1.tick_params(axis='y', labelcolor="black")

                fig.tight_layout()  # otherwise the right y-label is slightly clipped
                ax1.set_axisbelow(True)
                ax1.grid(linestyle='-', linewidth=0.5)
                plt.show()
                plt.savefig(
                    outfiledirectory + '\\' + 'Fliessgeschwindigkeiten entlang der tiefsten Rinne' + " TS" + str(
                        i) + '.' + plotformat)
            else:
                ax1.axhline(y=maximalgeschwindigkeit, color='r', linestyle='-', label='Maximalgeschwindigkeit', linewidth=0.8)
                ax1.legend(loc='best', frameon=False,fontsize=6)
                ax1.tick_params(axis='y', labelcolor="black")

                fig.tight_layout()  # otherwise the right y-label is slightly clipped
                ax1.set_axisbelow(True)
                ax1.grid(linestyle='-', linewidth=0.5)
                plt.show()
                plt.savefig(outfiledirectory + '\\' + 'Fliessgeschwindigkeiten entlang der tiefsten Rinne' + " TS" + str(i) + '.' + plotformat)

            #Cummulative distribution of water depth
            dis = thalwegarrayprozentdepth[:,13]
            depth = thalwegarrayprozentdepth[:,4]

            fig, ax1 = plt.subplots()

            color = 'tab:blue'
            plt.title('Kumulative Verteilungskurve der Wassertiefen entlang des Talwegs'+ "\n TS: " + str(i) + " und Abfluss: " + str(round(thalwegarraydeltaz[:, 8][0],1)) + r" $m^3/s$",fontsize = 10)
            ax1.set_xlabel('Wassertiefe (m)', color="black", fontsize=6)
            ax1.set_ylabel('Anteil (%)', color="black", fontsize=6)
            ax1.plot(depth, dis, label= 'Wassertiefe', color=color, linewidth=0.8)
            if mindestwassertiefe_generell == 0 and mindestwassertiefe_untiefen == 0:
                ax1.legend(loc='best', frameon=False, fontsize=6)
                ax1.tick_params(axis='y', labelcolor="black")
                fig.tight_layout()  # otherwise the right y-label is slightly clipped
                ax1.set_axisbelow(True)
                ax1.grid(linestyle='-', linewidth=0.5)
                plt.show()
                plt.savefig(
                    outfiledirectory + '\\' + 'Kumulative Verteilungskurve der Wassertiefen entlang des Talwegs' + " TS" + str(
                        i) + '.' + plotformat)
            elif mindestwassertiefe_generell > 0 and mindestwassertiefe_untiefen == 0:
                ax1.axvline(x=mindestwassertiefe_generell, color='r', linestyle='-', label='Generelle Mindestwassertiefe', linewidth=0.8)
                ax1.legend(loc='best', frameon=False,fontsize=6)
                ax1.tick_params(axis='y', labelcolor="black")
                fig.tight_layout()  # otherwise the right y-label is slightly clipped
                ax1.set_axisbelow(True)
                ax1.grid(linestyle='-', linewidth=0.5)
                plt.show()
                plt.savefig(outfiledirectory + '\\' + 'Kumulative Verteilungskurve der Wassertiefen entlang des Talwegs' + " TS" + str(i) + '.' + plotformat)

            elif mindestwassertiefe_generell == 0 and mindestwassertiefe_untiefen > 0:
                ax1.axvline(x=mindestwassertiefe_untiefen, color='r', linestyle='dashed', label='Mindestwassertiefe für Untiefen', linewidth=0.8)
                ax1.legend(loc='best', frameon=False,fontsize=6)
                ax1.tick_params(axis='y', labelcolor="black")
                fig.tight_layout()  # otherwise the right y-label is slightly clipped
                ax1.set_axisbelow(True)
                ax1.grid(linestyle='-', linewidth=0.5)
                plt.show()
                plt.savefig(outfiledirectory + '\\' + 'Kumulative Verteilungskurve der Wassertiefen entlang des Talwegs' + " TS" + str(i) + '.' + plotformat)

            else:
                ax1.axvline(x=mindestwassertiefe_generell, color='r', linestyle='-', label='Generelle Mindestwassertiefe', linewidth=0.8)
                ax1.axvline(x=mindestwassertiefe_untiefen, color='r', linestyle='dashed', label='Mindestwassertiefe für Untiefen', linewidth=0.8)
                ax1.legend(loc='best', frameon=False,fontsize=6)
                ax1.tick_params(axis='y', labelcolor="black")
                fig.tight_layout()  # otherwise the right y-label is slightly clipped
                ax1.set_axisbelow(True)
                ax1.grid(linestyle='-', linewidth=0.5)
                plt.show()
                plt.savefig(outfiledirectory + '\\' + 'Kumulative Verteilungskurve der Wassertiefen entlang des Talwegs' + " TS" + str(i) + '.' + plotformat)

            #Cummulative distribution of flow velocity
            dis = thalwegarrayprozentvelocity[:,13]
            velocity = thalwegarrayprozentvelocity[:,5]

            fig, ax1 = plt.subplots()

            color = 'tab:cyan'
            plt.title('Kumulative Verteilungskurve der Fliessgeschwindigkeiten entlang des Talwegs'+ "\n TS: " + str(i) + " und Abfluss: " + str(round(thalwegarraydeltaz[:, 8][0],1)) + r" $m^3/s$",fontsize = 10)
            ax1.set_xlabel('Fliessgeschwindigkeit (m/s)', color="black", fontsize=6)
            ax1.set_ylabel('Anteil (%)', color="black", fontsize=6)
            ax1.plot(velocity, dis, label= 'Fliessgeschwindigkeit', color=color, linewidth=0.8)
            if maximalgeschwindigkeit == 0:
                ax1.legend(loc='best', frameon=False, fontsize=6)
                ax1.tick_params(axis='y', labelcolor="black")

                fig.tight_layout()  # otherwise the right y-label is slightly clipped
                ax1.set_axisbelow(True)
                ax1.grid(linestyle='-', linewidth=0.5)
                plt.show()
                plt.savefig(
                    outfiledirectory + '\\' + 'Kumulative Verteilungskurve der Fliessgeschindigkeiten entlang des Talwegs' + " TS" + str(
                        i) + '.' + plotformat)
            else:
                ax1.axvline(x=maximalgeschwindigkeit, color='r', linestyle='-', label='Maximalgeschwindigkeit', linewidth=0.8)
                ax1.legend(loc='best', frameon=False,fontsize=6)
                ax1.tick_params(axis='y', labelcolor="black")

                fig.tight_layout()  # otherwise the right y-label is slightly clipped
                ax1.set_axisbelow(True)
                ax1.grid(linestyle='-', linewidth=0.5)
                plt.show()
                plt.savefig(outfiledirectory + '\\' + 'Kumulative Verteilungskurve der Fliessgeschindigkeiten entlang des Talwegs' + " TS" + str(i) + '.' + plotformat)

            #Cummulative distribution of fall height
            dis = thalwegarrayprozentdeltaz[:,13]
            deltaz = thalwegarrayprozentdeltaz[:,12]

            fig, ax1 = plt.subplots()

            color = 'tab:brown'
            plt.title('Kumulative Verteilungskurve der Hindernishöhen entlang des Talwegs'+ "\n TS: " + str(i) + " und Abfluss: " + str(round(thalwegarraydeltaz[:, 8][0],1)) + r" $m^3/s$",fontsize = 10)
            ax1.set_xlabel('Hindernishöhe (m)', color="black",fontsize=6)
            ax1.set_ylabel('Anteil (%)', color="black", fontsize=6)
            ax1.plot(deltaz, dis, label= 'Hindernishöhe', color=color, linewidth=0.8)
            if maximalabsturzhoehe == 0:
                ax1.legend(loc='best', frameon=False, fontsize=6)
                ax1.tick_params(axis='y', labelcolor="black")

                fig.tight_layout()  # otherwise the right y-label is slightly clipped
                ax1.set_axisbelow(True)
                ax1.grid(linestyle='-', linewidth=0.5)
                plt.show()
                plt.savefig(
                    outfiledirectory + '\\' + 'Kumulative Verteilungskurve der Hindernishöhen entlang des Talwegs' + " TS" + str(
                        i) + '.' + plotformat)
            else:
                ax1.axvline(x=maximalabsturzhoehe, color='r', linestyle='-', label='Maximale Sprunghöhe', linewidth=0.8)
                ax1.legend(loc='best', frameon=False, fontsize=6)
                ax1.tick_params(axis='y', labelcolor="black")

                fig.tight_layout()  # otherwise the right y-label is slightly clipped
                ax1.set_axisbelow(True)
                ax1.grid(linestyle='-', linewidth=0.5)
                plt.show()
                plt.savefig(outfiledirectory+ '\\'+ 'Kumulative Verteilungskurve der Hindernishöhen entlang des Talwegs' + " TS" + str(i)+ '.' + plotformat)

            # Frequency and length of the underruns of the general and minimum water depths for shallows along the deepest channel
            # Percentage of the length of the minimum water depth not reached.
            meankmgplot = meankmg
            meankmuplot = meankmu
            lengthkmgplot = lengthkmg
            lengthkmuplot = lengthkmu

            fig, ax1 = plt.subplots()

            color = 'tab:blue'
            color1 = 'tab:red'
            plt.title('Häufigkeit und Länge der Unterschreitungen der generellen und \n der Mindestwassertiefe für Untiefen entlang der tiefsten Rinne' + "\n TS: " + str(i) + " und Abfluss: " + str(round(thalwegarraydeltaz[:, 8][0], 1)) + r" $m^3/s$", fontsize=10)
            ax1.set_xlabel('Mittlere Tiefe der Unterschreitung (m)', color="black", fontsize=6)
            ax1.set_ylabel('Länge (m)', color="black", fontsize=6)
            if maximallpassuntiefen == 0 and maximallpasseinzelfall == 0:
                print("Keine maximale Länge für passierbare Untiefen und Einzelfälle definiert")
                if mindestwassertiefe_generell == 0 and mindestwassertiefe_untiefen == 0:
                    ax1.set_axisbelow(True)
                    ax1.grid(linestyle='-', linewidth=0.5)
                    ax1.scatter(meankmuplot, lengthkmuplot, s=4, label="Points", color=color1)
                    ax1.legend(loc='best', frameon=False, fontsize=6)
                    # plt.figtext(0.80, 0.08, 'n = ' + str(len(lengthkmgplot)) + "/" + str(len(lengthkmuplot)), color='black', weight='roman',size='x-small')
                    ax1.text(0.0, 1.04, 'n = ' + str(len(lengthkmgplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmgplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm', color='blue', horizontalalignment='left',
                             verticalalignment='top', transform=ax1.transAxes, fontsize=6)
                    ax1.text(1.0, 1.04,
                             'n = ' + str(len(lengthkmuplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmuplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm',
                             color='red', horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes,
                             fontsize=6)
                    ax1.tick_params(axis='y', labelcolor="black")

                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.show()
                    plt.savefig(
                        outfiledirectory + '\\' + 'Häufigkeit und Länge der Tiefenunterschreitungen entlang der tiefsten Rinne' + " TS" + str(
                            i) + '.' + plotformat)
                elif mindestwassertiefe_generell == 0 and mindestwassertiefe_untiefen > 0.0:
                    ax1.set_axisbelow(True)
                    ax1.grid(linestyle='-', linewidth=0.5)
                    ax1.scatter(meankmuplot, lengthkmuplot, s=4,
                                label="Mindestwassertiefe für Unftiefen <" + str(mindestwassertiefe_untiefen) + 'm',
                                color=color1)
                    ax1.legend(loc='best', frameon=False, fontsize=6)
                    # plt.figtext(0.80, 0.08, 'n = ' + str(len(lengthkmgplot)) + "/" + str(len(lengthkmuplot)), color='black', weight='roman',size='x-small')
                    ax1.text(1.0, 1.04,
                             'n = ' + str(len(lengthkmuplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmuplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm',
                             color='red', horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes,
                             fontsize=6)
                    ax1.tick_params(axis='y', labelcolor="black")

                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.show()
                    plt.savefig(
                        outfiledirectory + '\\' + 'Häufigkeit und Länge der Tiefenunterschreitungen entlang der tiefsten Rinne' + " TS" + str(
                            i) + '.' + plotformat)
                elif mindestwassertiefe_generell > 0.0 and mindestwassertiefe_untiefen == 0:
                    ax1.set_axisbelow(True)
                    ax1.grid(linestyle='-', linewidth=0.5)
                    ax1.scatter(meankmgplot, lengthkmgplot, s=4,
                                label="Generelle Mindestwassertiefe <" + str(mindestwassertiefe_generell) + 'm',
                                color=color)
                    ax1.legend(loc='best', frameon=False, fontsize=6)
                    # plt.figtext(0.80, 0.08, 'n = ' + str(len(lengthkmgplot)) + "/" + str(len(lengthkmuplot)), color='black', weight='roman',size='x-small')
                    ax1.text(0.0, 1.04, 'n = ' + str(len(lengthkmgplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmgplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm', color='blue', horizontalalignment='left',
                             verticalalignment='top', transform=ax1.transAxes, fontsize=6)
                    ax1.tick_params(axis='y', labelcolor="black")

                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.show()
                    plt.savefig(
                        outfiledirectory + '\\' + 'Häufigkeit und Länge der Tiefenunterschreitungen entlang der tiefsten Rinne' + " TS" + str(
                            i) + '.' + plotformat)
                else:
                    ax1.set_axisbelow(True)
                    ax1.grid(linestyle='-', linewidth=0.5)
                    ax1.scatter(meankmgplot, lengthkmgplot, s=4,
                                label="Generelle Mindestwassertiefe <" + str(mindestwassertiefe_generell) + 'm',
                                color=color)
                    ax1.scatter(meankmuplot, lengthkmuplot, s=4,
                                label="Mindestwassertiefe für Unftiefen <" + str(mindestwassertiefe_untiefen) + 'm',
                                color=color1)
                    ax1.legend(loc='best', frameon=False, fontsize=6)
                    # plt.figtext(0.80, 0.08, 'n = ' + str(len(lengthkmgplot)) + "/" + str(len(lengthkmuplot)), color='black', weight='roman',size='x-small')
                    ax1.text(0.0, 1.04, 'n = ' + str(len(lengthkmgplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmgplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm', color='blue', horizontalalignment='left',
                             verticalalignment='top', transform=ax1.transAxes, fontsize=6)
                    ax1.text(1.0, 1.04,
                             'n = ' + str(len(lengthkmuplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmuplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm',
                             color='red', horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes,
                             fontsize=6)
                    ax1.tick_params(axis='y', labelcolor="black")

                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.show()
                    plt.savefig(
                        outfiledirectory + '\\' + 'Häufigkeit und Länge der Tiefenunterschreitungen entlang der tiefsten Rinne' + " TS" + str(
                            i) + '.' + plotformat)

            elif maximallpassuntiefen == 0 and maximallpasseinzelfall > 0:
                if mindestwassertiefe_generell == 0 and mindestwassertiefe_untiefen == 0:
                    ax1.set_axisbelow(True)
                    ax1.grid(linestyle='-', linewidth=0.5)
                    ax1.scatter(meankmuplot, lengthkmuplot, s=4, label="Points", color=color1, zorder=2)
                    ax1.axhline(y=maximallpasseinzelfall, color='r', linestyle='dashed',
                                label='maximale Länge passierbarer Einzelfälle',
                                linewidth=0.8, zorder=1)
                    ax1.legend(loc='best', frameon=False, fontsize=6)
                    # plt.figtext(0.80, 0.08, 'n = ' + str(len(lengthkmgplot)) + "/" + str(len(lengthkmuplot)), color='black', weight='roman',size='x-small')
                    ax1.text(0.0, 1.04, 'n = ' + str(len(lengthkmgplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmgplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm', color='blue', horizontalalignment='left',
                             verticalalignment='top', transform=ax1.transAxes, fontsize=6)
                    ax1.text(1.0, 1.04,
                             'n = ' + str(len(lengthkmuplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmuplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm',
                             color='red', horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes,
                             fontsize=6)
                    ax1.tick_params(axis='y', labelcolor="black")

                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.show()
                    plt.savefig(
                        outfiledirectory + '\\' + 'Häufigkeit und Länge der Tiefenunterschreitungen entlang der tiefsten Rinne' + " TS" + str(
                            i) + '.' + plotformat)
                elif mindestwassertiefe_generell == 0 and mindestwassertiefe_untiefen > 0.0:
                    ax1.set_axisbelow(True)
                    ax1.grid(linestyle='-', linewidth=0.5)
                    ax1.scatter(meankmuplot, lengthkmuplot, s=4,
                                label="Mindestwassertiefe für Unftiefen <" + str(mindestwassertiefe_untiefen) + 'm',
                                color=color1, zorder=2)
                    ax1.axhline(y=maximallpasseinzelfall, color='r', linestyle='dashed',
                                label='maximale Länge passierbarer Einzelfälle',
                                linewidth=0.8, zorder=1)
                    ax1.legend(loc='best', frameon=False, fontsize=6)
                    # plt.figtext(0.80, 0.08, 'n = ' + str(len(lengthkmgplot)) + "/" + str(len(lengthkmuplot)), color='black', weight='roman',size='x-small')
                    ax1.text(1.0, 1.04,
                             'n = ' + str(len(lengthkmuplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmuplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm',
                             color='red', horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes,
                             fontsize=6)
                    ax1.tick_params(axis='y', labelcolor="black")

                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.show()
                    plt.savefig(
                        outfiledirectory + '\\' + 'Häufigkeit und Länge der Tiefenunterschreitungen entlang der tiefsten Rinne' + " TS" + str(
                            i) + '.' + plotformat)
                elif mindestwassertiefe_generell > 0.0 and mindestwassertiefe_untiefen == 0:
                    ax1.set_axisbelow(True)
                    ax1.grid(linestyle='-', linewidth=0.5)
                    ax1.scatter(meankmgplot, lengthkmgplot, s=4,
                                label="Generelle Mindestwassertiefe <" + str(mindestwassertiefe_generell) + 'm',
                                color=color, zorder=2)
                    ax1.axhline(y=maximallpasseinzelfall, color='r', linestyle='dashed',
                                label='maximale Länge passierbarer Einzelfälle',
                                linewidth=0.8, zorder=1)
                    ax1.legend(loc='best', frameon=False, fontsize=6)
                    # plt.figtext(0.80, 0.08, 'n = ' + str(len(lengthkmgplot)) + "/" + str(len(lengthkmuplot)), color='black', weight='roman',size='x-small')
                    ax1.text(0.0, 1.04, 'n = ' + str(len(lengthkmgplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmgplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm', color='blue', horizontalalignment='left',
                             verticalalignment='top', transform=ax1.transAxes, fontsize=6)
                    ax1.tick_params(axis='y', labelcolor="black")

                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.show()
                    plt.savefig(
                        outfiledirectory + '\\' + 'Häufigkeit und Länge der Tiefenunterschreitungen entlang der tiefsten Rinne' + " TS" + str(
                            i) + '.' + plotformat)
                else:
                    ax1.set_axisbelow(True)
                    ax1.grid(linestyle='-', linewidth=0.5)
                    ax1.scatter(meankmgplot, lengthkmgplot, s=4,
                                label="Generelle Mindestwassertiefe <" + str(mindestwassertiefe_generell) + 'm',
                                color=color, zorder=2)
                    ax1.scatter(meankmuplot, lengthkmuplot, s=4,
                                label="Mindestwassertiefe für Unftiefen <" + str(mindestwassertiefe_untiefen) + 'm',
                                color=color1, zorder=2)
                    ax1.axhline(y=maximallpasseinzelfall, color='r', linestyle='dashed',
                                label='maximale Länge passierbarer Einzelfälle', linewidth=0.8, zorder=1)
                    ax1.legend(loc='best', frameon=False, fontsize=6)
                    # plt.figtext(0.80, 0.08, 'n = ' + str(len(lengthkmgplot)) + "/" + str(len(lengthkmuplot)), color='black', weight='roman',size='x-small')
                    ax1.text(0.0, 1.04, 'n = ' + str(len(lengthkmgplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmgplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm', color='blue', horizontalalignment='left',
                             verticalalignment='top', transform=ax1.transAxes, fontsize=6)
                    ax1.text(1.0, 1.04,
                             'n = ' + str(len(lengthkmuplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmuplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm',
                             color='red', horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes,
                             fontsize=6)
                    ax1.tick_params(axis='y', labelcolor="black")

                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.show()
                    plt.savefig(
                        outfiledirectory + '\\' + 'Häufigkeit und Länge der Tiefenunterschreitungen entlang der tiefsten Rinne' + " TS" + str(
                            i) + '.' + plotformat)
            elif maximallpassuntiefen > 0 and maximallpasseinzelfall == 0:
                if mindestwassertiefe_generell == 0 and mindestwassertiefe_untiefen == 0:
                    ax1.set_axisbelow(True)
                    ax1.grid(linestyle='-', linewidth=0.5)
                    ax1.scatter(meankmuplot, lengthkmuplot, s=4, label="Points", color=color1, zorder=2)
                    ax1.axhline(y=maximallpassuntiefen, color='r', linestyle='-',
                                label='maximale Länge passierbarer Untiefen',
                                linewidth=0.8, zorder=1)
                    ax1.legend(loc='best', frameon=False, fontsize=6)
                    # plt.figtext(0.80, 0.08, 'n = ' + str(len(lengthkmgplot)) + "/" + str(len(lengthkmuplot)), color='black', weight='roman',size='x-small')
                    ax1.text(0.0, 1.04, 'n = ' + str(len(lengthkmgplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmgplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm', color='blue', horizontalalignment='left',
                             verticalalignment='top', transform=ax1.transAxes, fontsize=6)
                    ax1.text(1.0, 1.04,
                             'n = ' + str(len(lengthkmuplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmuplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm',
                             color='red', horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes,
                             fontsize=6)
                    ax1.tick_params(axis='y', labelcolor="black")

                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.show()
                    plt.savefig(
                        outfiledirectory + '\\' + 'Häufigkeit und Länge der Tiefenunterschreitungen entlang der tiefsten Rinne' + " TS" + str(
                            i) + '.' + plotformat)
                elif mindestwassertiefe_generell == 0 and mindestwassertiefe_untiefen > 0.0:
                    ax1.set_axisbelow(True)
                    ax1.grid(linestyle='-', linewidth=0.5)
                    ax1.scatter(meankmuplot, lengthkmuplot, s=4,
                                label="Mindestwassertiefe für Unftiefen <" + str(mindestwassertiefe_untiefen) + 'm',
                                color=color1, zorder=2)
                    ax1.axhline(y=maximallpassuntiefen, color='r', linestyle='-',
                                label='maximale Länge passierbarer Untiefen',
                                linewidth=0.8, zorder=1)
                    ax1.legend(loc='best', frameon=False, fontsize=6)
                    # plt.figtext(0.80, 0.08, 'n = ' + str(len(lengthkmgplot)) + "/" + str(len(lengthkmuplot)), color='black', weight='roman',size='x-small')
                    ax1.text(1.0, 1.04,
                             'n = ' + str(len(lengthkmuplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmuplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm',
                             color='red', horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes,
                             fontsize=6)
                    ax1.tick_params(axis='y', labelcolor="black")

                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.show()
                    plt.savefig(
                        outfiledirectory + '\\' + 'Häufigkeit und Länge der Tiefenunterschreitungen entlang der tiefsten Rinne' + " TS" + str(
                            i) + '.' + plotformat)
                elif mindestwassertiefe_generell > 0.0 and mindestwassertiefe_untiefen == 0:
                    ax1.set_axisbelow(True)
                    ax1.grid(linestyle='-', linewidth=0.5)
                    ax1.scatter(meankmgplot, lengthkmgplot, s=4,
                                label="Generelle Mindestwassertiefe <" + str(mindestwassertiefe_generell) + 'm',
                                color=color, zorder=2)
                    ax1.axhline(y=maximallpassuntiefen, color='r', linestyle='-',
                                label='maximale Länge passierbarer Untiefen',
                                linewidth=0.8, zorder=1)
                    ax1.legend(loc='best', frameon=False, fontsize=6)
                    # plt.figtext(0.80, 0.08, 'n = ' + str(len(lengthkmgplot)) + "/" + str(len(lengthkmuplot)), color='black', weight='roman',size='x-small')
                    ax1.text(0.0, 1.04, 'n = ' + str(len(lengthkmgplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmgplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm', color='blue', horizontalalignment='left',
                             verticalalignment='top', transform=ax1.transAxes, fontsize=6)
                    ax1.tick_params(axis='y', labelcolor="black")

                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.show()
                    plt.savefig(
                        outfiledirectory + '\\' + 'Häufigkeit und Länge der Tiefenunterschreitungen entlang der tiefsten Rinne' + " TS" + str(
                            i) + '.' + plotformat)
                else:
                    ax1.set_axisbelow(True)
                    ax1.grid(linestyle='-', linewidth=0.5)
                    ax1.scatter(meankmgplot, lengthkmgplot, s=4,
                                label="Generelle Mindestwassertiefe <" + str(mindestwassertiefe_generell) + 'm',
                                color=color, zorder=2)
                    ax1.scatter(meankmuplot, lengthkmuplot, s=4,
                                label="Mindestwassertiefe für Unftiefen <" + str(mindestwassertiefe_untiefen) + 'm',
                                color=color1, zorder=2)
                    ax1.axhline(y=maximallpassuntiefen, color='r', linestyle='-',
                                label='maximale Länge passierbarer Untiefen', linewidth=0.8, zorder=1)
                    ax1.legend(loc='best', frameon=False, fontsize=6)
                    # plt.figtext(0.80, 0.08, 'n = ' + str(len(lengthkmgplot)) + "/" + str(len(lengthkmuplot)), color='black', weight='roman',size='x-small')
                    ax1.text(0.0, 1.04, 'n = ' + str(len(lengthkmgplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmgplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm', color='blue', horizontalalignment='left',
                             verticalalignment='top', transform=ax1.transAxes, fontsize=6)
                    ax1.text(1.0, 1.04,
                             'n = ' + str(len(lengthkmuplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmuplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm',
                             color='red', horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes,
                             fontsize=6)
                    ax1.tick_params(axis='y', labelcolor="black")

                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.show()
                    plt.savefig(
                        outfiledirectory + '\\' + 'Häufigkeit und Länge der Tiefenunterschreitungen entlang der tiefsten Rinne' + " TS" + str(
                            i) + '.' + plotformat)
            elif maximallpassuntiefen > 0 and maximallpasseinzelfall > 0:
                if mindestwassertiefe_generell == 0 and mindestwassertiefe_untiefen == 0:
                    ax1.set_axisbelow(True)
                    ax1.grid(linestyle='-', linewidth=0.5)
                    ax1.scatter(meankmuplot, lengthkmuplot, s=4, label="Points", color=color1, zorder=2)
                    ax1.axhline(y=maximallpassuntiefen, color='r', linestyle='-',
                                label='maximale Länge passierbarer Untiefen',
                                linewidth=0.8, zorder=1)
                    ax1.axhline(y=maximallpasseinzelfall, color='r', linestyle='dashed',
                                label='maximale Länge passierbarer Einzelfälle',
                                linewidth=0.8, zorder=1)
                    ax1.legend(loc='best', frameon=False, fontsize=6)
                    # plt.figtext(0.80, 0.08, 'n = ' + str(len(lengthkmgplot)) + "/" + str(len(lengthkmuplot)), color='black', weight='roman',size='x-small')
                    ax1.text(0.0, 1.04, 'n = ' + str(len(lengthkmgplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmgplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm', color='blue', horizontalalignment='left',
                             verticalalignment='top', transform=ax1.transAxes, fontsize=6)
                    ax1.text(1.0, 1.04,
                             'n = ' + str(len(lengthkmuplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmuplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm',
                             color='red', horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes,
                             fontsize=6)
                    ax1.tick_params(axis='y', labelcolor="black")

                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.show()
                    plt.savefig(
                        outfiledirectory + '\\' + 'Häufigkeit und Länge der Tiefenunterschreitungen entlang der tiefsten Rinne' + " TS" + str(
                            i) + '.' + plotformat)
                elif mindestwassertiefe_generell == 0 and mindestwassertiefe_untiefen > 0.0:
                    ax1.set_axisbelow(True)
                    ax1.grid(linestyle='-', linewidth=0.5)
                    ax1.scatter(meankmuplot, lengthkmuplot, s=4,
                                label="Mindestwassertiefe für Unftiefen <" + str(mindestwassertiefe_untiefen) + 'm',
                                color=color1, zorder=2)
                    ax1.axhline(y=maximallpassuntiefen, color='r', linestyle='-',
                                label='maximale Länge passierbarer Untiefen',
                                linewidth=0.8, zorder=1)
                    ax1.axhline(y=maximallpasseinzelfall, color='r', linestyle='dashed',
                                label='maximale Länge passierbarer Einzelfälle',
                                linewidth=0.8, zorder=1)
                    ax1.legend(loc='best', frameon=False, fontsize=6)
                    # plt.figtext(0.80, 0.08, 'n = ' + str(len(lengthkmgplot)) + "/" + str(len(lengthkmuplot)), color='black', weight='roman',size='x-small')
                    ax1.text(1.0, 1.04,
                             'n = ' + str(len(lengthkmuplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmuplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm',
                             color='red', horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes,
                             fontsize=6)
                    ax1.tick_params(axis='y', labelcolor="black")

                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.show()
                    plt.savefig(
                        outfiledirectory + '\\' + 'Häufigkeit und Länge der Tiefenunterschreitungen entlang der tiefsten Rinne' + " TS" + str(
                            i) + '.' + plotformat)
                elif mindestwassertiefe_generell > 0.0 and mindestwassertiefe_untiefen == 0:
                    ax1.set_axisbelow(True)
                    ax1.grid(linestyle='-', linewidth=0.5)
                    ax1.scatter(meankmgplot, lengthkmgplot, s=4,
                                label="Generelle Mindestwassertiefe <" + str(mindestwassertiefe_generell) + 'm',
                                color=color, zorder=2)
                    ax1.axhline(y=maximallpassuntiefen, color='r', linestyle='-',
                                label='maximale Länge passierbarer Untiefen',
                                linewidth=0.8, zorder=1)
                    ax1.axhline(y=maximallpasseinzelfall, color='r', linestyle='dashed',
                                label='maximale Länge passierbarer Einzelfälle',
                                linewidth=0.8, zorder=1)
                    ax1.legend(loc='best', frameon=False, fontsize=6)
                    # plt.figtext(0.80, 0.08, 'n = ' + str(len(lengthkmgplot)) + "/" + str(len(lengthkmuplot)), color='black', weight='roman',size='x-small')
                    ax1.text(0.0, 1.04, 'n = ' + str(len(lengthkmgplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmgplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm', color='blue', horizontalalignment='left',
                             verticalalignment='top', transform=ax1.transAxes, fontsize=6)
                    ax1.tick_params(axis='y', labelcolor="black")

                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.show()
                    plt.savefig(
                        outfiledirectory + '\\' + 'Häufigkeit und Länge der Tiefenunterschreitungen entlang der tiefsten Rinne' + " TS" + str(
                            i) + '.' + plotformat)
                else:
                    ax1.set_axisbelow(True)
                    ax1.grid(linestyle='-', linewidth=0.5)
                    ax1.axhline(y=maximallpassuntiefen, color='r', linestyle='-',
                                label='maximale Länge passierbarer Untiefen', linewidth=0.8, zorder=1)
                    ax1.axhline(y=maximallpasseinzelfall, color='r', linestyle='dashed',
                                label='maximale Länge passierbarer Einzelfälle',
                                linewidth=0.8, zorder=1)
                    ax1.scatter(meankmgplot, lengthkmgplot, s=4,
                                label="Generelle Mindestwassertiefe <" + str(mindestwassertiefe_generell) + 'm',
                                color=color, zorder=2)
                    ax1.scatter(meankmuplot, lengthkmuplot, s=4,
                                label="Mindestwassertiefe für Unftiefen <" + str(mindestwassertiefe_untiefen) + 'm',
                                color=color1, zorder=2)
                    ax1.legend(loc='best', frameon=False, fontsize=6)
                    # plt.figtext(0.80, 0.08, 'n = ' + str(len(lengthkmgplot)) + "/" + str(len(lengthkmuplot)), color='black', weight='roman',size='x-small')
                    ax1.text(0.0, 1.04, 'n = ' + str(len(lengthkmgplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmgplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm', color='blue', horizontalalignment='left',
                             verticalalignment='top', transform=ax1.transAxes, fontsize=6)
                    ax1.text(1.0, 1.04,
                             'n = ' + str(len(lengthkmuplot)) + '/' + str(len(thalwegarraydeltaz[:,4])) + ' | ' + 'd = ' + str(round(sum(lengthkmuplot), 1))+ '/' + str(round(max(thalwegarraydeltaz[:,10]),1)) + 'm',
                             color='red', horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes,
                             fontsize=6)
                    ax1.tick_params(axis='y', labelcolor="black")

                    fig.tight_layout()  # otherwise the right y-label is slightly clipped
                    plt.show()
                    plt.savefig(
                        outfiledirectory + '\\' + 'Häufigkeit und Länge der Tiefenunterschreitungen entlang der tiefsten Rinne' + " TS" + str(
                            i) + '.' + plotformat)


            plt.close("all")

            # Make lists for each attribute in order to be able to display several values of the same attribute in the same graphic.
            waterdepthlist.append(thalwegarraydeltaz[:, 4])
            flowvelocitylist.append(thalwegarraydeltaz[:, 5])
            inflowdischargelist.append(thalwegarraydeltaz[:, 7])
            dischargelist.append(thalwegarraydeltaz[:, 8])
            laengelist.append(thalwegarraydeltaz[:, 10])

            waterdepthprozentdislist.append(thalwegarrayprozentdepth[:, 13])
            waterdepthcumlist.append(thalwegarrayprozentdepth[:, 4])
            flowvelocityprozentdislist.append(thalwegarrayprozentvelocity[:, 13])
            flowvelocitycumlist.append(thalwegarrayprozentvelocity[:, 5])
            deltazprozentdislist.append(thalwegarrayprozentdeltaz[:, 13])
            deltazcumlist.append(thalwegarrayprozentdeltaz[:, 12])

            lengthkmglist.append(lengthkmg)
            lengthkmulist.append(lengthkmu)
            dischargelengthkmglist.append(lengthkmgl)
            dischargelengthkmulist.append(lengthkmul)


        except nx.NetworkXNoPath:
            print('Kein Pfad für TS ' + str(i))
            continue
        except nx.NetworkXError:
            print('Networkx Fehler für TS ' + str(i))
            continue
        except nx.NodeNotFound:
            print('Source- or Target Node nicht gefunden für TS ' + str(i))
            continue
        except KeyError:
            print('TS ' + str(i) + ' ist nicht verfuegbar')
            continue
    # Make the attributes of the function available for other functions
    thalweg.dischargelist = dischargelist
    thalweg.waterdepthlist = waterdepthlist
    thalweg.flowvelocitylist = flowvelocitylist
    thalweg.inflowdischargelist = inflowdischargelist
    thalweg.laengelist = laengelist
    thalweg.waterdepthprozentdislist = waterdepthprozentdislist
    thalweg.waterdepthcumlist = waterdepthcumlist
    thalweg.flowvelocityprozentdislist = flowvelocityprozentdislist
    thalweg.flowvelocitycumlist = flowvelocitycumlist
    thalweg.deltazprozentdislist = deltazprozentdislist
    thalweg.deltazcumlist = deltazcumlist
    thalweg.lengthkmglist = lengthkmglist
    thalweg.lengthkmulist = lengthkmulist
    thalweg.dischargelengthkmglist = dischargelengthkmglist
    thalweg.dischargelengthkmulist = dischargelengthkmulist
    thalweg.allpaths = allpaths
    thalweg.thalwegtimestepslist = thalwegtimestepslistnew
# END - finding the thalweg for different timesteps
# ------------------------------------------------------------------------------
def thalwegdischarge(outfiledirectory, mindestwassertiefe_generell, mindestwassertiefe_untiefen, maximallpassuntiefen, maximallpasseinzelfall, maximalgeschwindigkeit, maximalabsturzhoehe, plotformat):
    """
    # Function definition:
        In this function, different thalwege are plotted for different outflows in one graphic per attribute.
    # Input (1) needed:
        filepath for the outputdirectory
    # Input (2) needed:
        general minimum water depth for fish
    # Input (3) needed:
        minimum water depth for shallows for fish
    # Input (4) needed:
        maximum velocity for fish
    # Input (5) needed:
        maximum fall height for fish, e.g. at threshold buoats
    # Input (6) needed:
        memory format for graphics e.g. pdf or png
    # Input (7) needed:
        timesteps
    # Output created:
        returns graphs with the cumulative distribution of water depths, flow velocities and fall heights along the thalweg.
        As well as the water depths and flow velocities in the longitudinal profile of this deepest channel.
        In addition, the frequency and length of the undershoots of the minimum water depth.
    """
    # Loads the lists for the plots from the Thalweg function.
    dischargelist = thalweg.dischargelist
    waterdepthlist = thalweg.waterdepthlist
    flowvelocitylist = thalweg.flowvelocitylist
    inflowdischargelist = thalweg.inflowdischargelist
    laengelist = thalweg.laengelist
    waterdepthprozentdislist = thalweg.waterdepthprozentdislist
    waterdepthcumlist = thalweg.waterdepthcumlist
    flowvelocityprozentdislist = thalweg.flowvelocityprozentdislist
    flowvelocitycumlist = thalweg.flowvelocitycumlist
    deltazprozentdislist = thalweg.deltazprozentdislist
    deltazcumlist = thalweg.deltazcumlist
    lengthkmglist = thalweg.lengthkmglist
    lengthkmulist = thalweg.lengthkmulist
    dischargelengthkmglist = thalweg.dischargelengthkmglist
    dischargelengthkmulist = thalweg.dischargelengthkmulist

    thalwegtimestepslist = thalweg.thalwegtimestepslist
    anzahltimesteps = len(thalwegtimestepslist)
    timesteps = np.arange(0, anzahltimesteps, 1)

    ###Plots###
    #Water Depth
    fig, ax1 = plt.subplots()
    color = iter(cm.rainbow(np.linspace(0, 0.8, anzahltimesteps)))

    for t in timesteps:
        dis = laengelist[t]
        depth = waterdepthlist[t]
        c = next(color)
        plt.title('Wassertiefen entlang der tiefsten Rinne', fontsize=10)
        ax1.set_xlabel('Distanz (m)', color="black", fontsize=6)
        ax1.set_ylabel('Wassertiefe (m)', color="black", fontsize=6)
        ax1.plot(dis, depth, label='Wassertiefe '+ str(round(dischargelist[t][0],1)) + r" $m^3/s$", color=c, linewidth=0.8)
        ax1.legend(loc='best', frameon=False, fontsize=6)
        ax1.tick_params(axis='y', labelcolor="black")
        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        ax1.set_axisbelow(True)
        ax1.grid(linestyle='-', linewidth=0.5)
        plt.show()
    if mindestwassertiefe_generell == 0 and mindestwassertiefe_untiefen == 0:
        ax1.legend(loc='best', frameon=False, fontsize=6)
        plt.savefig(outfiledirectory + '\\' + 'Verschiedene Abflüsse Wassertiefen entlang der tiefsten Rinne' + " TS" + str(thalwegtimestepslist) + '.' + plotformat)
        plt.close("all")
        print("Keine generelle Mindestwassertiefe und keine Mindestwassertiefe für Untiefen definiert")
        print("Plot Wassertiefen entlang der tiefsten Rinne vollendet")
    elif mindestwassertiefe_generell > 0.0 and mindestwassertiefe_untiefen == 0:
        ax1.axhline(y=mindestwassertiefe_generell, color='r', linestyle='-', label='Generelle Mindestwassertiefe', linewidth=0.8)
        ax1.legend(loc='best', frameon=False, fontsize=6)
        plt.savefig(outfiledirectory + '\\' + 'Verschiedene Abflüsse Wassertiefen entlang der tiefsten Rinne' + " TS" + str(thalwegtimestepslist) + '.' + plotformat)
        plt.close("all")
        print("Keine Mindestwassertiefe für Untiefen definiert")
        print("Plot Wassertiefen entlang der tiefsten Rinne vollendet")
    elif mindestwassertiefe_generell == 0 and mindestwassertiefe_untiefen > 0:
        ax1.axhline(y=mindestwassertiefe_untiefen, color='r', linestyle='dashed', label='Mindestwassertiefe für Untiefen',linewidth=0.8)
        ax1.legend(loc='best', frameon=False, fontsize=6)
        plt.savefig(outfiledirectory + '\\' + 'Verschiedene Abflüsse Wassertiefen entlang der tiefsten Rinne' + " TS" + str(thalwegtimestepslist) + '.' + plotformat)
        plt.close("all")
        print("Keine generelle Mindestwassertiefe definiert")
        print("Plot Wassertiefen entlang der tiefsten Rinne vollendet")
    else:
        ax1.axhline(y=mindestwassertiefe_generell, color='r', linestyle='-', label='Generelle Mindestwassertiefe', linewidth=0.8)
        ax1.axhline(y=mindestwassertiefe_untiefen, color='r', linestyle='dashed', label='Mindestwassertiefe für Untiefen',linewidth=0.8)
        ax1.legend(loc='best', frameon=False, fontsize=6)
        plt.savefig(outfiledirectory + '\\' + 'Verschiedene Abflüsse Wassertiefen entlang der tiefsten Rinne' + " TS" + str(thalwegtimestepslist) + '.' + plotformat)
        plt.close("all")
        print("Plot Wassertiefen entlang der tiefsten Rinne vollendet")

    #Flow Velocity
    fig, ax1 = plt.subplots()
    color = iter(cm.rainbow(np.linspace(0, 0.8, anzahltimesteps)))
    for t in timesteps:
        dis = laengelist[t]
        velocity = flowvelocitylist[t]
        c = next(color)
        plt.title('Fliessgeschwindigkeiten entlang der tiefsten Rinne', fontsize=10)
        ax1.set_xlabel('Distanz (m)', color="black", fontsize=6)
        ax1.set_ylabel('Geschwindigkeit (m/s)', color="black", fontsize=6)
        ax1.plot(dis, velocity, label='Geschwindigkeit '+ str(round(dischargelist[t][0],1)) + r" $m^3/s$", color=c, linewidth=0.8)
        ax1.legend(loc='best', frameon=False, fontsize=6)
        ax1.tick_params(axis='y', labelcolor="black")
        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        ax1.set_axisbelow(True)
        ax1.grid(linestyle='-', linewidth=0.5)
        plt.show()
    if maximalgeschwindigkeit == 0:
        ax1.legend(loc='best', frameon=False, fontsize=6)
        plt.savefig(outfiledirectory + '\\' + 'Verschiedene Abflüsse Fliessgeschwindigkeiten entlang der tiefsten Rinne' + " TS" + str(thalwegtimestepslist) + '.' + plotformat)
        plt.close("all")
        print("Keine Maximalgeschwindigkeit definiert")
        print("Plot Fliessgeschwindigkeiten entlang der tiefsten Rinne vollendet")
    else:
        ax1.axhline(y=maximalgeschwindigkeit, color='r', linestyle='-', label='Maximalgeschwindigkeit', linewidth=0.8)
        ax1.legend(loc='best', frameon=False, fontsize=6)
        plt.savefig(outfiledirectory + '\\' + 'Verschiedene Abflüsse Fliessgeschwindigkeiten entlang der tiefsten Rinne' + " TS" + str(thalwegtimestepslist) + '.' + plotformat)
        plt.close("all")
        print("Plot Fliessgeschwindigkeiten entlang der tiefsten Rinne vollendet")

    #Cummulative distribution of water depth
    fig, ax1 = plt.subplots()
    color = iter(cm.rainbow(np.linspace(0, 0.8, anzahltimesteps)))
    for t in timesteps:
        disw = waterdepthprozentdislist[t]
        depth = waterdepthcumlist[t]
        c = next(color)
        plt.title('Kumulative Verteilungskurve der Wassertiefen entlang des Talwegs', fontsize=10)
        ax1.set_xlabel('Wassertiefe (m)', color="black", fontsize=6)
        ax1.set_ylabel('Anteil (%)', color="black", fontsize=6)
        ax1.plot(depth, disw, label='Wassertiefe '+ str(round(dischargelist[t][0],1)) + r" $m^3/s$", color=c, linewidth=0.8)
        ax1.legend(loc='best', frameon=False, fontsize=6)
        ax1.tick_params(axis='y', labelcolor="black")
        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        ax1.set_axisbelow(True)
        ax1.grid(linestyle='-', linewidth=0.5)
        plt.show()
    if mindestwassertiefe_generell == 0 and mindestwassertiefe_untiefen == 0:
        ax1.legend(loc='best', frameon=False, fontsize=6)
        plt.savefig(outfiledirectory + '\\' + 'Verschiedene Abflüsse Kumulative Verteilungskurve der Wassertiefen entlang des Talwegs' + " TS" + str(thalwegtimestepslist) + '.' + plotformat)
        plt.close("all")
        print("Keine generelle Mindestwassertiefe und keine Mindestwassertiefe für Untiefen definiert")
        print("Plot Kumulative Verteilung der Wassertiefen vollendet")
    elif mindestwassertiefe_generell > 0.0 and mindestwassertiefe_untiefen == 0:
        ax1.axvline(x=mindestwassertiefe_generell, color='r', linestyle='-', label='Generelle Mindestwassertiefe', linewidth=0.8)
        ax1.legend(loc='best', frameon=False, fontsize=6)
        plt.savefig(outfiledirectory + '\\' + 'Verschiedene Abflüsse Kumulative Verteilungskurve der Wassertiefen entlang des Talwegs' + " TS" + str(thalwegtimestepslist) + '.' + plotformat)
        plt.close("all")
        print("Keine Mindestwassertiefe für Untiefen definiert")
        print("Plot Kumulative Verteilung der Wassertiefen vollendet")
    elif mindestwassertiefe_generell == 0 and mindestwassertiefe_untiefen > 0:
        ax1.axvline(x=mindestwassertiefe_untiefen, color='r', linestyle='dashed', label='Mindestwassertiefe für Untiefen', linewidth=0.8)
        ax1.legend(loc='best', frameon=False, fontsize=6)
        plt.savefig(outfiledirectory + '\\' + 'Verschiedene Abflüsse Kumulative Verteilungskurve der Wassertiefen entlang des Talwegs' + " TS" + str(thalwegtimestepslist) + '.' + plotformat)
        plt.close("all")
        print("Keine generelle Mindestwassertiefe definiert")
        print("Plot Kumulative Verteilung der Wassertiefen vollendet")
    else:
        ax1.axvline(x=mindestwassertiefe_generell, color='r', linestyle='-', label='Generelle Mindestwassertiefe', linewidth=0.8)
        ax1.axvline(x=mindestwassertiefe_untiefen, color='r', linestyle='dashed', label='Mindestwassertiefe für Untiefen', linewidth=0.8)
        ax1.legend(loc='best', frameon=False, fontsize=6)
        plt.savefig(outfiledirectory + '\\' + 'Verschiedene Abflüsse Kumulative Verteilungskurve der Wassertiefen entlang des Talwegs' + " TS" + str(thalwegtimestepslist) + '.' + plotformat)
        plt.close("all")
        print("Plot Kumulative Verteilung der Wassertiefen vollendet")


    #Cummulative distribution of flow velocity
    fig, ax1 = plt.subplots()
    color = iter(cm.rainbow(np.linspace(0, 0.8, anzahltimesteps)))
    for t in timesteps:
        disv = flowvelocityprozentdislist[t]
        velocity = flowvelocitycumlist[t]
        c = next(color)
        plt.title('Kumulative Verteilungskurve der Fliessgeschwindigkeiten entlang des Talwegs', fontsize=10)
        ax1.set_xlabel('Fliessgeschwindigkeit (m/s)', color="black", fontsize=6)
        ax1.set_ylabel('Anteil (%)', color="black", fontsize=6)
        ax1.plot(velocity, disv, label='Fliessgeschwindigkeit '+ str(round(dischargelist[t][0],1)) + r" $m^3/s$", color=c, linewidth=0.8)
        ax1.legend(loc='best', frameon=False, fontsize=6)
        ax1.tick_params(axis='y', labelcolor="black")
        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        ax1.set_axisbelow(True)
        ax1.grid(linestyle='-', linewidth=0.5)
        plt.show()
    if maximalgeschwindigkeit == 0:
        ax1.legend(loc='best', frameon=False, fontsize=6)
        plt.savefig(outfiledirectory + '\\' + 'Verschiedene Abflüsse Kumulative Verteilungskurve der Fliessgeschindigkeiten entlang des Talwegs' + " TS" + str(thalwegtimestepslist) + '.' + plotformat)
        plt.close("all")
        print("Keine Maximalgeschwindigkeit definiert")
        print("Plot Kumulative Verteilung der Fliessgeschwindigkeiten vollendet")
    else:
        ax1.axvline(x=maximalgeschwindigkeit, color='r', linestyle='-', label='Maximalgeschwindigkeit', linewidth=0.8)
        ax1.legend(loc='best', frameon=False, fontsize=6)
        plt.savefig(outfiledirectory + '\\' + 'Verschiedene Abflüsse Kumulative Verteilungskurve der Fliessgeschindigkeiten entlang des Talwegs' + " TS" + str(thalwegtimestepslist) + '.' + plotformat)
        plt.close("all")
        print("Plot Kumulative Verteilung der Fliessgeschwindigkeiten vollendet")

    #Cummulative distribution of fall height
    fig, ax1 = plt.subplots()
    color = iter(cm.rainbow(np.linspace(0, 0.8, anzahltimesteps)))
    for t in timesteps:
        disd = deltazprozentdislist[t]
        deltaz = deltazcumlist[t]
        c = next(color)
        plt.title('Kumulative Verteilungskurve der Hindernishöhen entlang des Talwegs', fontsize=10)
        ax1.set_xlabel('Hindernishöhe (m)', color="black", fontsize=6)
        ax1.set_ylabel('Anteil (%)', color="black", fontsize=6)
        ax1.plot(deltaz, disd, label='Hindernishöhe '+ str(round(dischargelist[t][0],1)) + r" $m^3/s$", color=c, linewidth=0.8)
        ax1.legend(loc='best', frameon=False, fontsize=6)
        ax1.tick_params(axis='y', labelcolor="black")
        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        ax1.set_axisbelow(True)
        ax1.grid(linestyle='-', linewidth=0.5)
        plt.show()
    if maximalabsturzhoehe == 0:
        ax1.legend(loc='best', frameon=False, fontsize=6)
        plt.savefig(outfiledirectory + '\\' + 'Verschiedene Abflüsse Kumulative Verteilungskurve der Hindernishöhen entlang des Talwegs' + " TS" + str(thalwegtimestepslist) + '.' + plotformat)
        plt.close("all")
        print("Keine maximale Sprunghöhe definiert")
        print("Plot Kumulative Verteilung der Hindernishöhen vollendet")
    else:
        ax1.axvline(x=maximalabsturzhoehe, color='r', linestyle='-', label='Maximale Sprunghöhe', linewidth=0.8)
        ax1.legend(loc='best', frameon=False, fontsize=6)
        plt.savefig(outfiledirectory + '\\' + 'Verschiedene Abflüsse Kumulative Verteilungskurve der Hindernishöhen entlang des Talwegs' + " TS" + str(thalwegtimestepslist) + '.' + plotformat)
        plt.close("all")
        print("Plot Kumulative Verteilung der Hindernishöhen vollendet")

    #Frequency and length of the underruns of the general and minimum water depths for shallows along the deepest channel
    #Percentage of the length of the minimum water depth not reached.
    fig, ax1 = plt.subplots()
    tableinfos = []
    for t in timesteps:
        lengthkmgarray = np.asarray(lengthkmglist[t])
        lengthkmuarray = np.asarray(lengthkmulist[t])
        if lengthkmgarray.size > 0 and lengthkmuarray.size == 0:
            meankmgplot = lengthkmglist[t]
            lengthkmgplot = [dischargelengthkmglist[t]]*len(meankmgplot)
            color = 'tab:blue'
            color1 = 'tab:red'
            plt.title('Häufigkeit und Länge der Unterschreitungen der generellen und \n der Mindestwassertiefe für Untiefen entlang der tiefsten Rinne', fontsize=10)
            ax1.set_xlabel('Abfluss ($m^3/s$)', color="black", fontsize=6)
            ax1.set_ylabel('Länge (m)', color="black", fontsize=6)
            ax1.set_axisbelow(True)
            ax1.grid(linestyle='-', linewidth=0.5)
            ax1.scatter(lengthkmgplot, meankmgplot, s=4, color=color, zorder=2)
            tableinfoslist = [round(dischargelist[t][0],1),(len(lengthkmgplot),len(waterdepthlist[t])),(round(sum(meankmgplot), 1),round(max(laengelist[t]),1)),0,0]
            tableinfos.append(tableinfoslist)
            ax1.tick_params(axis='y', labelcolor="black")
            fig.tight_layout()  # otherwise the right y-label is slightly clipped
            plt.show()
        elif lengthkmgarray.size == 0 and lengthkmuarray.size > 0:
            meankmuplot = lengthkmulist[t]
            lengthkmuplot = [dischargelengthkmulist[t]]*len(meankmuplot)
            color = 'tab:blue'
            color1 = 'tab:red'
            plt.title('Häufigkeit und Länge der Unterschreitungen der generellen und \n der Mindestwassertiefe für Untiefen entlang der tiefsten Rinne', fontsize=10)
            ax1.set_xlabel('Abfluss ($m^3/s$)', color="black", fontsize=6)
            ax1.set_ylabel('Länge (m)', color="black", fontsize=6)
            ax1.set_axisbelow(True)
            ax1.grid(linestyle='-', linewidth=0.5)
            ax1.scatter(lengthkmuplot, meankmuplot, s=4, color=color1, zorder=2)
            tableinfoslist = [round(dischargelist[t][0],1),0,0,(len(lengthkmuplot),len(waterdepthlist[t])),(round(sum(meankmuplot), 1),round(max(laengelist[t]),1))]
            tableinfos.append(tableinfoslist)
            ax1.tick_params(axis='y', labelcolor="black")
            fig.tight_layout()  # otherwise the right y-label is slightly clipped
            plt.show()
        elif lengthkmgarray.size == 0 and lengthkmuarray.size == 0:
            meankmuplot = lengthkmulist[t]
            lengthkmuplot = [dischargelengthkmulist[t]]*len(meankmuplot)
            plt.title('Häufigkeit und Länge der Unterschreitungen der generellen und \n der Mindestwassertiefe für Untiefen entlang der tiefsten Rinne', fontsize=10)
            ax1.set_axisbelow(True)
            ax1.grid(linestyle='-', linewidth=0.5)
            ax1.set_xlabel('Abfluss ($m^3/s$)', color="black", fontsize=6)
            ax1.set_ylabel('Länge (m)', color="black", fontsize=6)
            tableinfoslist = [round(dischargelist[t][0],1),0,0,0,0]
            tableinfos.append(tableinfoslist)
            ax1.tick_params(axis='y', labelcolor="black")
            fig.tight_layout()  # otherwise the right y-label is slightly clipped
            plt.show()
        else:
            meankmgplot = lengthkmglist[t]
            meankmuplot = lengthkmulist[t]
            lengthkmgplot = [dischargelengthkmglist[t]]*len(meankmgplot)
            lengthkmuplot = [dischargelengthkmulist[t]]*len(meankmuplot)
            color = 'tab:blue'
            color1 = 'tab:red'
            plt.title('Häufigkeit und Länge der Unterschreitungen der generellen und \n der Mindestwassertiefe für Untiefen entlang der tiefsten Rinne', fontsize=10)
            ax1.set_xlabel('Abfluss ($m^3/s$)', color="black", fontsize=6)
            ax1.set_ylabel('Länge (m)', color="black", fontsize=6)
            ax1.set_axisbelow(True)
            ax1.grid(linestyle='-', linewidth=0.5)
            ax1.scatter(lengthkmgplot, meankmgplot, s=4, color=color, zorder=2)
            ax1.scatter(lengthkmuplot, meankmuplot, s=4, color=color1, zorder=2)
            tableinfoslist = [round(dischargelist[t][0],1),(len(lengthkmgplot),len(waterdepthlist[t])),(round(sum(meankmgplot), 1),round(max(laengelist[t]),1)),(len(lengthkmuplot),len(waterdepthlist[t])),(round(sum(meankmuplot), 1),round(max(laengelist[t]),1))]
            tableinfos.append(tableinfoslist)
            ax1.tick_params(axis='y', labelcolor="black")
            fig.tight_layout()  # otherwise the right y-label is slightly clipped
            plt.show()
    if maximallpassuntiefen == 0 and maximallpasseinzelfall == 0:
        print("Keine maximale Länge für passierbare Untiefen und Einzelfälle definiert")
    elif maximallpassuntiefen > 0 and maximallpasseinzelfall == 0:
        ax1.axhline(y=maximallpassuntiefen, color='r', linestyle='-', label='maximale Länge passierbarer Untiefen',
                    linewidth=0.8, zorder=1)
    elif maximallpassuntiefen == 0 and maximallpasseinzelfall > 0:
        ax1.axhline(y=maximallpasseinzelfall, color='r', linestyle='dashed', label='maximale Länge passierbarer Einzelfälle',
                    linewidth=0.8, zorder=1)
    elif maximallpassuntiefen > 0 and maximallpasseinzelfall > 0:
        ax1.axhline(y=maximallpasseinzelfall, color='r', linestyle='dashed', label='maximale Länge passierbarer Einzelfälle',
                    linewidth=0.8, zorder=1)
        ax1.axhline(y=maximallpassuntiefen, color='r', linestyle='-', label='maximale Länge passierbarer Untiefen',
                    linewidth=0.8, zorder=1)
    ax1.legend(loc='best', frameon=False, fontsize=6)
    columns = ('Abfluss', 'Anzahl Längen\nzwischen '+str(mindestwassertiefe_generell) + 'm' + ' und ' + str(mindestwassertiefe_untiefen)+'m', 'Summierte Längen [m]\nzwischen '+str(mindestwassertiefe_generell) + 'm' + ' und ' + str(mindestwassertiefe_untiefen)+'m', 'Anzahl Längen\nkleiner '+str(mindestwassertiefe_untiefen) + 'm' , 'Summierte Längen [m]\nkleiner '+str(mindestwassertiefe_untiefen)+'m')
    rows = thalwegtimestepslist
    colors = ["w", 'tab:blue', 'tab:blue', 'tab:red', 'tab:red']
    cell_text = tableinfos
    plt.table(cellText=cell_text,rowLabels=rows,colLabels=columns,fontsize=6,loc='bottom',bbox=[0.0,-0.5,1,0.3],colColours=colors)
    #Adjust layout to make room for the table:
    plt.subplots_adjust(bottom=0.3)
    plt.savefig(outfiledirectory + '\\' + 'Verschiedene Abflüsse Häufigkeit und Länge der Tiefenunterschreitungen entlang der tiefsten Rinne' + " TS" + str(thalwegtimestepslist) + '.' + plotformat)
    plt.close("all")
    print("Plot Häufigkeit und Länge der Unterchreitung vollendet")
    #Creating a dataframe object from listoftuples
    tablehl = pd.DataFrame(tableinfos,columns = ['Abfluss', 'Häufigkeit<'+str(mindestwassertiefe_generell) + " (blau)", 'Länge<'+str(mindestwassertiefe_generell)+"(m)"+ " (blau)", 'Häufigkeit<'+str(mindestwassertiefe_untiefen)+ " (rot)", 'Länge<'+str(mindestwassertiefe_untiefen)+"(m)"+ " (rot)"])
    tablehl.to_csv(outfiledirectory + '/' + "Häufigkeit und Länge der Tiefenunterschreitungen bei allen untersuchten Timesteps und Abflüssen" + '.csv', index=True, header=True, sep=',')
# END - plotting the thalweg for different timesteps
# ------------------------------------------------------------------------------
def thalwegcentrality(inflowdat, outflowdat, maximalgeschwindigkeit, maximalabsturzhoehe, mindestwassertiefe_einzelfall1, sourcenode, targetnode):
    print("Talwegcentrality")
    thalwegtimestepslist = startcentrality.timestepsint1
    thalwegtimestepslistnew = np.array([], dtype=int)
    #bedingung = startcentrality.bedingung
    #talweg = startcentrality.talweg
    waterdepthlist = []
    allpaths = []

    # Starts the loop to output the thalweg as csv for each desired timestep
    for i in thalwegtimestepslist:
        try:
            #print("Berechnungen für den Timestep " + str(i))
            H = nx.Graph(createGraph.G)
            emd = createGraph.emd
            Attribute_nodes = createGraph.Attribute_nodes
            Depth_TS = createGraph.Depth_TS
            Velocity_TS = createGraph.Velocity_TS
            Wse_TS = createGraph.Wse_TS
            nodesarray = createGraph.nodesarray
            depth = dict(zip(Attribute_nodes, Depth_TS[("TS" + str(i))]))
            nx.set_node_attributes(H, depth, 'Water_Depth')

            # Change the weight of the edges
            nx.get_node_attributes(H, 'Water_Depth')
            H.edges(data=True)

            # add velocity to nodes
            velocity = dict(zip(Attribute_nodes, Velocity_TS[("TS" + str(i))]))
            nx.set_node_attributes(H, velocity, 'Flow_Velocity')

            # add wse to nodes
            wse = dict(zip(Attribute_nodes, Wse_TS[("TS" + str(i))]))
            nx.set_node_attributes(H, wse, 'WSE')

            # Delte edges with maximalgeschwindigkeit
            for u, v, f in H.edges(data=True):
                if maximalgeschwindigkeit != 0:
                    if abs(((H.nodes[u]['Flow_Velocity']) + (
                    H.nodes[v]['Flow_Velocity'])) / 2) < maximalgeschwindigkeit:
                        f['weightv'] = abs(((H.nodes[u]['Flow_Velocity']) + (H.nodes[v]['Flow_Velocity'])) / 2)
                    else:
                        f['weightv'] = 10000000000000.0
                elif maximalgeschwindigkeit == 0:
                    f['weightv'] = (((H.nodes[u]['Flow_Velocity']) + (H.nodes[v]['Flow_Velocity'])) / 2)

            # Delete edges with a weight of 1000000000000 otherwise the path will take the path with the fewest nodes at the edge because the distances are larger and it has fewer nodes.
            selected_edgesv = [(u, v) for u, v, f in H.edges(data=True) if f['weightv'] == 10000000000000.0]
            H.remove_edges_from(selected_edgesv)

            # Delte edges with maximalabsturzhoehe
            for u, v, g in H.edges(data=True):
                if maximalabsturzhoehe != 0:
                    if abs(((H.nodes[u]['WSE']) - (H.nodes[v]['WSE'])) / math.sqrt(((((H.nodes[u]['pos'])[0]) - ((H.nodes[v]['pos'])[0])) ** 2) + (((((H.nodes[u]['pos'])[1]) - ((H.nodes[v]['pos'])[1]))) ** 2))) < maximalabsturzhoehe:
                        g['weighta'] = abs(((H.nodes[u]['WSE']) - (H.nodes[v]['WSE'])) / math.sqrt(((((H.nodes[u]['pos'])[0]) - ((H.nodes[v]['pos'])[0])) ** 2) + (((((H.nodes[u]['pos'])[1]) - ((H.nodes[v]['pos'])[1]))) ** 2)))
                    else:
                        g['weighta'] = 10000000000000.0
                elif maximalabsturzhoehe == 0:
                    g['weighta'] = abs(((H.nodes[u]['WSE']) - (H.nodes[v]['WSE'])) / math.sqrt(((((H.nodes[u]['pos'])[0]) - ((H.nodes[v]['pos'])[0])) ** 2) + (((((H.nodes[u]['pos'])[1]) - ((H.nodes[v]['pos'])[1]))) ** 2)))
            # Delete edges with a weight of 1000000000000 otherwise the path will take the path with the fewest nodes at the edge because the distances are larger and it has fewer nodes.
            selected_edgesw = [(u, v) for u, v, g in H.edges(data=True) if g['weighta'] == 10000000000000.0]
            H.remove_edges_from(selected_edgesw)

            # 1/depth = x and then x/depth = y and again y/depth and so on to give the edges the right weighting.
            for u, v, d in H.edges(data=True):
                # d['weight'] = ((G.nodes[u]['Water_Depth'])+(G.nodes[v]['Water_Depth']))/2
                if mindestwassertiefe_einzelfall1 != 0:
                    if (H.nodes[u]['Water_Depth']) > 0.005 and (H.nodes[v]['Water_Depth']) > 0.005 and (
                            (H.nodes[u]['Water_Depth']) + (
                    H.nodes[v]['Water_Depth']) / 2) >= mindestwassertiefe_einzelfall1:
                        zwg1 = 1 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg2 = zwg1 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg3 = zwg2 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg4 = zwg3 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg5 = zwg4 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg6 = zwg5 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg7 = zwg6 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg8 = zwg7 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg9 = zwg8 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg10 = zwg9 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg11 = zwg10 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg12 = zwg11 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg13 = zwg12 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg14 = zwg13 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg15 = zwg14 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg16 = zwg15 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        d['weight'] = zwg16 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    else:
                        d['weight'] = 10000000000000.0
                elif mindestwassertiefe_einzelfall1 == 0:
                    if (H.nodes[u]['Water_Depth']) > 0.005 and (H.nodes[v]['Water_Depth']) > 0.005 and (
                            (H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth']) / 2) >= 0.005:
                        zwg1 = 1 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg2 = zwg1 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg3 = zwg2 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg4 = zwg3 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg5 = zwg4 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg6 = zwg5 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg7 = zwg6 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg8 = zwg7 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg9 = zwg8 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg10 = zwg9 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg11 = zwg10 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg12 = zwg11 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg13 = zwg12 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg14 = zwg13 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg15 = zwg14 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        zwg16 = zwg15 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                        d['weight'] = zwg16 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    else:
                        d['weight'] = 10000000000000.0

            # Delete edges with a weight of 1000000000000 otherwise the path will take the path with the fewest nodes at the edge because the distances are larger and it has fewer nodes.
            selected_edges = [(u, v) for u, v, d in H.edges(data=True) if d['weight'] == 10000000000000.0]
            H.remove_edges_from(selected_edges)
            #H.edges(data=True)


            # nx.info(H,pos=emd)
            #print("Folgende Informationen zum angepassten Netzwerk")
            #print(nx.info(H))
            global nodesremoveint
            nodesremoveint = startcentrality.nodesremoveint
            H.remove_nodes_from(nodesremoveint)
            # Find shortest path
            allpaths1 = list(islice(nx.shortest_simple_paths(H, sourcenode, targetnode, weight="weight"), 1))
            allpaths.append(allpaths1)
            allpaths2 = np.asarray(allpaths1)
            thalwegtimestepslistnew = np.append(thalwegtimestepslistnew, i)

            # Adding waterdepth to nodesarray
            depthttimestep = Depth_TS[("TS" + str(i))]
            depthttimesteptranspose = depthttimestep.transpose()
            depthtoarray = np.insert(nodesarray, 4, values=depthttimesteptranspose, axis=1)

            # Adding flow velocity to nodesarray
            velocitytimestep = Velocity_TS[("TS" + str(i))]
            velocitytimesteptranspose = velocitytimestep.transpose()
            velocitytoarray = np.insert(depthtoarray, 5, values=velocitytimesteptranspose, axis=1)

            # Adding wse to nodesarray
            wsetimestep = Wse_TS[("TS" + str(i))]
            wsetimesteptranspose = wsetimestep.transpose()
            wsetoarray = np.insert(velocitytoarray, 6, values=wsetimesteptranspose, axis=1)

            # Adding inflow to nodesarray
            inflow = inflowdat
            inflow1 = np.loadtxt(inflow)
            dischargein = inflow1[:, 1]
            dischargein1 = np.full((1, len(nodesarray)), dischargein[i])

            dischargeintoarray = np.insert(wsetoarray, 7, values=dischargein1, axis=1)

            # Adding outflow to nodesarray
            outflow = outflowdat
            outflow1 = np.loadtxt(outflow)
            discharge = outflow1[:, 1]
            discharge1 = np.full((1, len(nodesarray)), discharge[i])

            dischargetoarray = np.insert(dischargeintoarray, 8, values=discharge1, axis=1)

            # Leave only selected nodes in array test
            idx = allpaths1
            thalwegarray = dischargetoarray[np.in1d(dischargetoarray[:, 0], idx)]

            # Adding the drawing position
            test = thalwegarray.astype(int)
            test1 = test[:, 0]
            test2 = allpaths1[0]
            test3 = np.asarray(test2)

            position = np.array([])
            for element in test1:
                result = np.where(test3 == element)
                position = np.append(position, result)

            positiontranspose = position.transpose()

            thalwegarraypos = np.insert(thalwegarray, 9, values=positiontranspose, axis=1)

            # Sort the array with the thalweg nodes
            thalwegarraysort = np.sort(thalwegarraypos.view('f8,f8,f8,f8,f8,f8,f8,f8,f8,f8'), order=['f9'],
                                       axis=0).view(np.float)

            # Adding laenge and distance between nodes to nodesarray
            laenge = np.array([0.0])
            distancexy = np.array([0.0])

            index = len(thalwegarraysort) - 1
            e = 0
            while e < index:
                distance12 = math.sqrt(((thalwegarraysort[:, 1][e] - thalwegarraysort[:, 1][e + 1]) ** 2) + (
                            (thalwegarraysort[:, 2][e] - thalwegarraysort[:, 2][e + 1]) ** 2))
                distanceall = laenge[e] + distance12
                laenge = np.append(laenge, distanceall)
                distancexy = np.append(distancexy, distance12)

                e = e + 1

            thalwegarraylaenge = np.insert(thalwegarraysort, 10, values=laenge, axis=1)
            thalwegarraydis = np.insert(thalwegarraylaenge, 11, values=distancexy, axis=1)

            # Adding deltaz to nodesarray
            deltaz = np.array([0.0])
            index1 = len(thalwegarraysort) - 1
            o = 0
            while o < index1:
                distancedeltaz = math.sqrt(((thalwegarraysort[:, 1][o] - thalwegarraysort[:, 1][o + 1]) ** 2) + (
                            (thalwegarraysort[:, 2][o] - thalwegarraysort[:, 2][o + 1]) ** 2))
                deltazwse = abs(thalwegarraysort[:, 6][o] - thalwegarraysort[:, 6][o + 1]) / distancedeltaz
                deltaz = np.append(deltaz, deltazwse)
                o = o + 1

            thalwegarraydeltaz = np.insert(thalwegarraydis, 12, values=deltaz, axis=1)

            waterdepthlist.append(thalwegarraydeltaz[:, 4])

        except nx.NetworkXNoPath:
            print('Kein Pfad für TS ' + str(t))
            continue
        except nx.NetworkXError:
            print('Networkx Fehler für TS ' + str(t))
            continue
        except nx.NodeNotFound:
            print('Source- or Target Node nicht gefunden für TS ' + str(t))
            continue
        except KeyError:
            print('TS ' + str(t) + ' ist nicht verfuegbar')
            continue

    thalwegcentrality.waterdepthlist = waterdepthlist
    thalwegcentrality.allpaths = allpaths
    thalwegcentrality.thalwegtimestepslist = thalwegtimestepslistnew
# END - finding the thalweg for different timesteps
# ------------------------------------------------------------------------------
def centralityanalysis(outfiledirectory,maximalgeschwindigkeit, maximalabsturzhoehe,mindestwassertiefe_untiefen, mindestwassertiefe_einzelfall,):
    #Centrality Analysis
    #Betweenness-Zentralität 1/depth und Load Zentralität mit Gewichtung
    #Betweenness-Zentralität (betweenness centrality):
    #Ein Knoten hat einen hohen Betweenness-Wert, wenn dieser Knoten Bestandteil besonders vieler kürzester Wege ist und die jeweiligen Paare wenige andere kürzeste Wege haben, auf der der Knoten nicht enthalten ist.
    #Für jedes Paar von Knoten wird daher der Anteil an kürzesten Wegen zwischen ihnen berechnet, die v enthalten. Diese Anteile werden für alle Paare von Knoten aufsummiert um die Betweennesszentralität von v zu berechnen.
    waterdepthlist = thalwegcentrality.waterdepthlist
    centralitytimestepslist = thalwegcentrality.thalwegtimestepslist
    centralitytimestepslist1 = centralitytimestepslist.tolist()
    centralitytimestepslistnew = np.array([], dtype=int)
    anzahltimesteps = len(centralitytimestepslist)
    timesteps = np.arange(0, anzahltimesteps, 1)
    #bedingung = startcentrality.bedingung
    #talweg = startcentrality.talweg
    #global maxdepthnodes


    for t in centralitytimestepslist:
        try:
            H = nx.Graph(createGraph.G)
            emd = createGraph.emd
            Attribute_nodes = createGraph.Attribute_nodes
            Depth_TS = createGraph.Depth_TS
            Velocity_TS = createGraph.Velocity_TS
            Wse_TS = createGraph.Wse_TS
            nodesarray = createGraph.nodesarray
            allpaths = thalwegcentrality.allpaths
            depth = dict(zip(Attribute_nodes, Depth_TS[("TS" + str(t))]))
            nx.set_node_attributes(H, depth, 'Water_Depth')
            #Change the weight of the edges
            nx.get_node_attributes(H, 'Water_Depth')
            H.edges(data=True)

            # Add velocity to nodes
            velocity = dict(zip(Attribute_nodes, Velocity_TS[("TS" + str(t))]))
            nx.set_node_attributes(H, velocity, 'Flow_Velocity')

            # Add wse to nodes
            wse = dict(zip(Attribute_nodes, Wse_TS[("TS" + str(t))]))
            nx.set_node_attributes(H, wse, 'WSE')

            # Change the weight of the edges
            nx.get_node_attributes(H, 'Water_Depth')
            nx.get_node_attributes(H, 'Flow_Velocity')
            nx.get_node_attributes(H, 'WSE')
            nx.get_node_attributes(H, 'pos')

            H.edges(data=True)

            # Delete edges with maximalgeschwindigkeit
            for u, v, d in H.edges(data=True):
                if maximalgeschwindigkeit != 0:
                    if (((H.nodes[u]['Flow_Velocity']) + (H.nodes[v]['Flow_Velocity'])) / 2) < maximalgeschwindigkeit:
                        d['weightv'] = (((H.nodes[u]['Flow_Velocity']) + (H.nodes[v]['Flow_Velocity'])) / 2)
                    else:
                        d['weightv'] = 10000000000000.0
                elif maximalgeschwindigkeit == 0:
                    d['weightv'] = (((H.nodes[u]['Flow_Velocity']) + (H.nodes[v]['Flow_Velocity'])) / 2)

            # Delete edges with a weight of 1000000000000 otherwise the path will take the path with the fewest nodes at the edge because the distances are larger and it has fewer nodes.
            selected_edgesv = [(u, v) for u, v, e in H.edges(data=True) if e['weightv'] == 10000000000000.0]
            H.remove_edges_from(selected_edgesv)

            # Delete edges with maximalabsturzhoehe
            for u, v, d in H.edges(data=True):
                if maximalabsturzhoehe != 0:
                    if ((H.nodes[u]['WSE']) - (H.nodes[v]['WSE'])) / (
                            math.sqrt((((H.nodes[u]['pos'])[0]) - ((H.nodes[v]['pos'])[0])) ** 2) + (
                            (((H.nodes[u]['pos'])[1]) - ((H.nodes[v]['pos'])[1])) ** 2)) < maximalabsturzhoehe:
                        d['weightw'] = ((H.nodes[u]['WSE']) - (H.nodes[v]['WSE'])) / (
                                math.sqrt((((H.nodes[u]['pos'])[0]) - ((H.nodes[v]['pos'])[0])) ** 2) + (
                                (((H.nodes[u]['pos'])[1]) - ((H.nodes[v]['pos'])[1])) ** 2))
                    else:
                        d['weightw'] = 10000000000000.0
                elif maximalabsturzhoehe == 0:
                    d['weightw'] = ((H.nodes[u]['WSE']) - (H.nodes[v]['WSE'])) / (
                            math.sqrt((((H.nodes[u]['pos'])[0]) - ((H.nodes[v]['pos'])[0])) ** 2) + (
                            (((H.nodes[u]['pos'])[1]) - ((H.nodes[v]['pos'])[1])) ** 2))

            # Delete edges with a weight of 1000000000000 otherwise the path will take the path with the fewest nodes at the edge because the distances are larger and it has fewer nodes.
            selected_edgesw = [(u, v) for u, v, e in H.edges(data=True) if e['weightw'] == 10000000000000.0]
            H.remove_edges_from(selected_edgesw)
            tindex = centralitytimestepslist1.index(t)
            maxdepthnodes = min(waterdepthlist[tindex])

            # 1/depth = x and then x/depth = y and again y/depth and so on to give the edges the right weighting.
            for u, v, d in H.edges(data=True):
                # d['weight'] = ((G.nodes[u]['Water_Depth'])+(G.nodes[v]['Water_Depth']))/2
                if (H.nodes[u]['Water_Depth']) > 0.005 and (H.nodes[v]['Water_Depth']) > 0.005 and (
                        (H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth']) / 2) >= mindestwassertiefe_einzelfall:
                    """
                    zwg1 = 1 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    zwg2 = zwg1 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    zwg3 = zwg2 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    zwg4 = zwg3 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    zwg5 = zwg4 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    zwg6 = zwg5 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    zwg7 = zwg6 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    zwg8 = zwg7 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    zwg9 = zwg8 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    zwg10 = zwg9 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    zwg11 = zwg10 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    zwg12 = zwg11 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    zwg13 = zwg12 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    zwg14 = zwg13 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    zwg15 = zwg14 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    zwg16 = zwg15 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    d['weightd'] = zwg16 / (((H.nodes[u]['Water_Depth']) + (H.nodes[v]['Water_Depth'])) / 2)
                    """
                    d['weightd'] = ((1 / ((H.nodes[u]['Water_Depth']))) + ((1 / (H.nodes[v]['Water_Depth'])))) / 2
                else:
                    d['weightd'] = 10000000000000.0

            # Delete edges with a weight of 1000000000000 otherwise the path will take the path with the fewest nodes at the edge because the distances are larger and it has fewer nodes.
            selected_edgesd = [(u, v) for u, v, e in H.edges(data=True) if e['weightd'] == 10000000000000.0]
            H.remove_edges_from(selected_edgesd)
            H.edges(data=True)

            tindex = centralitytimestepslist1.index(t)
            maxdepthnodes = min(waterdepthlist[tindex])
            maxdepth = Depth_TS[
                ("TS" + str(t))]  # change the time step in the bracket so that each time step is executed
            ## timestep only with the dry nodes
            global water_pointsbenetzt
            global water_pointseinzelfall
            global water_pointsuntiefen
            maxdepthrraysmall = np.asarray(maxdepth)
            maxdepth_transposesmall = maxdepthrraysmall.transpose()
            water_pointssmall = np.where(
                maxdepth_transposesmall < mindestwassertiefe_einzelfall)  # nodes die unter der Tiefe des tiefsten Nodes des Talwegs liegen
            water_pointsbenetzt = np.where(maxdepth_transposesmall > 0.0)  # nodes die benetzt sind

            water_pointsuntiefen = np.where(
                maxdepth_transposesmall < mindestwassertiefe_untiefen)  # nodes die unter der Tiefe für Untiefen liegen
            water_pointseinzelfall = np.where(
                maxdepth_transposesmall < mindestwassertiefe_einzelfall)  # nodes die unter der Tiefe für Untiefen liegen

            #print("Prozent der benetzten Fläche liegen zwischen der Mindestwassertiefe für Untiefen " + str(mindestwassertiefe_untiefen) + " und der Mindestwassertiefe für Einzelfälle " + str(mindestwassertiefe_einzelfall))

            water_points_array_small = np.asarray(water_pointssmall)
            water_points_transpose_small = water_points_array_small.transpose()
            # Matrix
            smallvaluelist = water_points_array_small.tolist()
            # List
            smalldepthnodes = smallvaluelist[0]
            # Increase each number in the list by one because the nodes begin with 1 instead of 0
            smalldepthnodes = [int(n + 1) for n in smalldepthnodes]
            H.remove_nodes_from(smalldepthnodes)

            nodesremoveint = startcentrality.nodesremoveint
            H.remove_nodes_from(nodesremoveint)

            # Adding waterdepth to nodesarray
            depthttimestep = Depth_TS[("TS" + str(t))]
            depthttimesteptranspose = depthttimestep.transpose()
            depthtoarray = np.insert(nodesarray, 4, values=depthttimesteptranspose, axis=1)
            """
            #Only ausgewählte nodes
            idxc = H.nodes
            centralityarray = depthtoarray[np.in1d(depthtoarray[:, 0], idxc)]
            len(idxc)
            """

            # Adding flowvelocity to nodesarray
            velocitytimestep = Velocity_TS[("TS" + str(t))]
            velocitytimesteptranspose = velocitytimestep.transpose()
            velocitytoarray = np.insert(depthtoarray, 5, values=velocitytimesteptranspose, axis=1)
            # Only ausgewählte nodes
            idxc = H.nodes
            centralityarray = velocitytoarray[np.in1d(velocitytoarray[:, 0], idxc)]
            len(idxc)

            # nx.draw_networkx(H, emd)
            # plt.axis('equal')
            # plt.show()
            num = len(H.nodes)
            num1 = round((num / 4) * 2)
            if (num1 % 2) == 0:
                nk = num1 / 2
                nkint = int(nk)
            else:
                nk = (num1 + 1) / 2
                nkint = int(nk)
            print(
                "Aufgrund ihrer Angaben zur Mindestwassertiefe, maximalen Geschwindigkeit und maximalen Hindernishöhe wurden Nodes und Edges, die ihre Angaben nicht erfüllen im TS " + str(
                    t) + " gelöscht")
            print("Folgende Informationen zum neuen Netzwerk beim TS " + str(t))
            print(nx.info(H))
            # Betweeness centrality
            print("Betweeness Zentralität für TS" + str(t) + " wird berechnet")
            betweenness_centrality = nx.betweenness_centrality(H, weight="weightd", normalized=True,
                                                               k=400)  # k so setzen, dass es ein Bruchteil der Nodes ist oder alle nodes berücksichtigen
            sorted_betweenness_centrality = sorted(betweenness_centrality.items(), key=lambda kv: kv[0], reverse=False)
            sbc = 1
            sorted_betweenness_centrality1 = [x[sbc] for x in sorted_betweenness_centrality]
            # Adding betweenness_centrality to centralityarray
            betweenness_centralityarray = np.insert(centralityarray, 6, values=sorted_betweenness_centrality1, axis=1)
            print("Betweeness Zentralität Berechnung abgeschlossen")

            # Betweeness centrality subset
            # betweenness_centrality_subset = nx.betweenness_centrality_subset(H, sources=[10429], targets=[18895],normalized=False, weight="weightd")
            print("Betweeness Zentralität Subset für TS" + str(t) + " wird berechnet")
            betweenness_centrality_subset = nx.betweenness_centrality_subset(H, sources=startcentrality.sourcenodesbetint, targets=startcentrality.targetnodesbetint, normalized=False, weight="weightd")
            sorted_betweenness_centrality_subset = sorted(betweenness_centrality_subset.items(), key=lambda kv: kv[0],
                                                          reverse=False)
            sbcs = 1
            sorted_betweenness_centrality_subset1 = [x[sbcs] for x in sorted_betweenness_centrality_subset]
            # Adding betweenness_centrality to centralityarray
            betweenness_centrality_subsetarray = np.insert(betweenness_centralityarray, 7,
                                                           values=sorted_betweenness_centrality_subset1, axis=1)

            print("Betweeness Zentralität Subset abgeschlossen")

            # Load centrality
            """
            nc = len(allpaths[timesteps[tindex]][0])+1
            print("Load Zentralität für TS"  + str(t) +  " wird berechnet")
            load_centrality = nx.load_centrality(H, weight="weightd", cutoff=nc)
            sorted_load_centrality = sorted(load_centrality.items(), key=lambda kv: kv[0], reverse=False)
            lc = 1
            load_centrality1 = [x[lc] for x in sorted_load_centrality]
            #Adding betweenness_centrality to centralityarray
            load_centralityarray = np.insert(betweenness_centrality_subsetarray, 7, values=load_centrality1, axis=1)
            print("Load Zentralität Berechnung abgeschlossen")
            """

            deptharray = betweenness_centrality_subsetarray[:, 4]
            velocityarray = betweenness_centrality_subsetarray[:, 5]
            depthneigung = 9.81 * deptharray
            depthsqrt = np.sqrt(depthneigung)
            froude = velocityarray / depthsqrt
            froudetranspose = froude.transpose()
            froudearray = np.insert(betweenness_centrality_subsetarray, 8, values=froudetranspose, axis=1)

            velocitydepthratio = velocityarray / deptharray
            velocitydepthratiotranspose = velocitydepthratio.transpose()
            velocitydepthratioarray = np.insert(froudearray, 9, values=velocitydepthratiotranspose, axis=1)

            # Output to csv test
            # np.savetxt(outfiledirectory + "\\" + "betweenness_centrality" + str(t) + '.csv', X=load_centralityarray,fmt='%.8f', delimiter=',')
            float_formatter = lambda x: "%.8f" % x
            np.set_printoptions(formatter={'float_kind': float_formatter})
            matrix1 = pd.DataFrame(data=velocitydepthratioarray,
                                   columns=["Node", "X", "Y", "Z", "Wassertiefe", "Fliessgeschwindigkeit",
                                            "betweeness_centrality", "betweeness_centrality_subset", "Froude",
                                            "velocitydepthratio"])
            matrix1.to_csv(outfiledirectory + '\\' + 'Zentralitätsanalyse-TS' + str(t) + '.csv', index=True,
                           header=True, sep=',')
            print("Zentralitätsanalyse für TS " + str(t) + " abgeschlossen und Outputs gespeichert")
        except nx.NetworkXError:
            print('Networkx Fehler für TS ' + str(t))
            continue
        except nx.NodeNotFound:
            print("für TS " + str(
                t) + " wurden Nodes aus den angegebenen Subset-Nodes nicht gefunden, da sie eine zum geringe Wassertiefe aufweisen")
            continue
        except KeyError:
            print('TS ' + str(t) + ' ist nicht verfuegbar')
            continue
# END - finding the most important nodes in the network
# ------------------------------------------------------------------------------


# ---------------------------------____GUI____------------------------------
top = Tk()
top.title("Fischdurchgängigkeit Tool")
top.geometry("1200x640")
frame1 = Frame(top,height=720,width=600,bg = "LightSkyBlue1", borderwidth=2)
frame1.place(x=0, y=30)
frame2 = Frame(top,height=30,width=600,bg = "LightSkyBlue2", borderwidth=2)
frame2.place(x=0, y=0)
frame3 = Frame(top,height=720,width=600,bg = "PaleTurquoise1", borderwidth=2)
frame3.place(x=600, y=30)
frame2 = Frame(top,height=30,width=600,bg = "PaleTurquoise2", borderwidth=2)
frame2.place(x=600, y=0)
color1 = "LightSkyBlue1"
color2 = "PaleTurquoise1"

T1 = Label(top, text='Massenkonservierung und Suche der Source- und Targetnodes',fg="black", font=12,bg = "LightSkyBlue2")
T1.place(x=0, y=0)
#Frame(top,height=100,width=100,bg = "RED", borderwidth=2)

def clickedinflow():
    messagebox.showinfo('inflow Datei Pfad', 'Querprofil, an dem das Wasser ins Modell fliesst\nWichtig: Gleiche timesteps wie .sol files')
def inflowdat():
    inflowdat_filename = filedialog.askopenfilename(initialdir = os.getcwd(), title="Wähle inflow file", filetypes=(("inflow file", "*.dat"), ("all files", "*.*")))
    L101.delete(1, END)  # Remove current text in entry
    L101.insert(0, inflowdat_filename)  # Insert the 'path'
    inflowdat.path = inflowdat_filename
L100 = Label(top, text='inflow.dat File Pfad:',bg = color1)
L100.place(x=10, y=40)
L101 = Entry(top, text="", width=75)
L101.place(x=10, y=60)
B100 = Button(top, text='Browse', command=inflowdat)
B100.place(x=480, y=57)
Binflow = Button(top, text="Info",command=clickedinflow)
Binflow.place(x=550, y=57)

def clickedoutflow():
    messagebox.showinfo('outflow Datei Pfad', 'Querprofil, an dem das Wasser aus dem Modell fliesst\nWichtig: Gleiche timesteps wie .sol files')
def outflowdat():
    outflowdat_filename = filedialog.askopenfilename(initialdir = os.getcwd(), title="Wähle outflow file", filetypes=(("outflow file", "*.dat"), ("all files", "*.*")))
    L201.delete(1, END)  # Remove current text in entry
    L201.insert(0, outflowdat_filename)  # Insert the 'path'
    outflowdat.path = outflowdat_filename
L200 = Label(top, text='outflow.dat File Pfad:',bg = color1)
L200.place(x=10, y=80)
L201 = Entry(top, text="", width=75)
L201.place(x=10, y=100)
B200 = Button(top, text='Browse', command=outflowdat)
B200.place(x=480, y=97)
Boutflow = Button(top, text="Info",command=clickedoutflow)
Boutflow.place(x=550, y=97)

def clickedout():
    messagebox.showinfo('Output Ordner', 'Ordner, für Sicherung des Massenkonservierungsplots')
def outfiledirectory():
    outfiledirectory_filename = filedialog.askdirectory(initialdir = os.getcwd(), title="Wähle den Output Ordner")
    L301.delete(1, END)  # Remove current text in entry
    L301.insert(0, outfiledirectory_filename)  # Insert the 'path'
    outfiledirectory.path = outfiledirectory_filename
L300 = Label(top, text='Output Ordner:',bg = color1)
L300.place(x=10, y=120)
L301 = Entry(top, text="", width=75)
L301.place(x=10, y=140)
B300 = Button(top, text='Browse', command=outfiledirectory)
B300.place(x=480, y=137)
Boutdir = Button(top, text="Info",command=clickedout)
Boutdir.place(x=550, y=137)

def clickedprozent():
    messagebox.showinfo('Prozentuale Abweichung inflow und outflow', 'Dadurch werden die timesteps für die angegebene Abweichung zwischen inflow und outflow definiert')
L2500 = Label(top,  text="Prozentuale Abweichung inflow und outflow",bg = color1)
L2500.place(x=10, y=160)
var = IntVar()
var.set(1)
data=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
cb = Combobox(top, values=data)
cb.place(x=10, y=180)
B2500 = Button(top, text="Info",command=clickedprozent)
B2500.place(x=200, y=180)

def startmassconservation():
    print("Start Massenkonservierung")
    inflow = L101.get()
    outflow = L201.get()
    outfolder = L301.get()
    prozentabstand = int(cb.get())
    massconservation(inflow, outflow, prozentabstand,'pdf', outfolder)
    print("Massenkonservierung abgeschlossen")
S1 = Button(top, text='START MASSENKONSERVIERUNG',command=startmassconservation,background="SteelBlue1")
S1.place(x=200, y=220)


#-------------------------------------------------------------------------------------------------


def clickededg():
    messagebox.showinfo('Edge Datei Pfad', 'Edge File mit allen Edges')
def edgfile():
    edgfile_filename = filedialog.askopenfilename(initialdir = os.getcwd(), title="Wähle .edg file", filetypes=(("edg file", "*.edg"), ("all files", "*.*")))
    L401.delete(1, END)  # Remove current text in entry
    L401.insert(0, edgfile_filename)  # Insert the 'path'
    edgfile.path = edgfile_filename
L400 = Label(top, text='.edg File Pfad:',bg = color1)
L400.place(x=10, y=240)
L401 = Entry(top, text="", width=75)
L401.place(x=10, y=260)
B400 = Button(top, text='Browse', command=edgfile)
B400.place(x=480, y=257)
Bedg = Button(top, text="Info",command=clickededg)
Bedg.place(x=550, y=257)

def clickedtwodm():
    messagebox.showinfo('2dm Datei Pfad', 'File mit Informationen zur Gerinnetopographie')
def twodm():
    twodm_filename = filedialog.askopenfilename(initialdir = os.getcwd(), title="Wähle .2dm file", filetypes=(("2dm file", "*.2dm"), ("all files", "*.*")))
    L501.delete(1, END)  # Remove current text in entry
    L501.insert(0, twodm_filename)  # Insert the 'path'
    twodm.path = twodm_filename
L500 = Label(top, text='.2dm File Pfad:',bg = color1)
L500.place(x=10, y=280)
L501 = Entry(top, text="", width=75)
L501.place(x=10, y=300)
B500 = Button(top, text='Browse', command=twodm)
B500.place(x=480, y=297)
Btwodm = Button(top, text="Info",command=clickedtwodm)
Btwodm.place(x=550, y=297)

def clickeddepth():
    messagebox.showinfo('ndsdepthsol Datei Pfad', 'Simulierte Wassertiefen\nWichtig: Gleiche timesteps wie .dat files')
def ndsdepthsol():
    ndsdepthsol_filename = filedialog.askopenfilename(initialdir = os.getcwd(), title="Wähle nds_depth.sol file", filetypes=(("sol file", "*.sol"), ("all files", "*.*")))
    L601.delete(1, END)  # Remove current text in entry
    L601.insert(0, ndsdepthsol_filename)  # Insert the 'path'
    ndsdepthsol.path = ndsdepthsol_filename
L600 = Label(top, text='nds_depth.sol File Pfad:',bg = color1)
L600.place(x=10, y=320)
L601 = Entry(top, text="", width=75)
L601.place(x=10, y=340)
B600 = Button(top, text='Browse', command=ndsdepthsol)
B600.place(x=480, y=337)
Bdepth = Button(top, text="Info",command=clickeddepth)
Bdepth.place(x=550, y=337)

def clickedvelocity():
    messagebox.showinfo('nds_abs_velocity Datei Pfad', 'Simulierte Fliessgeschwindigkeiten\nWichtig: Gleiche timesteps wie .dat files')
def ndsabsvelocitysol():
    ndsabsvelocitysol_filename = filedialog.askopenfilename(initialdir = os.getcwd(), title="Wähle nds_abs_velocity.sol file", filetypes=(("sol file", "*.sol"), ("all files", "*.*")))
    L701.delete(1, END)  # Remove current text in entry
    L701.insert(0, ndsabsvelocitysol_filename)  # Insert the 'path'
    ndsabsvelocitysol.path = ndsabsvelocitysol_filename
L700 = Label(top, text='nds_abs_velocity.sol File Pfad:',bg = color1)
L700.place(x=10, y=360)
L701 = Entry(top, text="", width=75)
L701.place(x=10, y=380)
B700 = Button(top, text='Browse', command=ndsabsvelocitysol)
B700.place(x=480, y=377)
Bvelocity = Button(top, text="Info",command=clickedvelocity)
Bvelocity.place(x=550, y=377)

def clickedwse():
    messagebox.showinfo('nds_wse Datei Pfad', 'Simulierte Wasserspiegel\nWichtig: Gleiche timesteps wie .dat files')
def ndswsesol():
    ndswsesol_filename = filedialog.askopenfilename(initialdir = os.getcwd(), title="Wähle nds_wse.sol Datei", filetypes=(("sol file", "*.sol"), ("all files", "*.*")))
    L801.delete(1, END)  # Remove current text in entry
    L801.insert(0, ndswsesol_filename)  # Insert the 'path'
    ndswsesol.path = ndswsesol_filename
L800 = Label(top, text='nds_wse.sol File Pfad:',bg = color1)
L800.place(x=10, y=400)
L801 = Entry(top, text="", width=75)
L801.place(x=10, y=420)
B800 = Button(top, text='Browse', command=ndswsesol)
B800.place(x=480, y=417)
Bwse = Button(top, text="Info",command=clickedwse)
Bwse.place(x=550, y=417)


def clicked1():
    messagebox.showinfo('Timestep', 'Zeitschritt für den die Berechnungen durchgeführt werden')
L900 = Label(top,  text="Timestep",bg = color1)
L900.place(x=10, y=460)
E900 = Entry(top, text="", bd =5, width=35)
E900.place(x=110, y=460)
E900.insert(END, 0)
L901 = Label(top,  text="Example: 33",bg = color1)
L901.place(x=370, y=460)
B900 = Button(top, text="Info",command=clicked1)
B900.place(x=550, y=460)

def clicked2():
    messagebox.showinfo('Remove Nodes', 'Nodes, welche aus dem Netzwerk entfernt werden\nWichtig: Leerschlag nach Koma')
L10000 = Label(top,  text="Remove Nodes",bg = color1)
L10000.place(x=10, y=490)
E10000 = Entry(top, bd =5, text="", width=35)
E10000.place(x=110, y=490)
E10000.insert(END, 0)
L10001 = Label(top,  text="Example: 19632, 28, 6398, 1",bg = color1)
L10001.place(x=370, y=490)
B10000 = Button(top, text="Info",command=clicked2)
B10000.place(x=550, y=490)

def clicked301():
    messagebox.showinfo('Sourcenodes', 'Nodes, aus denen der beste Node als Startpunkt für den simulierten Talweg berechnet wird\nWichtig: Leerschlag nach Koma')
L1000 = Label(top,  text="Sourcenodes",bg = color1)
L1000.place(x=10, y=520)
E1000 = Entry(top, bd =5, text="", width=35)
E1000.place(x=110, y=520)
E1000.insert(END, 0)
L1001 = Label(top,  text="Example: 33465, 10, 5430, 5",bg = color1)
L1001.place(x=370, y=520)
B1000 = Button(top, text="Info",command=clicked301)
B1000.place(x=550, y=520)

def clicked302():
    messagebox.showinfo('Targetnodes', 'Nodes, aus denen der beste Node als Endpunkt für den simulierten Talweg berechnet wird\nWichtig: Leerschlag nach Koma')
L1100 = Label(top,  text="Targetnodes",bg = color1)
L1100.place(x=10, y=550)
E1100 = Entry(top, bd =5, text="", width=35)
E1100.place(x=110, y=550)
E1100.insert(END, 0)
L1101 = Label(top,  text="Example: 31257, 39, 6531, 4",bg = color1)
L1101.place(x=370, y=550)
B1100 = Button(top, text="Info",command=clicked302)
B1100.place(x=550, y=550)

def startsourcetargetnodes():
    print("Start Suche Sourcenodes und Targetnodes")
    edgfilname = L401.get()
    meshfilename = L501.get()
    depthfilename = L601.get()
    velocityfilname = L701.get()
    wsefilname = L801.get()
    createGraph(edgfilname, meshfilename, depthfilename, velocityfilname, wsefilname)
    timestep = int(E900.get())
    nodess = E1000.get().split(", ")
    nodessint = list(int(l) for l in nodess)
    startsourcetargetnodes.nodessint = nodessint
    nodest = E1100.get().split(", ")
    nodestint = list(int(l) for l in nodest)
    startsourcetargetnodes.nodestint = nodestint
    nodesremove = E10000.get().split(", ")
    nodesremoveint = list(int(l) for l in nodesremove)
    startsourcetargetnodes.nodesremoveint = nodesremoveint
    sourcenodes(timestep)
    targetnodes(timestep)
    print("Suche nach Suche Sourcenodes und Targetnodes abgeschlossen")
S2 = Button(top, text='START SUCHE SOURCENODES UND TARGETNODES',command=startsourcetargetnodes,background="SteelBlue1")
S2.place(x=150, y=590)


#-------------------------------------------------------------------------------------------------


T2 = Label(top, text='Hydraulische Durchgängigkeitsparameter entlang des Talwegs und Zentralitätsanalyse',fg="black", font=12,bg = "PaleTurquoise2")
T2.place(x=600, y=0)

def clicked4():
    messagebox.showinfo('Generelle Mindestwassertiefe', 'Mindestwassertiefe, die grundsätzlich vorhanden sein muss, sodass sich Fische ohne beachtlichen Energieaufwand bewegen können')
L1200 = Label(top,  text="Generelle Mindestwassertiefe [m]",bg = color2)
L1200.place(x=610, y=40)
E1200 = Entry(top, bd =5)
E1200.place(x=850, y=40)
E1200.insert(END, 0)
L1201 = Label(top,  text="Beispiel: 0.4",bg = color2)
L1201.place(x=1000, y=40)
B1200 = Button(top, text="Info",command=clicked4)
B1200.place(x=1150, y=37)

def clicked5():
    messagebox.showinfo('Mindestwassertiefe für Untiefen', 'Unterschreitungen der generellen Mindestwassertiefe sind in räumlich begrenzten Abschnitten mit natürlichen Untiefen möglich')
L1300 = Label(top,  text="Mindestwassertiefe für Untiefen [m]",bg = color2)
L1300.place(x=610, y=70)
E1300 = Entry(top, bd =5)
E1300.place(x=850, y=70)
E1300.insert(END, 0)
L1301 = Label(top,  text="Beispiel: 0.2",bg = color2)
L1301.place(x=1000, y=70)
B1300 = Button(top, text="Info",command=clicked5)
B1300.place(x=1150, y=67)

def clicked14():
    messagebox.showinfo('Mindestwassertiefe für Einzelfälle', 'Unterschreitungen der Mindestwassertiefe für Untiefen sind als Sonderfälle über kurze Distanzen möglich')
L2800 = Label(top,  text="Mindestwassertiefe für Einzelfälle [m]",bg = color2)
L2800.place(x=610, y=100)
E2800 = Entry(top, bd =5)
E2800.place(x=850, y=100)
E2800.insert(END, 0)
L2801 = Label(top,  text="Beispiel: 0.07",bg = color2)
L2801.place(x=1000, y=100)
B2800 = Button(top, text="Info",command=clicked14)
B2800.place(x=1150, y=97)

def clicked6():
    messagebox.showinfo('Maximale Länge passierbarer Untiefen', 'Streckenlänge, die die Länge der Unterschreitung der generellen Mindestwassertiefe definiert')
L1400 = Label(top,  text="Maximale Länge passierbarer Untiefen [m]",bg = color2)
L1400.place(x=610, y=130)
E1400 = Entry(top, bd =5)
E1400.place(x=850, y=130)
E1400.insert(END, 0)
L1401 = Label(top,  text="Beispiel: 30",bg = color2)
L1401.place(x=1000, y=130)
B1400 = Button(top, text="Info",command=clicked6)
B1400.place(x=1150, y=127)

def clicked17():
    messagebox.showinfo('Maximale Länge passierbarer Einzelfälle', 'Streckenlänge, die die Länge der Unterschreitung der Mindestwassertiefe für Untiefen definiert')
L2900 = Label(top,  text="Maximale Länge passierbarer Einzelfälle [m]",bg = color2)
L2900.place(x=610, y=160)
E2900 = Entry(top, bd =5)
E2900.place(x=850, y=160)
E2900.insert(END, 0)
L2901 = Label(top,  text="Beispiel: 5",bg = color2)
L2901.place(x=1000, y=160)
B2900 = Button(top, text="Info",command=clicked17)
B2900.place(x=1150, y=157)


def clicked7():
    messagebox.showinfo('Maximale Geschwindigkeit', 'Maximale Geschwindigkeiten, bei denen eine Wanderung der Fische noch möglich ist')
L1500 = Label(top,  text="Maximale Geschindigkeit [m/s]",bg = color2)
L1500.place(x=610, y=190)
E1500 = Entry(top, bd =5)
E1500.place(x=850, y=190)
E1500.insert(END, 0)
L1501 = Label(top,  text="Beispiel: 1.4",bg = color2)
L1501.place(x=1000, y=190)
B1500 = Button(top, text="Info",command=clicked7)
B1500.place(x=1150, y=187)

def clicked8():
    messagebox.showinfo('Maximale Sprunghöhe', 'Maximale Sprunghöhen, bei denen eine Wanderung der Fische noch möglich ist')
L1600 = Label(top,  text="Maximale Sprunghöhe [m]",bg = color2)
L1600.place(x=610, y=220)
E1600 = Entry(top, bd =5)
E1600.place(x=850, y=220)
E1600.insert(END, 0)
L1601 = Label(top,  text="Beispiel: 0.5",bg = color2)
L1601.place(x=1000, y=220)
B1600 = Button(top, text="Info",command=clicked8)
B1600.place(x=1150, y=217)

def clicked9():
    messagebox.showinfo('Timesteps', 'Berechnung erfolgt für die angegebenen Zeitschritte')
L1700 = Label(top,  text="Timesteps",bg = color2)
L1700.place(x=610, y=250)
E1700 = Entry(top, bd =5)
E1700.place(x=850, y=250)
E1700.insert(END, 0)
L1701 = Label(top,  text="Example: 33, 44, 52",bg = color2)
L1701.place(x=1000, y=250)
B1700 = Button(top, text="Info",command=clicked9)
B1700.place(x=1150, y=247)

def clicked10():
    messagebox.showinfo('Sourcenode', 'Startpunkt für die Simulation des Talweges')
L1800 = Label(top,  text="Sourcenode",bg = color2)
L1800.place(x=610, y=280)
E1800 = Entry(top, bd =5)
E1800.place(x=850, y=280)
E1800.insert(END, 0)
L1801 = Label(top,  text="Example: 35879",bg = color2)
L1801.place(x=1000, y=280)
B1800 = Button(top, text="Info",command=clicked10)
B1800.place(x=1150, y=277)

def clicked11():
    messagebox.showinfo('Targetnode', 'Endpunkt für die Simulation des Talweges')
L1900 = Label(top,  text="Targetnode",bg = color2)
L1900.place(x=610, y=310)
E1900 = Entry(top, bd =5)
E1900.place(x=850, y=310)
E1900.insert(END, 0)
L1901 = Label(top,  text="Example: 13267",bg = color2)
L1901.place(x=1000, y=310)
B1900 = Button(top, text="Info",command=clicked11)
B1900.place(x=1150, y=307)

def clicked12():
    messagebox.showinfo('Plot Outputformat', 'Format der Ausgabeplots, bspw. .pdf oder .png und weitere')
L2000 = Label(top,  text="Plot Outputformat",bg = color2)
L2000.place(x=610, y=340)
E2000 = Entry(top, bd =5)
E2000.place(x=850, y=340)
E2000.insert(END, "pdf")
L2001 = Label(top,  text='Example: pdf',bg = color2)
L2001.place(x=1000, y=340)
B2000 = Button(top, text="Info",command=clicked12)
B2000.place(x=1150, y=337)

def clickedout1():
    messagebox.showinfo('Output Ordner', 'Ordner, in dem die Plots und Textfiles der Talweganalyse gespeichert werden')
def outfiledirectory1():
    outfiledirectory_filename1 = filedialog.askdirectory(initialdir = os.getcwd(), title="Wähle Output Ordner")
    L2101.delete(1, END)  # Remove current text in entry
    L2101.insert(0, outfiledirectory_filename1)  # Insert the 'path'
    outfiledirectory1.path = outfiledirectory_filename1
L2100 = Label(top, text='Output File Pfad:',bg = color2)
L2100.place(x=610, y=360)
L2101 = Entry(top, text="", width=75)
L2101.place(x=610, y=380)
B2100 = Button(top, text='Browse', command=outfiledirectory1)
B2100.place(x=1080, y=377)
Boutdir1 = Button(top, text="Info",command=clickedout1)
Boutdir1.place(x=1150, y=377)

def clickedbedingung():
    messagebox.showinfo('Talweg mit oder ohne Bedingungen berechnen', 'Auswahl, ob die angegebenen Ansprüche als Restriktionen in die Berechnung einfliessen (mit) oder ob sie einfach in den Plots abgebildet werden (ohne)\nWichtig: Bei der Durchführung mit der Auswahl mitBedingungen ist es möglic, dass es keinen Talweg gibt, da die Bedingungen keinen Talweg zulassen')
L2600 = Label(top,  text="Talweg mit oder ohne Bedingungen berechnen",bg = color2)
L2600.place(x=610, y=410)
var = StringVar()
var.set("ohneBedingungen")
data=("ohneBedingungen", "mitBedingungen")
cb1 = Combobox(top, values=data)
cb1.place(x=610, y=430)
B2600 = Button(top, text="Info",command=clickedbedingung)
B2600.place(x=800, y=430)

def clickedweg():
    messagebox.showinfo('Tiefster oder kürzester Talweg berechnen', 'Auswahl, ob tiefster oder kürzester Talweg berechnet wird\nWichtig: Die Gewichtung ist entscheidend und wird im Benutzerhandbuch erklärt')
L2700 = Label(top,  text="Tiefster oder kürzester Talweg berechnen",bg = color2)
L2700.place(x=910, y=410)
var1 = StringVar()
var1.set("tiefster")
data=("tiefsterTalweg", "kuerzesterTalweg")
cb2 = Combobox(top, values=data)
cb2.place(x=910, y=430)
B2700 = Button(top, text="Info",command=clickedweg)
B2700.place(x=1100, y=430)

def startthalweg():
    print("Start Hydraulische Durchgängigkeitsanalyse")
    edgfilname = L401.get()
    meshfilename = L501.get()
    depthfilename = L601.get()
    velocityfilname = L701.get()
    wsefilname = L801.get()
    createGraph(edgfilname, meshfilename, depthfilename, velocityfilname, wsefilname)
    outfolder1 = L2101.get()
    inflow = L101.get()
    outflow = L201.get()
    mindestwassertiefe_generell = float(E1200.get())
    mindestwassertiefe_untiefen = float(E1300.get())
    mindestwassertiefe_einzelfall = float(E2800.get())
    maximallpassuntiefen = float(E1400.get())
    maximallpasseinzelfall = float(E2900.get())
    maximalgeschwindigkeit = float(E1500.get())
    maximalabsturzhoehe = float(E1600.get())
    sourcenode = int(E1800.get())
    targetnode = int(E1900.get())
    plotformat = E2000.get()
    timesteps = E1700.get().split(", ")
    timestepsint = list(int(l) for l in timesteps)
    startthalweg.timestepsint = timestepsint
    bedingung = cb1.get()
    startthalweg.bedingung = bedingung
    talweg = cb2.get()
    startthalweg.talweg = talweg
    nodesremove = E10000.get().split(", ")
    nodesremoveint = list(int(l) for l in nodesremove)
    startthalweg.nodesremoveint = nodesremoveint
    thalweg(outfolder1, inflow, outflow, mindestwassertiefe_generell, mindestwassertiefe_untiefen, mindestwassertiefe_einzelfall, maximallpassuntiefen, maximallpasseinzelfall, maximalgeschwindigkeit, maximalabsturzhoehe, sourcenode, targetnode, plotformat)
    print("Verschiedene Talwege werden in je einer Grafik pro Attribut geplottet")
    thalwegdischarge(outfolder1, mindestwassertiefe_generell, mindestwassertiefe_untiefen, maximallpassuntiefen, maximallpasseinzelfall, maximalgeschwindigkeit, maximalabsturzhoehe, plotformat)
    print("Hydraulische Durchgängigkeitsanalyse abgeschlossen und Plots gespeichert")
S1 = Button(top, text='START HYDRAULISCHE DURCHGÄNGIGKEITSANALYSE',command=startthalweg,background="SteelBlue1")
S1.place(x=750, y=460)


#-------------------------------------------------------------------------------------------------


def clicked13():
    messagebox.showinfo('betweenness_centrality_subset Sourcenodes', 'Startnodes für die Zentralitätsanalyse')
L2200 = Label(top,  text="betweenness_centrality_subset Sourcenodes",bg = "PaleTurquoise1")
L2200.place(x=610, y=490)
E2200 = Entry(top, bd =5)
E2200.place(x=850, y=490)
E2200.insert(END, 0)
L2201 = Label(top,  text='Example: 33465, 10, 5430',bg = color2)
L2201.place(x=1000, y=490)
B2200 = Button(top, text="Info",command=clicked13)
B2200.place(x=1150, y=490)


def clicked13():
    messagebox.showinfo('betweenness_centrality_subset Targetnodes', 'Endnodes für die Zentralitätsanalyse')
L2300 = Label(top,  text="betweenness_centrality_subset Targetnodes",bg = "PaleTurquoise1")
L2300.place(x=610, y=520)
E2300 = Entry(top, bd =5)
E2300.place(x=850, y=520)
E2300.insert(END, 0)
L2301 = Label(top,  text='Example: 31257, 38, 6531',bg = color2)
L2301.place(x=1000, y=520)
B2300 = Button(top, text="Info",command=clicked13)
B2300.place(x=1150, y=520)


def clickedout2():
    messagebox.showinfo('Output Ordner', 'Ordner in dem das Textfile der Zentralitätsanalyse gespeichert wird')
def outfiledirectory2():
    outfiledirectory_filename2 = filedialog.askdirectory(initialdir = os.getcwd(), title="Wähle den Output Ordner")
    L2401.delete(1, END)  # Remove current text in entry
    L2401.insert(0, outfiledirectory_filename2)  # Insert the 'path'
    outfiledirectory2.path = outfiledirectory_filename2
L2400 = Label(top, text='Output File Pfad:',bg = color2)
L2400.place(x=610, y=560)
L2401 = Entry(top, text="", width=75)
L2401.place(x=610, y=580)
B2400 = Button(top, text='Browse', command=outfiledirectory2)
B2400.place(x=1080, y=577)
Boutdir2 = Button(top, text="Info",command=clickedout2)
Boutdir2.place(x=1150, y=577)

def startcentrality():
    print("Start Zentralitätsanalyse")
    print("Dieser Vorgang kann mehrere Minuten dauern")
    edgfilname = L401.get()
    meshfilename = L501.get()
    depthfilename = L601.get()
    velocityfilname = L701.get()
    wsefilname = L801.get()
    createGraph(edgfilname, meshfilename, depthfilename, velocityfilname, wsefilname)
    inflow = L101.get()
    outflow = L201.get()
    maximalgeschwindigkeit = float(E1500.get())
    maximalabsturzhoehe = float(E1600.get())
    sourcenode = int(E1800.get())
    targetnode = int(E1900.get())
    timesteps1 = E1700.get().split(", ")
    timestepsint1 = list(int(k) for k in timesteps1)
    startcentrality.timestepsint1 = timestepsint1
    bedingung1 = cb1.get()
    startcentrality.bedingung = bedingung1
    talweg1 = cb2.get()
    startcentrality.talweg = talweg1
    mindestwassertiefe_einzelfall1 = float(E2800.get())
    nodesremove = E10000.get().split(", ")
    nodesremoveint = list(int(l) for l in nodesremove)
    startcentrality.nodesremoveint = nodesremoveint
    thalwegcentrality(inflow, outflow, maximalgeschwindigkeit, maximalabsturzhoehe, mindestwassertiefe_einzelfall1, sourcenode, targetnode)


    outfolder2 = L2401.get()
    maximalgeschwindigkeit1 = float(E1500.get())
    maximalabsturzhoehe1 = float(E1600.get())
    mindestwassertiefe_untiefen1 = float(E1300.get())
    mindestwassertiefe_einzelfall2 = float(E2800.get())
    sourcenodesbet = E2200.get().split(", ")
    sourcenodesbetint = list(int(k) for k in sourcenodesbet)
    startcentrality.sourcenodesbetint = sourcenodesbetint
    targetnodesbet = E2300.get().split(", ")
    targetnodesbetint = list(int(p) for p in targetnodesbet)
    startcentrality.targetnodesbetint = targetnodesbetint
    bedingung2 = cb1.get()
    startcentrality.bedingung = bedingung2
    talweg2 = cb2.get()
    startcentrality.talweg = talweg2
    centralityanalysis(outfolder2, maximalgeschwindigkeit1,maximalabsturzhoehe1,mindestwassertiefe_untiefen1,mindestwassertiefe_einzelfall2)
    print("Zentralitätsanalyse abgeschlossen")
S2 = Button(top, text='START ZENTRALITÄTSANALYSE',command=startcentrality,background="SteelBlue1")
S2.place(x=800, y=610)


def close():
    top.destroy()
C1 = Button(top, text='Close', command=close)
C1.place(x=1150, y=610)

top.mainloop()



