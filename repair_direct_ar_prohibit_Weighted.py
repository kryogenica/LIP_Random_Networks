

import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from gurobipy import abs_
from gurobipy import quicksum
import pandas as pd
import itertools as itools
from collections import defaultdict
import time

##precision parameter
epsilon = .001

#DEPRECATED:
#percentage of edge weights allowed for weight_mods to be 
#WM_perc = 1

##bound on how much a weight can be changed
WM_bound = 2

#this might need to be set in readdata. Right now a fixed constant
#bound on how much a new weight can be set to
NW_bound = 1

charsep='\t'

def get_key_from_value(dictionary, value):
    for key, val in dictionary.items():
        if value in val:
            return key

def readdata(fname,colorfile,xlinks=None):
    
    #load starting network
    GraphData = pd.read_csv(fname,sep=charsep,index_col=[0,1],header=None)
    #print(fname)
    #directed graph!    

    #remove selfloops
    GraphData=GraphData[GraphData.index.get_level_values(0)!=GraphData.index.get_level_values(1)]

    #data validation -- check for duplicate edges
    idx = GraphData.index
    if max(idx.duplicated()):
        print('Duplicate edges found. Returning garbage from readdata! \n')
        return []
    
    #directed graph!
    EdgeDict = GraphData.to_dict()[2]
    
    
    Edges,EdgeWeights = gp.multidict(EdgeDict)
    

    #make the node set
    Nodes = []
    for tup in Edges:
        if tup[0] not in Nodes:
            Nodes.append(tup[0])
        if tup[1] not in Nodes:
            Nodes.append(tup[1])
    
    ##read in color set
    ctable=pd.read_csv(colorfile,index_col=0,sep=charsep,header=None)        
    cdict = ctable.to_dict()[1]
    
    

    #set up a list of colorsets
    colorsets = []        
    colordict = {}
    for c in set(cdict.values()):
        C = [i for i in cdict.keys() if cdict[i]==c]
        colorsets.append(C)
        colordict[c] = C
    
    
     
    colorpairs =[]
    for C in colorsets:
        for p in C:
            
            #just in case the nodes are disconnected in a graph
            if p not in Nodes:
                Nodes.append(p)
                
        for p,q in itools.combinations(C,2):
            colorpairs.append((p,q))            
            
    nc_tuples = []
    outer_imbalance_dict = defaultdict(dict)
    inner_imbalance_dict = defaultdict(dict)
    support_num=0;
    for C,D in itools.combinations(colorsets,2):
        for p in C:
            p_color = get_key_from_value(colordict,p)
            for q in D:
                support_num = support_num + 1;
                q_color = get_key_from_value(colordict,q)
                inner_imbalance_dict[p][q]=[p_color,q_color]
                base_colors = list(colordict.keys())
                base_colors.remove(p_color)
                base_colors.remove(q_color)
                outer_imbalance_dict[p][q]=base_colors
                
                for c in colordict.keys():
                    nc_tuples.append((p,q,c))
                    nc_tuples.append((q,p,c))
                
                

    #NodePairs is all pairs of nodes in both directions as we're using 
    #directed graphs
    NodePairs = []
    for p in Nodes:
        for q in Nodes:
            if p != q:
                NodePairs.append((p,q))
    
    if xlinks!=None:
        prohibited = pd.read_csv(xlinks,sep=charsep,index_col=[0,1],header=None)
        no_access = pd.concat([GraphData,prohibited]);
        
        #directed graph!
        non_existing_EdgeDict = no_access.to_dict()[2]
        
        Edges_to_avoid = non_existing_EdgeDict.copy()
        Edges_to_avoid.update(EdgeDict)
        
        avoid_Edges,EdgeWeights = gp.multidict(Edges_to_avoid)
        
        NotE = {(p,q):1 for (p,q) in NodePairs if (p,q) not in avoid_Edges}
        NotEdges,NEWeights = gp.multidict(NotE)
        
    else:
        NotE = {(p,q):1 for (p,q) in NodePairs if (p,q) not in Edges}
        NotEdges,NEWeights = gp.multidict(NotE)
    
    return Nodes,Edges,EdgeWeights,colorpairs,colorsets,NotEdges,colordict,nc_tuples,outer_imbalance_dict,inner_imbalance_dict, support_num

    ##create a MIP
def CreateRMIP(Nodes,Edges,EdgeWeights,colorpairs,colorsets,outer_imbalance_dict,inner_imbalance_dict,support_num,env, \
               NotEdges,colordict,nc_tuples,HardFlag,\
               FixedEdges,FixedNonEdges,RMOnly,InDegOneFlag,StrongMinFlag,WeightFlag):
    
    
    rmip = gp.Model(name='RepairKnown-Directed',env=env)

    #initialize edge variables
    
    #these indicate whether an existing edge is removed or not
    remove_edge=rmip.addVars(Edges,vtype=GRB.BINARY,name='remove_edge')
    
    if WeightFlag:
        #these are modifications to the current weight, assumed integer
        weight_mods = rmip.addVars(Edges,lb=-GRB.INFINITY, vtype=GRB.INTEGER,name='weight_mods')
        weight_mods_abs = rmip.addVars(Edges,lb=0,ub=WM_bound, vtype=GRB.INTEGER,name='weight_mods_abs')
        max_weight = max(EdgeWeights.values())

    #these are variables representing potential positive or negative color balances
    #not used if HardFlag=True
    node_balance_pos = rmip.addVars(colorpairs,lb=0.0,vtype=GRB.CONTINUOUS,name='node_balance_pos')
    node_balance_neg = rmip.addVars(colorpairs,lb=0.0,vtype=GRB.CONTINUOUS,name='node_balance_neg')
    max_nodebalance = rmip.addVar(lb=0.0,vtype=GRB.CONTINUOUS,name='max_nodebalance')
    
    #binary indicaators if an edge is added or not
    add_edge=rmip.addVars(NotEdges,vtype=GRB.BINARY,name='add_edge')
    
    if WeightFlag:
        #if an edge is added, these variables represent the weight added
        #assuming this is integer, nonnegative
        new_weights=rmip.addVars(NotEdges,lb=0, ub=NW_bound, vtype=GRB.INTEGER,name="new_weights")
    
    
    #Bryant modified -- djp still needs to verify
    strict_balance = rmip.addVars(nc_tuples,vtype=GRB.BINARY,name='strict_balance')
    #auxiliary_var_1 = rmip.addVars(support_num,lb=-2,ub=2,vtype=GRB.SEMIINT,name='out_imbalance_one')
    #auxiliary_var_2 = rmip.addVars(support_num,lb=-2,ub=2,vtype=GRB.SEMIINT,name='out_imbalance_two')

    
    
    #dictionary to hold variables
    rvars = {'re':remove_edge,'nb_p':node_balance_pos,\
             'nb_n':node_balance_neg,'m_nb':max_nodebalance,\
                 'ae':add_edge,'sb':strict_balance}
    
    if WeightFlag:
        rvars['w_m']=weight_mods
        rvars['n_w']=new_weights
        rvars['w_m_a']=weight_mods_abs
        

    #constraint: colors in-balanced
    color_balance = []
    color_imbalance = []
    one_imbalance = []
    atleast_one = []
    indeg_one = []
    
    #bound constraints now controlled in variable declarations
    if WeightFlag:
        weight_abs_cons=[]
        weight_abs_cons.append(rmip.addConstrs((weight_mods_abs[i,j] == abs_(weight_mods[i,j]) \
                                               for [i,j] in Edges), name='weight_abs_bds'))

        
    #aux=[]
    n = len(Nodes)

    if InDegOneFlag and WeightFlag == False:
            indeg_one.append(rmip.addConstrs((sum((1-remove_edge[i,j]) for (i,j) in Edges if j == p) \
                                                + sum(add_edge[i,j] for (i,j) in NotEdges if j == p) >= 1 for p in Nodes), name='indeg_one'))

    if HardFlag:
        for D in colorsets:
            # color_balance.append(rmip.addConstrs(((sum((1-remove_edge[i,j]) for (i,j) in Edges \
            #                              if j == p and i in D) \
            #                         + sum(add_edge[i,j] for (i,j) in NotEdges \
            #                               if j == p and i in D) == \
            #                         sum((1-remove_edge[i,j]) for (i,j) in Edges \
            #                             if j == q and i in D) + \
            #                         sum(add_edge[i,j] for (i,j) in NotEdges \
            #                             if j == q and i in D)
            #                         ) for (p,q) in colorpairs), name='color_balance'))
            for (p,q) in colorpairs:                
                if WeightFlag:
                    #remove edge and add edge variables are not needed for the 
                    #weighted case the current way we do it
                    #relist_p = list(EdgeWeights[i,j]*(1-remove_edge[i,j]) for (i,j) in Edges if j == p and i in D)
                    relist_p = list(EdgeWeights[i,j] for (i,j) in Edges if j == p and i in D)
                    #Fcreating the list of weight_mods into p, D is the new weights into p
                    wmlist_p = list(weight_mods[i,j] for (i,j) in Edges if j == p)
                    nwlist_p = list(new_weights[i,j] for (i,j) in NotEdges if j == p)              
                    #relist_q = list(EdgeWeights[i,j]*(1-remove_edge[i,j]) for (i,j) in Edges if j == q and i in D)
                    relist_q = list(EdgeWeights[i,j] for (i,j) in Edges if j == q and i in D)
                    wmlist_q = list(weight_mods[i,j] for (i,j) in Edges if j == q)
                    nwlist_q = list(new_weights[i,j] for (i,j) in NotEdges if j == q)
                else:
                    #A and B are lists of removed and added edges into p
                    relist_p = list((1-remove_edge[i,j]) for (i,j) in Edges if j == p and i in D)
                    aelist_p = list(add_edge[i,j] for (i,j) in NotEdges if j == p and i in D)
                    #creating analogous lists for q
                    relist_q = list((1-remove_edge[i,j]) for (i,j) in Edges if j == q and i in D)
                    aelist_q = list(add_edge[i,j] for (i,j) in NotEdges if j == q and i in D)
                

                if WeightFlag: 
                    color_balance.append(rmip.addConstr(quicksum(relist_p) + quicksum(wmlist_p) + quicksum(nwlist_p)== \
                                          quicksum(relist_q) + quicksum(wmlist_q) + quicksum(nwlist_q), \
                                         name='color_balance'+str(p)+'_'+str(q)))
                        
                    
                else:
                    color_balance.append(rmip.addConstr((quicksum(relist_p) + quicksum(aelist_p) == quicksum(relist_q) + quicksum(aelist_q)), name='color_balance'+str(p)+'_'+str(q)))


                       
        #counter=0
        for C,D in itools.combinations(colorsets,2):
            for p in C:
                for q in D:
                    for c in colordict.keys():
                        #creating lists of removed and added edges into p
                        relist_p = list((1-remove_edge[i,j]) for (i,j) in Edges if j == p and i in colordict[c])
                        aelist_p = list(add_edge[i,j] for (i,j) in NotEdges if j == p and i in colordict[c])
                        if WeightFlag:
                        #creating the list of weight_mods into p, D is the new weights into p
                            relist_p = list(EdgeWeights[i,j]*(1-remove_edge[i,j]) for (i,j) in Edges if j == p and i in D)
                            wmlist_p = list(weight_mods[i,j] for (i,j) in Edges if j == p and i in colordict[c])
                            nwlist_p = list(new_weights[i,j] for (i,j) in NotEdges if j == p and i in colordict[c])

                        #a and b are lists of removed and addededges into q
                        #c is the list of weight mods into q, d is the new weights into q
                        relist_q = list((1-remove_edge[i,j]) for (i,j) in Edges if j == q and i in colordict[c])
                        aelist_q = list(add_edge[i,j] for (i,j) in NotEdges if j == q and i in colordict[c])
                        if WeightFlag:
                            relist_q = list(EdgeWeights[i,j]*(1-remove_edge[i,j]) for (i,j) in Edges if j == q and i in D)
                            wmlist_q = list(weight_mods[i,j] for (i,j) in Edges if j == q and i in colordict[c])
                            nwlist_q = list(new_weights[i,j] for (i,j) in NotEdges if j == q and i in colordict[c])
                        
                        
                        if WeightFlag:
                            color_imbalance.append(rmip.addConstr((quicksum(relist_p) + quicksum(wmlist_p) + quicksum(nwlist_p) >= \
                                                                   quicksum(relist_q) + quicksum(wmlist_q) + quicksum(nwlist_q) + \
                                                                   strict_balance[p,q,c] - n*max_weight*strict_balance[q,p,c]), name='imbalance_'+str(p)+'_'+str(q)+'_'+str(c)))
                            color_imbalance.append(rmip.addConstr((quicksum(relist_q) + quicksum(wmlist_q) + quicksum(nwlist_q) >= \
                                                                   quicksum(relist_p) + quicksum(wmlist_p) + quicksum(nwlist_p) + \
                                                                   strict_balance[q,p,c] - n*max_weight*strict_balance[p,q,c]), name='imbalance_'+str(q)+'_'+str(p)+'_'+str(c)))
                            
                        else:
                            color_imbalance.append(rmip.addConstr((quicksum(relist_p) + quicksum(aelist_p) >= \
                                                                   quicksum(relist_q) + quicksum(aelist_q) + \
                                                                   strict_balance[p,q,c] - n*strict_balance[q,p,c]), name='imbalance_'+str(p)+'_'+str(q)+'_'+str(c)))
                            color_imbalance.append(rmip.addConstr((quicksum(relist_q) + quicksum(aelist_q) >= \
                                                                   quicksum(relist_p) + quicksum(aelist_p) + \
                                                                   strict_balance[q,p,c] - n*strict_balance[p,q,c]), name='imbalance_'+str(q)+'_'+str(p)+'_'+str(c)))
                        
                        # color_imbalance.append(rmip.addConstrs(((quicksum((1-remove_edge[i,j]) for (i,j) in Edges \
                        #                                  if j == p and i in colordict[c]) \
                        #                             + quicksum(add_edge[i,j] for (i,j) in NotEdges \
                        #                                   if j == p and i in colordict[c]) >= \
                        #                             quicksum((1-remove_edge[i,j]) for (i,j) in Edges \
                        #                                 if j == q and i in colordict[c]) + \
                        #                             quicksum(add_edge[i,j] for (i,j) in NotEdges \
                        #                                 if j == q and i in colordict[c]) \
                        #                             + strict_balance[p,q,c] - \
                        #                             n*strict_balance[q,p,c]
                        #                             ) for c in colordict.keys()),name='imbalance_'+str(p)+'_'+str(q)))
                        
                            
                            
                        # color_imbalance.append(rmip.addConstrs(((quicksum((1-remove_edge[i,j]) for (i,j) in Edges \
                        #                                  if j == q and i in colordict[c]) \
                        #                             + quicksum(add_edge[i,j] for (i,j) in NotEdges \
                        #                                   if j == q and i in colordict[c]) >= \
                        #                             quicksum((1-remove_edge[i,j]) for (i,j) in Edges \
                        #                                 if j == p and i in colordict[c]) + \
                        #                             quicksum(add_edge[i,j] for (i,j) in NotEdges \
                        #                                 if j == p and i in colordict[c]) \
                        #                             + strict_balance[q,p,c] - \
                        #                             n*strict_balance[p,q,c]
                        #                             ) for c in colordict.keys()),name='imbalance_'+str(q)+'_'+str(p)))
                    
                        one_imbalance.append(rmip.addConstr((1 >= strict_balance[p,q,c] + \
                                                             strict_balance[q,p,c]) ,name='one_imbalance_'+str(p)+'_'+str(q)+'_'+str(c)))
                        
                    # one_imbalance.append(rmip.addConstrs(((1 >= strict_balance[p,q,i] + \
                    #                                      strict_balance[q,p,i]) for i in colordict.keys()),name='one_imbalance_'+str(p)+'_'+str(q)))
                    
                    if StrongMinFlag:
                        # this enforces a non-valid constraint but has stronger (although not complete) minimal properties
                        sblist = list(strict_balance[p,q,i]+strict_balance[q,p,i] for i in outer_imbalance_dict[p][q])
                    else:
                        # this enforces a valid constraint but still allows non-minimal                       
                        sblist = list(strict_balance[p,q,i] + strict_balance[q,p,i] for i in colordict.keys())
                    
                    if WeightFlag == False:
                        atleast_one.append(rmip.addConstr((quicksum(sblist) >= 1),name='atleast_one_'+str(p)+'_'+str(q)))
                    
                    
                    #this code will reject a correct solution 
                    #A = list(strict_balance[p,q,i] - strict_balance[q,p,i] for i in inner_imbalance_dict[p][q])
                    #rmip.addConstr(auxiliary_var_1[counter] == quicksum(A))
                    #rmip.addConstr(auxiliary_var_2[counter] == abs_(auxiliary_var_1[counter]))
                    
                    #B = list(strict_balance[p,q,i] + strict_balance[q,p,i] for i in outer_imbalance_dict[p][q])
                    #atleast_one.append(rmip.addConstr((quicksum(B) +\
                    #                   auxiliary_var_2[counter] >= 1),name='atleast_one_'+str(p)+'_'+str(q)))
                    #counter=counter+1
        
    else:
        for D in colorsets:
            color_balance.append(rmip.addConstrs((sum((1-remove_edge[i,j]) for (i,j) in Edges \
                                         if j == p and i in D) \
                                    + sum(add_edge[i,j] for (i,j) in NotEdges \
                                          if j == p and i in D) - \
                                    sum((1-remove_edge[i,j]) for (i,j) in Edges \
                                        if j == q and i in D) - \
                                    sum(add_edge[i,j] for (i,j) in NotEdges \
                                        if j == q and i in D)
                                         == \
                                    node_balance_pos[p,q] - node_balance_neg[p,q]\
                                    ) for (p,q) in colorpairs))

    FElist = []
    for (i,j) in FixedEdges:
        FElist.append(rmip.addConstr(remove_edge[i,j]==0))

    FNElist = []
    for (i,j) in FixedNonEdges:
        FNElist.append(rmip.addConstr(add_edge[i,j]==0))        

    if RMOnly:
        for (i,j) in NotEdges:
            FElist.append(rmip.addConstr(add_edge[i,j]==0))
    
    #keep track of edges/potential edges that are perturbed
    nodebalance_bounds_p = rmip.addConstrs((node_balance_pos[p,q] <= max_nodebalance \
                                          for (p,q) in colorpairs))
    nodebalance_bounds_n = rmip.addConstrs((node_balance_neg[p,q] <= max_nodebalance \
                                          for (p,q) in colorpairs))
    
      
          
    
    rcons={'cb':color_balance,'nb_b_p':nodebalance_bounds_p,\
           'nb_b_n':nodebalance_bounds_n,'FEl':FElist,'FNEl':FNElist,\
               'indeg_one':indeg_one}
        
    if WeightFlag:
        rcons['wa'] = weight_abs_cons

        
    return rmip,rcons,rvars,remove_edge,add_edge,node_balance_pos,node_balance_neg

def set_rmip(graphpath,colorpath,HardFlag,\
                 FixedEdges,FixedNonEdges,InDegOneFlag,RMOnly,prohibit,StrongMinFlag=False,WeightFlag=False):

    
    #create the inputs
    Nodes,Edges,EdgeWeights,ColorPairs,colorsets,NotEdges,colordict,nc_tuples,outer_imbalance_dict,inner_imbalance_dict,support_num = \
        readdata(graphpath,colorpath,prohibit)
    
    #set dictionary
    setdict = {'N':Nodes,'E':Edges,'CP':ColorPairs,'NE':NotEdges,'cd':colordict}
    
    if WeightFlag:
        setdict['EW'] = EdgeWeights
    
    #initialize an environment
    env = gp.Env()
    
    #create the model
    rmip,rcons,rvars,remove_edge,add_edge,node_balance_pos,node_balance_neg = CreateRMIP(Nodes,Edges,EdgeWeights,ColorPairs,colorsets,outer_imbalance_dict,inner_imbalance_dict,support_num,env,\
                   NotEdges,colordict,nc_tuples,HardFlag,FixedEdges,FixedNonEdges,RMOnly,InDegOneFlag,StrongMinFlag,WeightFlag)


    return rmip,rcons,rvars,setdict,colorsets,remove_edge,add_edge,node_balance_pos,node_balance_neg

def rmip_optimize(rmip,rcons,rvars,remove_edge,add_edge,node_balance_pos,\
                  node_balance_neg,rm_weight,add_weight,HardFlag,WeightFlag,\
                      wm_weight=1,nw_weight=1,bal_weight=1):
    
    #need objective
    if HardFlag:
        # #abs version
        # w = rmip.addVar(vtype=gp.GRB.INTEGER, name="w")  
        # rmip.addConstr(w >= rm_weight * gp.quicksum(remove_edge.select('*','*')) - add_weight * gp.quicksum(add_edge.select('*','*')))
        # rmip.addConstr(w >= - rm_weight * gp.quicksum(remove_edge.select('*','*')) + add_weight * gp.quicksum(add_edge.select('*','*')))
        
        if WeightFlag:
            weight_mods_abs  = rvars['w_m_a']
            new_weights = rvars['n_w']

#            obj = rm_weight*gp.quicksum(remove_edge.select('*','*')) + \
#                add_weight*gp.quicksum(add_edge.select('*','*')) + \
#                wm_weight*gp.quicksum(weight_mods_abs.select('*','*')) + \
#                nw_weight*gp.quicksum(new_weights.select('*','*'))
            
            obj = wm_weight*gp.quicksum(weight_mods_abs.select('*','*')) + \
                nw_weight*gp.quicksum(new_weights.select('*','*'))
        else:
            obj = rm_weight*gp.quicksum(remove_edge.select('*','*')) + \
                add_weight*gp.quicksum(add_edge.select('*','*'))
    
    else:
        obj = (epsilon + rm_weight)*(gp.quicksum(remove_edge.select('*','*'))) + \
            (epsilon + add_weight)*(gp.quicksum(add_edge.select('*','*'))) + \
            bal_weight*(gp.quicksum(node_balance_pos.select('*','*')) + \
               gp.quicksum(node_balance_neg.select('*','*'))) 
        
        
    # rmip.setObjective(w,GRB.MINIMIZE)
    rmip.setObjective(obj,GRB.MINIMIZE)

    #set the time limit -- not yet needed
    # rmip.setParam("TimeLimit", timelimit)

    #optimize
    startTime_Prime = time.time()
    rmip.optimize()
    executionTime = round(time.time() - startTime_Prime,5)

    return rmip,rcons,rvars,executionTime

#output file is fname
def solve_and_write(graphpath,colorpath,rm_weight,add_weight,fname,rmip,rcons,\
                       rvars,setdict,colorsets,remove_edge,add_edge,node_balance_pos,node_balance_neg,\
                       HardFlag=True,FixedEdges=[],FixedNonEdges=[],InDegOneFlag=True,\
                       RMOnly=False,prohibit=None,Save_info=True,NetX=False,WeightFlag=False,wm_weight=1,nw_weight=1):
    
    rmip,rcons,rvars,executionTime = rmip_optimize(rmip,rcons,rvars,remove_edge,\
                                                   add_edge,node_balance_pos,\
                                                   node_balance_neg,rm_weight,add_weight,\
                                                   HardFlag,WeightFlag,wm_weight,nw_weight,bal_weight=1)
    
    
    #find the edge removes
    E = setdict['E']
    NE = setdict['NE']

    cd = setdict['cd']

    re = rvars['re']
    ae = rvars['ae']
    sb = rvars['sb']
    
    
    if WeightFlag:
        EW = setdict['EW']
        w_m=rvars['w_m'] #weight modifications
        n_w=rvars['n_w'] #new weights
        w_m_a=rvars['w_m_a'] #absolute value of weight modificartions
    
    sumremovals = 0
    sumadds = 0
    idealnum=len(colorsets)
    feasible = (rmip.Status ==GRB.OPTIMAL)
    
    G_result = nx.DiGraph()

    if NetX==True:
        if feasible:
            for (i,j) in E:
                if WeightFlag:
                    if abs(EW[i,j] + w_m[i,j].x) > epsilon:
                        G_result.add_edge(i,j,weight=EW[i,j] +w_m[i,j].x)
                else:
                    if abs(re[i,j].x - 1) > epsilon:
                        G_result.add_edge(i, j)
    
            for (i,j) in NE:
                if WeightFlag:
                    if n_w[i,j].x > epsilon:
                        G_result.add_edge(i,j,weight=n_w[i,j].x)
                else:
                    if abs(ae[i,j].x - 1) < epsilon:
                        G_result.add_edge(i, j)
    
    if Save_info==True:
        outfname = fname+"directed.output.txt"
        f = open(outfname,"w")
        gname = fname+"directed.out.graph.txt"
        gf = open(gname,"w")
        
        if feasible:
            for (i,j) in E:
                if WeightFlag:
                    #sumremovals is the sum of absolute weight modification
                    sumremovals = sumremovals + w_m_a[i,j].x
                else:
                    if abs(re[i,j].x - 1) < epsilon:
                        sumremovals = sumremovals + 1

            for (i,j) in NE:
                if WeightFlag:
                    #sumadds is the sum of the new weights
                    sumadds = sumadds + n_w[i,j].x
                else:
                    if abs(ae[i,j].x - 1) < epsilon:
                        sumadds = sumadds + 1
        
    
                
        #print('Source Target Weight',file=gf)
        if WeightFlag:
            print(f'Total absolute edge modifications:\n{sumremovals}',file=f)
            print('Edges modified',file=f)            
        else:
            print(f'Total edges removed\n{sumremovals}',file=f)
            print('Edges removed',file=f)
        EdgesRemoved = []
        if feasible:
            for (i,j) in E:
                #print(f"weight mod:{i} {j} {w_m[i,j].x}")
                if WeightFlag:
                    if w_m_a[i,j].x > epsilon:
                        print(f'{i} {j} {w_m[i,j].x}',file=f)
                        EdgesRemoved.append((i,j))
                    #if abs(EW[i,j]+w_m[i,j].x) > epsilon:
                        #print(f'actual weight: {i} {j} {EW[i,j]+w_m[i,j].x}')
                else:
                    if abs(re[i,j].x - 1) < epsilon:
                        print(f'{i} {j}',file=f)
                        EdgesRemoved.append((i,j))
                    else:
                        print(f'{i} {j}',file=gf)
        if WeightFlag:
            print(f'Total weight of edges added\n{sumadds}',file=f)
            print('Edges added',file=f)            
        else:
            print(f'Total edges added\n{sumadds}',file=f)
            print('Edges added',file=f)
        EdgesAdded = []
        if feasible:
            for (i,j) in NE:
                if WeightFlag:
                    if n_w[i,j].x > epsilon:
                        print(f'{i} {j} {n_w[i,j].x}',file=f)
                        print(f'{i} {j} {n_w[i,j].x}',file=gf)   
                        EdgesAdded.append((i,j))
                else:
                    if abs(ae[i,j].x - 1) < epsilon:
                        print(f'{i} {j}',file=f)
                        print(f'{i} {j}',file=gf)   
                        EdgesAdded.append((i,j))
    
    
        CP = setdict['CP']
        m_nb = rvars['m_nb']
        nb_p = rvars['nb_p']
        nb_n = rvars['nb_n']
        if feasible:
            print(f'Maximum imbalance\n{m_nb.x}',file=f)            
        else:
            print("Maximum imbalance\n\n")
        print('Nonzero imbalances',file = f)    
        if feasible:
            for (i,j) in CP:
                imbalance = nb_p[i,j].x - nb_n[i,j].x
                if abs(imbalance) > epsilon:
                    print(f'{i} {j} {imbalance}',file=f)
                
        print('\nImbalances for each node and color',file=f)
        
        if feasible:        
            for C,D in itools.combinations(colorsets,2):
                for p in C:
                    for q in D:
                        print(f'Imbalances between {p} and {q}',file=f)
                        for i in cd:
                            if sb[p,q,i].x == 1 or sb[q,p,i].x == 1:
                                print(f'Color {i}',file=f)
    
        print("\n\n",end="",file=f)
        print("Input graph",file=f)
        GraphData = pd.read_csv(graphpath,sep=charsep,index_col=[0,1],header=None)
        GraphData.to_csv(f,sep=' ')
        
    
        print("Input colors",file = f)
        ctable=pd.read_csv(colorpath,index_col=0,sep=charsep,header=None)    
        ctable.to_csv(f,sep=' ')
        
        if prohibit!=None:
            print("Prohibited edges",file = f)
            prohibited = pd.read_csv(prohibit,sep=charsep,index_col=[0,1],header=None)
            prohibited.to_csv(f,sep=' ')
        
        f.close()
        gf.close()
    
        
    else:
        gname=[]; EdgesRemoved=[]; EdgesAdded=[]; outfname=[]; sumremovals=[]; sumadds=[];

    #return output file name and the number of partitions
    return gname,idealnum,EdgesRemoved,EdgesAdded,sumremovals,sumadds,outfname,rmip,rcons,rvars,G_result,executionTime

# ##main calls

# # root directory:
# workpath='/Users/bryant_avila/Projects/Network_repairs/'

# outpath=workpath+'outputs/celegans_LoS/Modularity_Maximization_LoS/'

# # list of directories to run files
# # dirlist = ['backward-chem','backward-gap','forward-chem','forward-gap']
# dirlist = ['celegans_LoS/Modularity_Maximization_LoS/']

# # alpha list
# rm_weights = [5]
# add_weights = [1]
# # for now gamma will be dependent on alpha

# # HardFlag = True if balancing must occur
# HardFlag = True

# for dp in dirlist:
#     dirpath = workpath + 'data/' + dp + '/'
#     graphfiles = [f for f in listdir(dirpath) if isfile(join(dirpath, f)) \
#                  and f.endswith('graph.txt')]

#     colorfiles = [f for f in listdir(dirpath) if isfile(join(dirpath, f)) \
#                  and f.endswith('colors.txt')]
        
#     for gf in graphfiles:
#         for cf in colorfiles:
#             for rm_weight in rm_weights:
#                 for add_weight in add_weights:
#                     gpath = dirpath+gf
#                     cpath = dirpath+cf
#                     if HardFlag:
#                         outfile = outpath+gf+'_'+cf+'.ar.balanced.'                    
#                     else:
#                         outfile = outpath+gf+'_'+cf+'a'+str(rm_weight)+'b'+str(add_weight) + '.ar.'
#                     write_one_solution(gpath,cpath,rm_weight,add_weight,outfile,HardFlag)
 

testpath = '/Users/phillips/Documents/test/DIRECTEDV1.graph.txt'
colorpath = '/Users/phillips/Documents/test/Names.colors.txt'
outpath = '/Users/phillips/Documents/test/out.txt'
HardFlag = True
InDegOneFlag=True
RMOnly = True
prohibit=None
#A,B,C,D,E,F,G,H,I = set_rmip(testpath,colorpath,HardFlag,[],[],InDegOneFlag,False,prohibit,True,False)
#solve_and_write(testpath,colorpath,1,1,outpath,A,B,C,D,E,F,G,H,I,HardFlag,[],[],InDegOneFlag,RMOnly,prohibit,Save_info=False,NetX=True)

#testwpath = '/Users/phillips/Documents/test/Cook_chem_back.graph.txt'
#wcolorpath = '/Users/phillips/Documents/test/Modularity_Maximization-LoS-Bryant_var_9.colors.txt'
#woutpath = '/Users/phillips/Documents/test/outw.txt'
testwpath = '/Users/phillips/Documents/test/small_test.txt'
wcolorpath = '/Users/phillips/Documents/test/small.colors.txt'
woutpath = '/Users/phillips/Documents/test/outw.txt'
A,B,C,D,E,F,G,H,I = set_rmip(testwpath,wcolorpath,HardFlag,[],[],InDegOneFlag,False,prohibit,WeightFlag=True,StrongMinFlag=False)
solve_and_write(testwpath,wcolorpath,1,1,woutpath,A,B,C,D,E,F,G,H,I,HardFlag,[],[],InDegOneFlag,RMOnly,prohibit,WeightFlag=True,Save_info=True,NetX=True)
