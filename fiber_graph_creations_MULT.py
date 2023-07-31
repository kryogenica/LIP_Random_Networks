
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:21:27 2023

@author: bryant_avila
"""
import repair_direct_ar_prohibit as rep
import random_color_balanced_graphs as rfibers
import numpy as np
import networkx as nx
from math import ceil
import shutil

number_of_nodes_lb=50
number_of_nodes_ub=100

number_of_colors_lb=5
number_of_colors_ub=10

equal_N_n_colors=10
norml_N_n_colors=10

path='/Users/bryant_avila/Projects/Network_repairs/random_fiber_symmetric_nets/'

edges_removed=[0.01,0.05,0.1,0.2,0.3,0.4,0.5]

weights = [[1,20],[1,10],[1,2],[1,1],[2,1],[10,1],[20,1]]

HardFlag=True; RMOnly=False; InDegOneFlag=True;

ppath = path+'G_prime.graph.txt'
ipath = path+'G_incre.graph.txt'
cpath = path+'G.colors.txt'
pbath = None
outfile = path

bad_solutions = 1
# %%create text file to keep record of all solutions
with open("random_fiber_symmetric_nets/log_IMP_test.txt", "a") as file:
    file.write('Number_of_nodes\tNumber_of_colors\tAdd_penalty\tRem_penalty\tB_selfloops\tDistribution\tColor_distribution\tId\tNum_edges_orig\tNum_loops\tPercent_edges_rem\tEdges_rem\tRandom_rem_org\tRandom_rem_rep\tEdges_rem_Rand\tRand_colors\tRand_time\tIncremental_rem_org\tIncremental_rem_rep\tEdges_rem_Inc\tIncr_colors\tIncre_time\n')

# %%
for N in range(number_of_nodes_lb,number_of_nodes_ub+1,50):
    for n_colors in range(number_of_colors_lb,number_of_colors_ub+1):
        P=1.2*(np.log(n_colors)/n_colors);
        for B_Sloops in [False,True]:
            for dis in ['equal','normal']:
                if dis=='equal':
                    flower = [N//n_colors]*(n_colors) 
                    for i in range(0,N%n_colors):
                        flower[i%n_colors] = flower[i%n_colors] + 1
                        
                    for net_num in range(0,equal_N_n_colors):
                        
                        # create Base graph
                        B_mat, B_colors = rfibers.random_base(n_colors,P,N,B_Sloops,True)
                        # create Total Space graph
                        G_mat, G_colors = rfibers.base_anthesis(B_mat,flower,N,n_colors,False)
                        
                        with open(cpath, 'w') as f:
                            for k, v in zip(range(0,N),G_colors):
                                f.write(str(k)+'\t'+str(v)+'\n')
                        
                        
                        Num_loops = B_mat.diagonal().sum()
                        # if n_colors!=rfibers.FindMP(nx.from_numpy_matrix(G_mat,create_using=nx.DiGraph))[0]:
                        #     impossibles_e = impossibles_e + 1
                            
                        #convert matrix into networkx and save as tsv file
                        G = nx.from_numpy_matrix(G_mat,create_using=nx.DiGraph)
                        
                        indices = np.where(G_mat == 1)
                        ones_count = len(indices[0])
                        G_mat_incremental = G_mat.copy()
                        holder = 0
                        for e_r in edges_removed:
                            
                            G_mat_prime = rfibers.remove_ones(G_mat.copy(),ceil(ones_count*e_r))
                            G_mat_incremental = rfibers.remove_ones(G_mat_incremental,ceil(ones_count*e_r)-holder)
                                                      
                            #convert matrix into networkx and save as tsv file
                            G_prime = nx.from_numpy_matrix(G_mat_prime,create_using=nx.DiGraph)
                            G_incre = nx.from_numpy_matrix(G_mat_incremental,create_using=nx.DiGraph)
                            edges_G_prime = set(G_prime.edges())
                            edges_G_incre = set(G_incre.edges())
                            nx.write_edgelist( G_prime, ppath, delimiter='\t', data=["weight"])
                            nx.write_edgelist( G_incre, ipath, delimiter='\t', data=["weight"])
                            
                            A,B,C,D,E,F,G,H,I = rep.set_rmip(ppath,cpath,HardFlag,[],[],InDegOneFlag,RMOnly,pbath)
                            R,S,T,U,V,W,X,Y,Z = rep.set_rmip(ipath,cpath,HardFlag,[],[],InDegOneFlag,RMOnly,pbath)
                            for option in range(0,7):
                                rm_weight = weights[option][0]
                                add_weight = weights[option][1]
                                _,_,_,_,_,_,_,_,_,_,G_prime_solu,executionTime_Prime = rep.solve_and_write(ppath,cpath,rm_weight,add_weight,outfile,A,B,C,D,E,F,G,H,I,HardFlag,[],[],InDegOneFlag,RMOnly,pbath,Save_info=False,NetX=True)
                                
                                _,_,_,_,_,_,_,_,_,_,G_incre_solu,executionTime_Incre = rep.solve_and_write(ipath,cpath,rm_weight,add_weight,outfile,R,S,T,U,V,W,X,Y,Z,HardFlag,[],[],InDegOneFlag,RMOnly,pbath,Save_info=False,NetX=True)
                                
                                
                                
                                edges_G_prime_solu = set(G_prime_solu.edges())
                                edges_G_incre_solu = set(G_incre_solu.edges())
                                
                                if len(edges_G_prime_solu)>0:
                                    n_colors_G_prime_solu,_,_ = rfibers.FindMP(G_prime_solu);
                                    if n_colors_G_prime_solu != n_colors:
                                        nx.write_gexf(G_prime_solu, path+"bad_solutions/"+str(n_colors)+"_example_"+str(bad_solutions)+"_colors_"+str(n_colors_G_prime_solu)+".gexf")
                                        shutil.copyfile(cpath, path+"bad_solutions/"+str(n_colors)+"_example_"+str(bad_solutions)+"_colors_"+str(n_colors_G_prime_solu)+".tsv")
                                        bad_solutions = bad_solutions +1
                                else:
                                    n_colors_G_prime_solu=0
                                    
                                if len(edges_G_incre_solu)>0:
                                    n_colors_G_incre_solu,_,_ = rfibers.FindMP(G_incre_solu);
                                    if n_colors_G_incre_solu != n_colors:
                                        nx.write_gexf(G_incre_solu, path+"bad_solutions/"+str(n_colors)+"_example_"+str(bad_solutions)+"_colors_"+str(n_colors_G_incre_solu)+".gexf")
                                        shutil.copyfile(cpath, path+"bad_solutions/"+str(n_colors)+"_example_"+str(bad_solutions)+"_colors_"+str(n_colors_G_incre_solu)+".tsv")
                                        bad_solutions = bad_solutions +1
                                else:
                                    n_colors_G_incre_solu=0
                                
                                
                                random_rem_org = len(edges_G_prime.difference(edges_G_prime_solu))
                                random_rem_rep = len(edges_G_prime.symmetric_difference(edges_G_prime_solu)) - random_rem_org
                                
                                incremental_rem_org = len(edges_G_incre.difference(edges_G_incre_solu))
                                incremental_rem_rep = len(edges_G_incre.symmetric_difference(edges_G_incre_solu)) - incremental_rem_org
                                
                                with open("random_fiber_symmetric_nets/log_IMP_test.txt", "a") as file:
                                    file.write(str(N)+'\t'+str(n_colors)+'\t'+str(add_weight)+'\t'+str(rm_weight)+'\t'+str(B_Sloops)+'\t'+dis+'\t'+str(flower)+'\t'+str(net_num)+'\t'+str(ones_count)+'\t'+str(Num_loops)+'\t'+str(e_r)+'\t'+str(ceil(ones_count*e_r))+'\t'+str(random_rem_org)+'\t'+str(random_rem_rep)+'\t'+str(len(edges_G_prime_solu))+'\t'+str(n_colors_G_prime_solu)+'\t'+str(executionTime_Prime)+'\t'+str(incremental_rem_org)+'\t'+str(incremental_rem_rep)+'\t'+str(len(edges_G_incre_solu))+'\t'+str(n_colors_G_incre_solu)+'\t'+str(executionTime_Incre)+'\n')
                            
                            holder = ceil(ones_count*e_r)
                        
                elif dis=='normal':
                    for fl in range(0,norml_N_n_colors):
                        
                        # create Base graph
                        B_mat, B_colors = rfibers.random_base(n_colors,P,N,B_Sloops,True)
                        #used to create the contraint random multiplicity of nodes in the general graph
                        flower = rfibers.constrained_sum_sample_pos(n_colors,N,list(np.diag(B_mat)))
                        # create Total Space graph
                        G_mat, G_colors = rfibers.base_anthesis(B_mat,flower,N,n_colors,False)

                        with open(cpath, 'w') as f:
                            for k, v in zip(range(0,N),G_colors):
                                f.write(str(k)+'\t'+str(v)+'\n')

                        Num_loops = B_mat.diagonal().sum()
                        
                        G = nx.from_numpy_matrix(G_mat,create_using=nx.DiGraph)
                        
                        indices = np.where(G_mat == 1)
                        ones_count = len(indices[0])
                        G_mat_incremental = G_mat.copy()
                        holder = 0
                        for e_r in edges_removed:
                            
                            G_mat_prime = rfibers.remove_ones(G_mat.copy(),ceil(ones_count*e_r))
                            G_mat_incremental = rfibers.remove_ones(G_mat_incremental,ceil(ones_count*e_r)-holder)
                                                      
                            #convert matrix into networkx and save as tsv file
                            G_prime = nx.from_numpy_matrix(G_mat_prime,create_using=nx.DiGraph)
                            G_incre = nx.from_numpy_matrix(G_mat_incremental,create_using=nx.DiGraph)
                            edges_G_prime = set(G_prime.edges())
                            edges_G_incre = set(G_incre.edges())
                            nx.write_edgelist( G_prime, ppath, delimiter='\t', data=["weight"])
                            nx.write_edgelist( G_incre, ipath, delimiter='\t', data=["weight"])
                            
                            A,B,C,D,E,F,G,H,I = rep.set_rmip(ppath,cpath,HardFlag,[],[],InDegOneFlag,RMOnly,pbath)
                            R,S,T,U,V,W,X,Y,Z = rep.set_rmip(ipath,cpath,HardFlag,[],[],InDegOneFlag,RMOnly,pbath)
                            for option in range(0,7):
                                rm_weight = weights[option][0]
                                add_weight = weights[option][1]
                            
                                _,_,_,_,_,_,_,_,_,_,G_prime_solu,executionTime_Prime = rep.solve_and_write(ppath,cpath,rm_weight,add_weight,outfile,A,B,C,D,E,F,G,H,I,HardFlag,[],[],InDegOneFlag,RMOnly,pbath,Save_info=False,NetX=True)
                                
                                _,_,_,_,_,_,_,_,_,_,G_incre_solu,executionTime_Incre = rep.solve_and_write(ipath,cpath,rm_weight,add_weight,outfile,R,S,T,U,V,W,X,Y,Z,HardFlag,[],[],InDegOneFlag,RMOnly,pbath,Save_info=False,NetX=True)
                            
                                edges_G_prime_solu = set(G_prime_solu.edges())
                                edges_G_incre_solu = set(G_incre_solu.edges())
                                
                                if len(edges_G_prime_solu)>0:
                                    n_colors_G_prime_solu,_,_ = rfibers.FindMP(G_prime_solu);
                                    if n_colors_G_prime_solu != n_colors:
                                        nx.write_gexf(G_prime_solu, path+"bad_solutions/"+str(n_colors)+"_example_"+str(bad_solutions)+"_colors_"+str(n_colors_G_prime_solu)+".gexf")
                                        shutil.copyfile(cpath, path+"bad_solutions/"+str(n_colors)+"_example_"+str(bad_solutions)+"_colors_"+str(n_colors_G_prime_solu)+".tsv")
                                        bad_solutions = bad_solutions +1
                                else:
                                    n_colors_G_prime_solu=0
                                    
                                if len(edges_G_incre_solu)>0:
                                    n_colors_G_incre_solu,_,_ = rfibers.FindMP(G_incre_solu);
                                    if n_colors_G_incre_solu != n_colors:
                                        nx.write_gexf(G_incre_solu, path+"bad_solutions/"+str(n_colors)+"_example_"+str(bad_solutions)+"_colors_"+str(n_colors_G_incre_solu)+".gexf")
                                        shutil.copyfile(cpath, path+"bad_solutions/"+str(n_colors)+"_example_"+str(bad_solutions)+"_colors_"+str(n_colors_G_incre_solu)+".tsv")
                                        bad_solutions = bad_solutions +1
                                else:
                                    n_colors_G_incre_solu=0
                                
                                random_rem_org = len(edges_G_prime.difference(edges_G_prime_solu))
                                random_rem_rep = len(edges_G_prime.symmetric_difference(edges_G_prime_solu)) - random_rem_org
                                
                                incremental_rem_org = len(edges_G_incre.difference(edges_G_incre_solu))
                                incremental_rem_rep = len(edges_G_incre.symmetric_difference(edges_G_incre_solu)) - incremental_rem_org
                                
                                with open("random_fiber_symmetric_nets/log_IMP_test.txt", "a") as file:
                                    file.write(str(N)+'\t'+str(n_colors)+'\t'+str(add_weight)+'\t'+str(rm_weight)+'\t'+str(B_Sloops)+'\t'+dis+'\t'+str(flower)+'\t'+str(fl)+'\t'+str(ones_count)+'\t'+str(Num_loops)+'\t'+str(e_r)+'\t'+str(ceil(ones_count*e_r))+'\t'+str(random_rem_org)+'\t'+str(random_rem_rep)+'\t'+str(len(edges_G_prime_solu))+'\t'+str(n_colors_G_prime_solu)+'\t'+str(executionTime_Prime)+'\t'+str(incremental_rem_org)+'\t'+str(incremental_rem_rep)+'\t'+str(len(edges_G_incre_solu))+'\t'+str(n_colors_G_incre_solu)+'\t'+str(executionTime_Incre)+'\n')
                            
                            holder = ceil(ones_count*e_r)