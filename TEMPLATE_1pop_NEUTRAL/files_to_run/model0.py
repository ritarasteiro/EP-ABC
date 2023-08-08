#TODO
from __future__ import print_function
# this line for I using the python 2.7 instead of 3,
# so some print function eg sep="\t" will not work, but which is very popular in python 3
import msprime
#import numpy as np
#import math
def model0(Length, theta, rou,n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,nanc,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,s1):
 population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=s1, initial_size=n0)
           ]
 demographic_events = [
        # Size changes to Ni at Ti
    msprime.PopulationParametersChange(time=t1, initial_size=n1,population_id=0),
    msprime.PopulationParametersChange(time=t2, initial_size=n2,population_id=0),
    msprime.PopulationParametersChange(time=t3, initial_size=n3,population_id=0),
    msprime.PopulationParametersChange(time=t4, initial_size=n4,population_id=0),
    msprime.PopulationParametersChange(time=t5, initial_size=n5,population_id=0),
    msprime.PopulationParametersChange(time=t6, initial_size=n6,population_id=0),
    msprime.PopulationParametersChange(time=t7, initial_size=n7,population_id=0),
    msprime.PopulationParametersChange(time=t8, initial_size=n8,population_id=0),
    msprime.PopulationParametersChange(time=t9, initial_size=n9,population_id=0),
    msprime.PopulationParametersChange(time=t10, initial_size=n10,population_id=0),
    msprime.PopulationParametersChange(time=t11, initial_size=nanc,population_id=0)
    ]

# Use the demography debugger to print out the demographic history
#    # that we have just described.
# dd = msprime.DemographyDebugger(
#    population_configurations=population_configurations,
#    demographic_events=demographic_events)
    # now try to use simulation
# dd.print_history()

 tree_seq = msprime.simulate( Ne= nanc,length = Length, recombination_rate = rou,
       mutation_rate = theta, population_configurations=population_configurations,demographic_events=demographic_events)
 outF1 = open("loci", "w")
 outF2 = open("position", "w")
 for variant in tree_seq.variants(as_bytes=True):
        #print (variant.genotypes,file=outF1)
        print(variant.position,file=outF2)
 for haplotype in tree_seq.haplotypes():
        print(haplotype,file=outF1)
 outF1.close()
 outF2.close()


#Ni=N0 #for sample size at time of changing simulators (all population)

model0(NBASES,MUT,REC,N0,N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,Nanc,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,SAMPLE)

