#!/usr/bin/env python3
import pandas as pd

import libcoolhext as lc

print("internal fluid => fluid_1")
print("external fluid => fluid_2")  

class Coolhext():
    data = pd.read_csv("input.csv")
    #data.columns # obtain the header names
    print(data)
    
    fgeom = pd.read_csv("fin_geom.csv")
    #fgeom.columns
    print(fgeom)
    
    tgeom = pd.read_csv("tube_geom.csv")
    #tgeom.columns
    print(tgeom)
    
    for index, row in data.iterrows():
        t_1_out = data.t_1_out[index]
        #print(data.fluid_1[index])
        if data.fluid_1[index] == "water":
            print("è acqua")
        if row.fluid_1 == "water": # altro modo
            print("è acqua")
        # mettere tutti i dati in un dict
        if pd.isna(t_1_out): # se è NaN
        #if row.t_1_out|m_dot_1 == "m_dot_1":  t_1_out_or_m_dot1
    	    print(index, "calcolo diretto con portata data")
    	    fluid_1 = libcoolhext.Fluid(row.t_1_in, None, row.m_dot_1, row.fluid_1, row.press_1_bara*1e5) 
    	    fluid_2 = libcoolhext.Fluid(row.t_2_in, None, None, row.fluid_2, row.press_2_bara*1e5)
    	    fluid_2.calcProp(row.t_2_in)
    	    fluid_2.m_dot = row.V_dot_2*fluid_2.rho
    	    hexc = libcoolhext.HeatExchanger(fluid_1, fluid_2, row.arrangement, row.control_vol_per_tube)
    	    hexc.run()
        else:
    	    print(index, "calcolo iterativo calcolando portata")
    	    
    def Help():
         print("## help for Coolhext ##")
         print("control_vol_per_tube:")
         print("\t = 0\t\t simple model - uniform overall heat tranfer coefficient")
         print("\t = [1..n]\t distribuited model - multiple heat tranfer coefficients")
Coolhext()
#Coolhext.Help()
