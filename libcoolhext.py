from math import log, pi, tanh

import CoolProp.CoolProp as CP
import numpy as np

mat_cond = {"copper": 330, "aluminium": 200, "prepainted_aluminium": 170} # W/(m K)
def DegToKelvin(t):
       return t + 273.15     

class Fluid():
    def __init__(self, t_i, t_o, m_dot, fluid, pressure):
        self.t_i = t_i 		# inlet temp
        self.t_o = t_o		# outlet temp
        self.m_dot = m_dot
        self.fluid = fluid        
        self.pressure = pressure

    def calcProp(self, t):
	# specific heat
        self.cp = CP.PropsSI('CPMASS','T',DegToKelvin(t),'P',self.pressure,self.fluid)
	# density
        self.rho = CP.PropsSI('D','T',DegToKelvin(t),'P',self.pressure,self.fluid)
	# thermal conductivity
        self.k =  CP.PropsSI('L','T',DegToKelvin(t),'P',self.pressure,self.fluid)	
	# viscosity
        self.nu =  CP.PropsSI('V','T',DegToKelvin(t),'P',self.pressure,self.fluid)
        # Prandtl
        self.Pr = self.cp*self.nu/self.k

class Tube():
    def __init__(self, d_ext, thickness, mat):
        self.d_ext = d_ext
        self.d_int = d_ext - 2*thickness
        self.thickness = thickness
        self.k = mat_cond[mat]
        self.A_ext_lin = pi*self.d_ext
        self.A_int_lin = pi*self.d_int
        self.A_cross_lin = 0.25*pi*self.d_int**2
        
    def HeatTransferCoefficientConvection(self):
        pass

    def HeatTransferCoefficientCondensation(self):
        pass
        
    def PressureDropConvection(self):
        pass     

    def PressureDropCondensation(self):
        pass

class Fin():
    def __init__(self, tube_pitch, row_pitch, fin_pitch, thickness, mat, tube_diam):
        self.tube_pitch = tube_pitch
        self.row_pitch = row_pitch
        self.collar_height = fin_pitch - thickness
        self.collar_diam = tube_diam + 2*thickness
        self.collar_area = pi*self.collar_diam*self.collar_height
        self.thickness = thickness
        self.k = mat_cond[mat]
        deq = (4 * tube_pitch * row_pitch / pi)**0.5 # equivalent fin diameter
        h = (deq - self.collar_diam)*0.5
        self.he = h * (1 + 0.35 * log(deq / self.collar_diam))
        

    def FinEfficiency(self, alpha, A_fin, A_ext): # A_fin A_ext  o meglio il loro rapporto si possono già definire su FIN (non cambia se è una sola aletta o un pacco alettato)
        m = (2 * alpha / (self.k_fin * self.thickness))**0.5
        eta = tanh(m * self.he) / (m * self.he)
        eta_avg = 1 - A_fin / A_ext * (1 - eta)
        return eta_avg
        
    def HeatTransferCoefficientConvection(self):
        pass

    def PressureDropConvection(self):
        pass        

class HexBasics():
    def __init__(self, fluid_1, fluid_2):
        self.fluid_1 = Fluid(fluid_1.t_i, fluid_1.t_o, fluid_1.m_dot, fluid_1.fluid, fluid_1.pressure)   ## creare qui l'istanza del fluido della cella ricordarsi di dividere per la cella, oppure prima di passarla alla cella (meglio)
        self.fluid_2 = Fluid(fluid_2.t_i, fluid_2.t_o, fluid_2.m_dot, fluid_2.fluid, fluid_2.pressure) ## creare qui l'istanza del fluido della cella
        self.t_1_in = fluid_1.t_i
        self.t_2_in = fluid_2.t_i
        self.t_1_out = None
        self.t_2_out = None
        self.U = None
        self.A_int = None # calcolarli da tube
        self.A_ext = None # calcolarlo da fin
        #self.fluid_1.calcProp(self.t_1_in) ## forse non serve?
        #self.fluid_2.calcProp(self.t_2_in) ## forse non serve?    
        
class HexCell(HexBasics):
    def __init__(self, fluid_1, fluid_2):
        HexBasics.__init__(self, fluid_1, fluid_2)
        
    def dtml(self):
        dtml = abs((self.t_1_in - self.t_2_out) - (self.t_1_out - self.t_2_in))/\
        	log((self.t_1_in - self.t_2_out)/(self.t_1_out - self.t_2_in))
        return dtml
        
    def NTU2(self):
        C_2 = fluid_2.cp*fluid_2.m_dot
        NTU2 = self.U*self.A_ext/C_2
        return NTU2
    
    def R(self):
        C_1 = fluid_1.cp*fluid_1.m_dot
        C_2 = fluid_2.cp*fluid_2.m_dot
        R = C_2/C_1
        return R
        
    def OverallHTC(self, alpha_int, alpha_ext, R_tube, R_fincollar):
        # TO DO
        U_int = 1 / (R_tube + R_fincollar + 1 / alpha_int + self.A_int / (self.A_ext * alpha_ext))
        U_ext = U_int*self.A_int/self.A_ext
        return U_ext     
        
    def P(self, arrangement, rows):
        # P refferred to fluid 2, unmixed (fin side)
        # R = C2/C1
        
        NTU2 = self.NTU2()
        R = self.R()

        if arrangement == "countercurrent" or \
          (arrangement == "crossflow-countercurrent" and rows > 4): 
            # counter flow
            P = (1-exp(-NTU2*(1-R)))/(1-R*exp(-NTU2*(1-R)))

        if arrangement == "cocurrent":
            # parallel flow (cocurrent flow)
            P = (1-exp(-(R+1)*NTU2))/(R+1)

        if (arrangement == "crossflow" or arrangement == "crossflow-countercurrent")\
            and rows == 1:
            # single pass, one row cross flow. Fluid 2 unmixed; Fluid 1 Mixed 
            # HEDH 1983; ESDU 1991
            P = 1/R*(1-exp(-R*(1-exp(-NTU2))))

        K = 1 - exp(-NTU2/rows)
        if arrangement == "crossflow": 
            # HEDH 1983 1.5.3-8. ESDU 1991. Fluid 2 unmixed; Fluid 1 Mixed 
            if rows == 2:
                # single pass, 2 rows   
                P = 1/R*(1-exp(-2*K*R)*(1+R*K**2))
            if rows == 3:
                # single pass, 3 rows
                P = 1/R*(1-(exp(3*K*R)/(1+R*K**2*(3-K)+3/2*R**2*K**4))**-1)
            if rows == 4:
                # single pass, 4 rows
                P = 1/R*(1-(exp(4*K*R)/(1+R*K**2*(6-4*K+K**2)+\
                    4*R**2*K**4*(2-K)+8/3*R**3*K**6))**-1)
            if rows > 4:
                # unmixed/unmixed cross flow
                # HEDH 1983 1.3.1-4 (modified sign NTU**-0.22); ESDU 1991
                P = 1- exp(NTU2**0.22*(exp(-R*NTU2**0.78)-1)/R)

        if arrangement == "crossflow-countercurrent" and (rows > 1 and rows <= 4) :
            # HEDH 1983 1.5.3-8. Fluid 2 unmixed
            if rows == 2: # 2 rows, 2 passes
                AA = K/2 + (1-K/2)*exp(2*K*R)

            if rows == 3: # 3 rows, 3 passes
                AA = K*(1-K/4-R*K*(1-K/2))*exp(K*R)+(1-K/2)**2*exp(3*K*R)

            if rows == 4: # 4 rows, 4 passes
                ## ESDU 1983 non è corretta
                # AA = K/2*(1-K/1+1/4*K**2) + K*(1-K/2)*(1-R/8*K*(1-K/2)*exp(2*K*R))+(1-K/2)**3*exp(4*K*R)
                ## VDI HEAT ATLAS
                AA = K/2*(1-K/2+1/4*K**2) + K*(1-K/2)*(1-2*K*R*(1-K/2))*exp(2*K*R)+(1-K/2)**3*exp(4*K*R)

            P = 1/R*(1-1/AA)
        return P  
        
    def EnergyBalance(self):
       t_1_avg = 0.5*(self.fluid_1.t_i + self.fluid_1.t_o) # l'istanza del fluido deve essere quella relativa alla cella e non a quella globale
       t_2_avg = 0.5*(self.fluid_2.t_i + self.fluid_2.t_o)
       self.fluid_1.calcProp(t_1_avg)
       self.fluid_2.calcProp(t_2_avg)
       alpha_int = self.tube.HeatTransferCoefficientConvection()
       alpha_ext = self.fin.HeatTransferCoefficientConvection()
       eta_avg = FinEfficiency(alpha_ext) # sistemare aree int ext
       alpha_ext_w =  alpha_ext * eta_avg)      
       self.U = OverallHTC(alpha_int, alpha_ext_w, R_tube, R_fincollar)
       P = self.P(arrangement, row_num)       
       theta =  P / self.NTU2()
       dtmax = abs(self.fluid_1.t_i - self.fluid_2.t_i)
       dtm = dtmax*theta
       self.q = self.U * self.A_ext * dtm
       self.fluid_2.t_o = self.fluid_2.t_i + P * (self.fluid_1.t_i - self.fluid_2.t_i)
       weight_newvalue = 0.95
       self.fluid_1.t_o = (self.fluid_1.t_i - self.R() * (self.fluid_2.t_o - self.fluid_2.t_i))*weight_newvalue + self.fluid_1.t_o*(1 - weight_newvalue)
              
        

class HeatExchanger(HexBasics): # ma serve fare inheritance da HexBasic o basta solo HexCell??
    def __init__(self, fluid_1, fluid_2, arrangement, control_vol_per_tube):   
        HexBasics.__init__(self, fluid_1, fluid_2)  # ma serve HexBasic o basta solo HexCell??
        self.control_vol_per_tube = control_vol_per_tube
        
        self.tube = Tube() # aggiungere proprietà tube, oppure crearle direttamente sul main e passarle alla classe (meglio), come il fluido
        self.fin = Fin()
        
    def run(self):
        if self.control_vol_per_tube == 0:
            # simple model
            i_max = 1 # num of rows
            j_max = 1 # num of tubes per row
            k_max = 1 # control volumes per tube
        else:
            #i_max = self.row_num # num of rows
            #j_max = self.tubePerRow_num # num of tubes per row
            k_max = self.control_vol_per_tube # control volumes per tube
            i_max = 4
            j_max = 1

        hxCells = []
        hxCells_j = []
        hxCells_k = []
                
        for k in range(0, k_max):
            hxCells_k.append(HexCell(self.fluid_1, self.fluid_2))  # l'istanza del fluido deve essere quella relativa alla cella e non a quella globale, modificare
        for j in range(0, j_max):
            hxCells_j.append(hxCells_k)
        for i in range(0, i_max):
            hxCells.append(hxCells_j)
            
        #print(hxCells)
        err_t = 1
        counter = 0
                        
        while err_t < 1e-6: ## ciclo temporale // sarebbe meglio un errore sulle temperature uscita globali, da inventarsi es. con numpy)
            for i in range(0, i_max):
                for j in range(0, j_max):
                    for k in range(0, k_max):
                        print("ora qui va fatto il calcolo") ## inizializzare prima
                        hxCells[i][j][k].EnergyBalance()
           ttfo = hxCells[i][j][k_max-1]    # inserire anche mod 2          
           err_t = abs(ttfo - self.fluid_1.t_o)
           counter += 1  
            

