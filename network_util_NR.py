import numpy as np 

# np.set_printoptions(precision=3, suppress=True)

# To-do's: 

# - To add: nonlinear update solver 

# - Test if connecting block definition can replace explicit definition of connect_list
# - Implement RC, CR and RCR blocks
# - Implement Migliavacca's code
# - To add graphical output routines
# - Test : inductances, diodes, time or volume varying capacitors (needs function input)
# - Implement input output format 
# - Add code to allow conversion in standard elements
# - Make code consistent with cplbc
# - Look up Sam's GUI 

# Programming tip - don't use mutable default arguments in classes!


# ------------ General LPN definitions --------------------
#==========================================================
#						LPN VARIABLES
#==========================================================

class LPNVariable:
	def __init__(self,value,units,name="NoName",vtype='ArbitraryVariable'):
		self.type = vtype
		self.value = value
		# Two generic units accepted : SI, cgs. Conversion for special values applied
		self.units = units  
		self.name = name 


class PressureVariable(LPNVariable): 

	def __init__(self,value,units='cgs',name='NoNamePressure'): 
		LPNVariable.__init__(self,value=value,units=units,name=name,vtype='Pressure')

	def convert_to_cgs(self): 
		if self.units == 'cgs':
			print "Variable: "+self.name+" already at cgs" 
		elif self.units == 'SI' : 
			self.value = self.value*1.0E5
			self.units = 'cgs'
		elif self.units == 'mmHg': 
			self.value = self.value*0.001333224
			self.units = 'cgs'
		else : 
			raise Exception("Units "+self.units+" not recognized")


	def convert_to_mmHg(self): 
		if self.units == 'cgs':
			self.value = self.value*750.06
			self.units = 'mmHg'
		elif self.units == 'SI' : 
			self.value = self.value*7.50*1E-3
			self.units = 'mmHg'
		elif self.units == 'mmHg': 
			print "Variable: "+self.name+" already at mmHg" 
		else : 
			raise Exception("Units "+self.units+" not recognized")


class FlowVariable(LPNVariable): 

	def __init__(self,value,units='cgs',name='NoNameFlow'): 
		LPNVariable.__init__(self,value=value,units=units,name=name,vtype='Flow')

	def convert_to_cgs(self): 
		if self.units == 'cgs':
			print "Variable: "+self.name+" already at cgs" 
		elif self.units == 'SI' : 
			self.value = self.value*1.0E-6
			self.units = 'cgs'
		elif self.units == 'Lpm' :  # litres per minute 
			self.value = self.value*16.6667
			self.units = 'cgs'			
		else : 
			raise Exception("Units "+self.units+" not recognized")


	def convert_to_Lpm(self): 
		if self.units == 'cgs':
			self.value = self.value/16.6667
			self.units = 'Lpm'
		elif self.units == 'SI' : 
			self.value = self.value/(16.6667*1.0E-6)
			self.units = 'Lpm'
		elif self.units == 'Lpm': 
			print "Variable: "+self.name+" already at Lpm" 
		else : 
			raise Exception("Units "+self.units+" not recognized")


#==========================================================
#						LPN WIRES
#==========================================================

# Wires connect circuit elements and junctions 
# They can only posses a single pressure and flow value (system variables)
# They can also only possess one element(or junction) at each end
class wire:
	def __init__(self,connecting_elements,Pval=0,Qval=0,name="NoNameWire",P_units='cgs',Q_units='cgs'): 
		self.name=name
		self.type='Wire'
		self.P = PressureVariable(value=Pval,units=P_units,name=name+"_P")
		self.Q = FlowVariable(value=Qval,units=Q_units,name=name+"_Q")
		if len(connecting_elements) >2 : 
			raise Exception('Wire cannot connect to more than two elements at a time. Use a junction LPN block')
		if type(connecting_elements)!=tuple : 
			raise Exception('Connecting elements to wire should be passed as a 2-tuple')
		self.connecting_elements = connecting_elements

		self.LPN_solution_ids = [None]*2  

	def initialize_PQ(self,y,P,Q):
		y[self.LPN_solution_ids[0]] = P 
		y[self.LPN_solution_ids[0]] = Q 


#==========================================================
#						LPN BLOCKS
#==========================================================


# -- GENERAL LPN BLOCK

class LPNBlock:
	def __init__(self,connecting_block_list=None,name = "NoName",flow_directions=None):
		if connecting_block_list == None: 
			connecting_block_list = []
		self.connecting_block_list = connecting_block_list
		self.num_connections = len(connecting_block_list)
		self.name = name
		self.neq = 2
		self.type="ArbitraryBlock"
		self.num_block_vars = 0
		self.connecting_wires_list = []
		if flow_directions == None : 
			self.flow_directions = [] # -1 : Inflow to block, +1 outflow from block 
		else : 
			self.flow_directions = flow_directions
		self.LPN_solution_ids = [] 
		self.emxcoe = []
		self.fmxcoe = []

		# Tangent matrix coes 
		self.demxcoe = [] # sum_k[ydot_k(d(E_ik)/dy_j)]
		self.dfmxcoe = [] # sum_k[y_k(d(F_ik)/dy_j)]
		self.dcmxcoe = [] # (d(C_i)/dy_j)

		self.eq_id_list = []

	def check_block_consistency(self): 
		return 

	def add_connecting_block(self,block,direction):
		# Direction = +1 if flow sent to block
		#			= -1 if flow recvd from block
		self.connecting_block_list.append(block)
		self.num_connections = len(self.connecting_block_list)
		self.flow_directions.append(direction)
		# print self.name


	def add_connecting_wire(self,new_wire):
		self.connecting_wires_list.append(new_wire)
		# self.flow_directions.append(direction)

	def local_eq_coe_def(self,args): 
		return 

	def eqids(self,wire_dict,local_eq): 
		# EqID returns variable's location in solution vector
		nwirevars = self.num_connections*2
		if local_eq < nwirevars: 
			vtype = local_eq%2 # 0 --> P, 1 --> Q
			wnum = int(local_eq/2)
			return wire_dict[self.connecting_wires_list[wnum]].LPN_solution_ids[vtype]
		else : 
			vnum = local_eq - nwirevars
			return self.LPN_solution_ids[vnum]

# -- JUNCTION
# Junction points between LPN blocks with specified directions of flow
class Junction(LPNBlock):
	def __init__(self,connecting_block_list=None,name="NoNameJunction",flow_directions=None):
		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
		self.type = "Junction"
		self.neq = self.num_connections


	def add_connecting_block(self,block,direction):
		self.connecting_block_list.append(block)
		self.num_connections = len(self.connecting_block_list) 
		self.neq = self.num_connections
		self.flow_directions.append(direction)
		# print self.name

	def local_eq_coe_def(self,args): 

		# Number of variables per tuple = 2*num_connections
		# Number of equations = num_connections-1 Pressure equations, 1 flow equation
		# Format : P1,Q1,P2,Q2,P3,Q3, .., Pn,Qm

		self.emxcoe = [(0,)*(2*self.num_connections)]*(self.num_connections)

		self.fmxcoe = [ (1.,)+(0,)*(2*i+1) + (-1,) + (0,)*(2*self.num_connections-2*i-3) for i in range(self.num_connections-1) ]

		tmp = (0,)
		for d in self.flow_directions[:-1]: 
			tmp+=(d,)
			tmp+=(0,)

		tmp += (self.flow_directions[-1],)
		self.fmxcoe.append(tmp)
		self.cveccoe = [0]*self.num_connections

		self.demxcoe = [(0,)*(2*self.num_connections)]*(self.num_connections)
		self.dfmxcoe = [(0,)*(2*self.num_connections)]*(self.num_connections)
		self.dcmxcoe = [(0,)*(2*self.num_connections)]*(self.num_connections)

# -- Resistance
class Resistance(LPNBlock):
	def __init__(self,R,connecting_block_list=None,name="NoNameResistance",flow_directions=None):
		
		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
		self.type = "Resistance"
		self.R = R

	def check_block_consistency(self): 
		if len(connecting_block_list) != 2:
			raise Exception("Resistance block can be connected only to two elements")

	def local_eq_coe_def(self,args): 

		# For resistors, the ordering is : (P_in,Q_in,P_out,Q_out)

		self.emxcoe = [(0,)*4]*2
		self.fmxcoe = [(1.,-1.*self.R,-1.,0),(0,1.,0,-1.)]
		self.cveccoe = [0]*2

		self.demxcoe = [(0,)*4]*2
		self.dfmxcoe = [(0,)*4]*2
		self.dcmxcoe = [(0,)*4]*2

# -- Flow dependent Resistance : delta_P = q*Rfunc(t,q)
class FlowDepResistance(LPNBlock):
	def __init__(self,Rfunc,connecting_block_list=None,name="NoNameFlowDepResistance",flow_directions=None):
		
		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
		self.type = "FlowDepResistance"
		self.Rfunc = Rfunc

	def check_block_consistency(self): 
		if len(connecting_block_list) != 2:
			raise Exception("FlowDepResistance block can be connected only to two elements")

	def linearize_func(self,t,q,func,eps_dq=1e-3): 

		# Input  : func(t,q), typically t = t_n+1 ie t_n + dt
		# Return : d/dq (func) at (t,q) using central differences

		dq = eps_dq*q 
		if abs(dq) < 1e-5:
			dq = 1e-5
		return (func(t,q+dq)-func(t,q-dq))/(2.*dq)

	def local_eq_coe_def(self,args): 
		# delta_P = q*Rfunc(t,q) , so linearization yields: 
		# delta_P(n+1) - q(n+1)*[Rfunc(t(n+1),q(n))+q(n)*dR] + q(n)^2*dR = 0

		t = args['Time']
		curr_y = args['Solution']
		# dt = args['Time step']
		# rho = args['rho']
		# alpha_f = 1/(1.0+rho)
		wire_dict = args['Wire dictionary']
		q = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1]]	

		# pi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[0]]	
		# po = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[0]]	

		dR = self.linearize_func(t,q,self.Rfunc)
		R = self.Rfunc(t,q)

		# c_lin = q*q*dR
		# self.emxcoe = [(0,)*4]*2
		# self.fmxcoe = [(1.,-1.*(self.Rfunc(t+dt*alpha_f,q)+q*dR),-1.,0),(0,1.,0,-1.)]
		# self.cveccoe = [c_lin,0]

		self.emxcoe = [(0,)*4]*2
		# self.fmxcoe = [(1./(R),-1.,-1./(R),0),(0,1.,0,-1.)]
		self.fmxcoe = [(1.,-1.*(R),-1.,0),(0,1.,0,-1.)]
		self.cveccoe = [0,0]

		self.demxcoe = [(0,)*4]*2
		# self.dfmxcoe = [(0,(pi-po)*(-dR/(R*R)),0,0),(0,)*4]
		self.dfmxcoe = [(0,-dR*q,0,0),(0,)*4]
		self.dcmxcoe = [(0,)*4]*2


# -- Pressure reference
class PressureRef(LPNBlock):
	def __init__(self,Pref,connecting_block_list=None,name="NoNamePressureRef",flow_directions=None):
		
		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
		self.type = "PressureRef"
		self.neq = 1
		self.Pref=Pref

	def check_block_consistency(self): 
		if len(connecting_block_list) != 1:
			raise Exception("PressureRef block can be connected only to one element")

	def local_eq_coe_def(self,args): 
		self.emxcoe = [(0,)]
		self.fmxcoe = [(1.,)]
		self.cveccoe = [-1.0*self.Pref]

		self.demxcoe = [(0,)]
		self.dfmxcoe = [(0,)]
		self.dcmxcoe = [(0,)]

# -- Unsteady P reference
class UnsteadyPressureRef(LPNBlock):
	def __init__(self,Pfunc,connecting_block_list=None,name="NoNamePressureRef",flow_directions=None):
	 
	 LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
	 self.type = "PressureRef"
	 self.neq = 1
	 self.Pfunc=Pfunc

	def check_block_consistency(self): 
	 if len(connecting_block_list) != 1:
	  raise Exception("PressureRef block can be connected only to one element")

	def local_eq_coe_def(self,args): 
	 t = args['Time']
	 self.emxcoe = [(0,0)]
	 self.fmxcoe = [(1.,0.)]
	 self.cveccoe = [-1.0*self.Pfunc(t)]

	 self.demxcoe = [(0,)]
	 self.dfmxcoe = [(0,)]
	 self.dcmxcoe = [(0,)]

# -- Flow reference
class UnsteadyFlowRef(LPNBlock):
	def __init__(self,Qfunc,connecting_block_list=None,name="NoNameFlowRef",flow_directions=None):
		
		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
		self.type = "FlowRef"
		self.neq = 1
		self.Qfunc=Qfunc

	def check_block_consistency(self): 
		if len(connecting_block_list) != 1:
			raise Exception("FlowRef block can be connected only to one element")

	def local_eq_coe_def(self,args): 
		t = args['Time']
		self.emxcoe = [(0,0)]
		self.fmxcoe = [(0,1.)]
		self.cveccoe = [-1.0*self.Qfunc(t)]

		self.demxcoe = [(0,)]
		self.dfmxcoe = [(0,)]
		self.dcmxcoe = [(0,)]

# -- Capacitance
class Capacitance(LPNBlock):
	def __init__(self,C,connecting_block_list=None,name="NoNameCapacitance",flow_directions=None):
		
		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
		self.type = "Capacitance"
		self.C = C

	def check_block_consistency(self): 
		if len(connecting_block_list) != 2:
			raise Exception("Capacitance block can be connected only to two elements")


	def local_eq_coe_def(self,args): 
		self.emxcoe = [(1.0*self.C,0,-1.0*self.C,0),(0,0,0,0)]
		self.fmxcoe = [(0,-1.0,0,0),(0,1.,0,-1.)]
		self.cveccoe = [0,0]

		self.demxcoe = [(0,)*4]*2
		self.dfmxcoe = [(0,)*4]*2
		self.dcmxcoe = [(0,)*4]*2



# -- RCL - constant resistor, capacitor, inductor - vessel representation
# Formulation includes additional variable : internal pressure proximal to capacitance.
# Trading off speed with enforced uniqueness
class RCLBlock(LPNBlock):
	def __init__(self,R,C,L,connecting_block_list=None,name="NoNameRCL",flow_directions=None):
		
		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
		self.type = "RCL"
		self.neq = 3
		self.num_block_vars = 1
		self.R = R
		self.C = C
		self.L = L

	def check_block_consistency(self): 
		if len(connecting_block_list) != 2:
			raise Exception("RCL block can be connected only to two elements")

	def local_eq_coe_def(self,args): 
		self.emxcoe = [(0,0,0,-self.L,0),      (0,0,0,0,-self.C),   (0,0,0,0,0)]
		self.fmxcoe = [(1.,-self.R,-1.,0,0), (0,1.,0,-1.,0),      (1.,-self.R,0,0,-1.)]
		self.cveccoe = [0,0,0]

		self.demxcoe = [(0,)*5]*3
		self.dfmxcoe = [(0,)*5]*3
		self.dcmxcoe = [(0,)*5]*3


# -- RC - constant resistor, capacitor - low inertia vessel 
class RCBlock(LPNBlock):
	def __init__(self,R,C,connecting_block_list=None,name="NoNameRC",flow_directions=None):
		
		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
		self.type = "RC"
		self.R = R
		self.C = C

	def check_block_consistency(self): 
		if len(connecting_block_list) != 2:
			raise Exception("RC block can be connected only to two elements")

	def local_eq_coe_def(self,args): 
		self.emxcoe = [(0,0,0,0),      (0,0,-self.C,0)]
		self.fmxcoe = [(1.0,-self.R,-1.0,0), (0,1.,0,-1.)]
		self.cveccoe = [0,0]

		self.demxcoe = [(0,)*4]*2
		self.dfmxcoe = [(0,)*4]*2
		self.dcmxcoe = [(0,)*4]*2


# -- RCR - constant RCR - outflow representation
# Formulation includes additional variable : internal pressure proximal to capacitance.
# Trading off speed with enforced uniqueness
class RCRBlock(LPNBlock):
	def __init__(self,Rp,C,Rd,connecting_block_list=None,name="NoNameRCR",flow_directions=None):
		
		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
		self.type = "RCR"
		self.neq = 3
		self.num_block_vars = 1
		self.Rp = Rp
		self.C = C
		self.Rd = Rd

	def check_block_consistency(self): 
		if len(connecting_block_list) != 2:
			raise Exception("RCR block can be connected only to two elements")

	def local_eq_coe_def(self,args): 
		self.emxcoe = [(0,0,0,0,0),      (0,0,0,0,-self.C),   (0,0,0,0,0)]
		self.fmxcoe = [(1.0,-self.Rp,-1.0,-self.Rd,0), (0,1.,0,-1.,0),      (1.,-self.Rp,0,0,-1.)]
		self.cveccoe = [0,0,0]

		self.demxcoe = [(0,)*5]*3
		self.dfmxcoe = [(0,)*5]*3
		self.dcmxcoe = [(0,)*5]*3



# -- Time Varying Capacitance
class TimeDependentCapacitance(LPNBlock):
	def __init__(self,Cfunc,connecting_block_list=None,name="NoNameTimeDependentCapacitance",flow_directions=None):
		
		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
		self.type = "TimeDependentCapacitance"
		self.Cfunc = Cfunc

	def check_block_consistency(self): 
		if len(connecting_block_list) != 2:
			raise Exception("Capacitance block can be connected only to two elements")


	def local_eq_coe_def(self,args): 
		t = args['Time']
		# dt = args['Time step']
		# rho = args['rho']
		# alpha_f = 1/(1.0+rho)		
		self.emxcoe = [(1.0*self.Cfunc(t),0,-1.0*self.Cfunc(t),0),(0,0,0,0)]
		self.fmxcoe = [(0,-1.0,0,0),(0,1.,0,-1.)]
		self.cveccoe = [0,0]

		self.demxcoe = [(0,)*4]*2
		self.dfmxcoe = [(0,)*4]*2
		self.dcmxcoe = [(0,)*4]*2


# -- Chamber Capacitance -- with direct prescription of pressure
class ChamberModel(LPNBlock):
	def __init__(self,Activefunc,Passivefunc,Pref,connecting_block_list=None,name="NoNameChamber",flow_directions=None):

		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)

		self.neq = 3
		self.type = "Chamber"
		self.num_block_vars = 1 # Chamber volume
		self.Activefunc = Activefunc
		self.Passivefunc = Passivefunc
		self.Pref = Pref


	def linearize_func(self,t,v,func,eps_dv=1e-4): 

		# Input  : func(t,v), typically t = t_n+1 ie t_n + dt
		# Return : d/dv (func) at (t,v) using central differences

		dv = eps_dv*v 
		return (func(t,v+dv)-func(t,v-dv))/(2.*dv)

	def linearize_func_t(self,t,v,func,dt=1e-3): 

		return (func(t+dt,v)-func(t-dt,v))/(2.*dt)

	def check_block_consistency(self): 
		if len(connecting_block_list) != 2:
			raise Exception("Chamber block can be connected only to two elements")
		return

	def initialize_volume(self,sol_vec,Vu=1.0,scale_fact=1.5):

		sol_vec[self.LPN_solution_ids[0]] = scale_fact*Vu

	def local_eq_coe_def(self,args): 
		t = args['Time']
		curr_y = args['Solution']
		# dt = args['Time step']
		# rho = args['rho']
		# alpha_f = 1/(1+rho)

		wire_dict = args['Wire dictionary']
		qi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1]]
		qo = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[1]]

		v = curr_y[self.LPN_solution_ids[0]] 

		# print self.name,' has volume: ',v
		if v<=0: 
			print "Time: ",t
			raise Exception("Zero chamber volume detected in chamber: ",self.name)

		a_lin = self.linearize_func(t,v,self.Activefunc)
		p_lin = self.linearize_func(t,v,self.Passivefunc)

		# c_from_lin = -(self.Activefunc(t,v) + self.Passivefunc(t,v)) + v*(a_lin+p_lin) 
		# print 'Volume ',v,'has pressure ',(self.Activefunc(t,v) + self.Passivefunc(t,v))/1333.33

		self.emxcoe = [(0,0,0,0,1.),     (0,0,0,0,0),     (0,0,0,0,0)]
		self.fmxcoe = [(0,-1.,0,1.,0), (-1.,0,1.,0,0),  (1.,0,0.,0.,0.)]
		self.cveccoe = [0,0,-(self.Activefunc(t,v)+self.Passivefunc(t,v))-self.Pref]

		# self.fmxcoe = [(0,-1.,0,1.,0), (-1.,0,1.,0,0),  (1.,0,0.,0,-(a_lin+p_lin))]
		# self.cveccoe = [0,0,c_from_lin+self.Pref]

		self.demxcoe = [(0,)*5,(0,)*5,(0,)*5]
		self.dfmxcoe = [(0,)*5,(0,)*5,(0,)*5]
		self.dcmxcoe = [(0,)*5,(0,)*5,(0,0,0,0,-(a_lin+p_lin))]

class pVADsphere(LPNBlock):
	def __init__(self, volume, Pa, h, Young, Poisson, threshold, connecting_block_list=None,name="NoNameVAD",flow_directions=None):
		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
		self.neq = 3
		self.type = "pVADsphere"
		self.num_block_vars = 1 # Chamber volume
		self.volume = volume 
		self.Pa = Pa
		self.h = h 
		self.Young = Young
		self.Poisson = Poisson 
		self.threshold = threshold 
		self.r = ((self.volume*2.)/((4./3.)*np.pi))**(1./3.)
		self.D = (self.Young*(self.h**3))/(12*(1-self.Poisson**2)) #di solito ~1.335*10**6 
		# self.D = 6.09*10**5
		# self.D = 7.7*10**4 
		print "D: ", self.D
		print "Radius: ", self.r
		print "Limit: ", 64.*self.D/self.r**3
		print "% Threshold: ", (1.-self.threshold)

	def check_block_consistency(self): 
		if len(connecting_block_list) != 2:
			raise Exception("pVAD block can be connected only to two elements")
		return

	def displacement(self, pi, t): 
		# return ((self.Pa-pi)*self.r**4)/(64.*self.D)
		w = ((self.Pa(t)-pi)*self.r**4)/(64.*self.D)
		return w
		# if abs(w) < (1.-self.threshold)*self.r: 
		# 	# print "AmeriCone Dream", w/self.r
		# 	return w 
		# else: 
		# 	# print "Soviet Sphere of Influence", w/self.r
		# 	c = 1./(self.r*(1.+self.threshold))
		# 	# return self.r*(1-self.threshold) + self.r*(1+self.threshold)*np.tanh(c*(w-self.r*(1-self.threshold)))
		# 	return self.r*(1-self.threshold)+self.r*self.threshold*np.tanh(c*(w-self.r*(1-self.threshold)))
		# return ((self.Pa(t)-pi)*self.r**4)/(64.*self.D)

	def displacement_for_plots(self, args): 
		# return ((self.Pa-pi)*self.r**4)/(64.*self.D)
		t = args['Time']
		curr_y = args['Solution']
		wire_dict = args['Wire dictionary']
		pi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[0]]
		# return ((self.Pa(t)-pi)*self.r**4)/(64.*self.D)
		return self.displacement(pi,t)

	def membrane_for_plots(self, args): 
		t = args['Time']
		curr_y = args['Solution']
		wire_dict = args['Wire dictionary']
		pi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[0]]
		a = self.displacement_for_plots(args)
		w_max = self.displacement(pi, t)
		return w_max*(1 - (self.r/a)**2)**2

	def fpi(self, pi, t): 
		return (1./3.)*np.pi*self.r**2*self.displacement(pi,t)
		# return (2./3.)*np.pi*self.r**2*self.displacement(pi, t)

	def linearize_pi(self,pi,t,eps): 
		return (self.fpi(pi+eps,t)-self.fpi(pi-eps,t))/(2*eps)

	def local_eq_coe_def(self,args): 
		t = args['Time']
		curr_y = args['Solution']

		wire_dict = args['Wire dictionary']
		# Previous time step flow in and out
		qi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1]]
		qo = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[1]]
		pi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[0]]
		po = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[0]]

		v = curr_y[self.LPN_solution_ids[0]] 

		# print "Pressure in: ", pi
		if (np.isnan(pi)): 
			raise Exception("Pressure in is NAN type")

		# Misc constants
		w = self.displacement(pi, t)

		# if abs(w) >= (1-self.threshold)*self.r: 
		# # 	print "The AmeriCone Dream", w/self.r
		# # else:
		# 	print "Soviet Sphere of Influence", w/self.r

		# args['Displacement'].append(w)
		eps = max(10**-5, 1e-8*pi)
		pi_lin = self.linearize_pi(pi,t,eps)

		# if abs(pi-self.Pa(t)) >= (64.*self.D)/(0.75*self.r**3): 
		# 	print "Exceeded soft limit"
		# 	max_r = 0.75**0.25*self.r 
		# 	w = (self.Pa(t)-pi)*max_r**4/(64.*self.D)

		# print "Central diff: ", pi_lin
		# cd = (self.fpi(pi+epsilon,t)-self.fpi(pi-epsilon,t))/(2*epsilon)

		# self.emxcoe = [(0.,0.,0.,0.,0.),	(0.,0.,0.,0.,0.),	(0.,0.,0.,0.,1.*self.volume)]
		# self.fmxcoe = [(1.,0.,-1.,0.,0.),	(pi_lin,0.,0.,0.,+1*self.volume),	(0.,-1.,0.,+1.,0.)] #same side of eqn
		# self.cveccoe = [0., -self.volume+self.fpi(pi,t)-pi_lin*pi, 0.]

		#				pi qi po qo V
		self.emxcoe = [(0.,0.,0.,0.,0.),	(0.,0.,0.,0.,0.),	(0.,0.,0.,0.,1.)]
		self.fmxcoe = [(1.,0.,-1.,0.,0.),	(pi_lin,0.,0.,0.,+1.),	(0.,-1.,0.,+1.,0.)] #same side of eqn
		self.cveccoe = [0., -1.*self.volume+self.fpi(pi,t)-pi_lin*pi, 0.]

		self.demxcoe = [(0,)*5,(0,)*5,(0,)*5]
		self.dfmxcoe = [(0,)*5,(0,0,0,0,0),(0,0,0,0,0)]
		self.dcmxcoe = [(0,)*5,(0,)*5,(0,)*5]



		# if (pi > self.Pa(t)): 
		# 	print "Pa < Pi = ", pi, "Displacement: ", w 

		# Annoying print statements to debug
		# if (pi-self.Pa >= (64*self.D/self.r**3)): 
		# if v < self.volume: 
		# 	print "Incomplete ejection at: ", pi 
		# if abs(pi-self.Pa(t)) >= (64.*self.D/self.r**3):
		# 	print 'Limit hit'
		# 	self.emxcoe[2] = (0.,0.,0.,0.,0.)

		# r = ((self.volume*2)/((4/3)*np.pi))**(1./3.)
		# D = (self.Young*(self.h**3))/(12*(1-self.Poisson**2))
		# C1 = (128/315)*(np.pi/((64*D)**2))
		# C2 = w**9 
		# C3 = 9*w**8
		# Pa_lin = C1*C2*(2*pi-2*self.Pa)
		# w_lin = C3*(self.Pa*2+pi*(pi-2*self.Pa))

		# self.emxcoe = [(0.,0.,0.,0.,0.),	(0.,0.,0.,0.,0.),	(0.,0.,0.,0.,1.)]
		# self.fmxcoe = [(1.,0.,-1.,0.,0.), (C1*C2*(pi-2*self.Pa),0.,0.,0.,-1),	(0.,1.,0.,-1.,0.)]
		# self.cveccoe = [0., Pa_lin+w_lin-self.volume+C1*C2*self.Pa**2, 0.]

		# if (pi-po >= (64*D/r**3)): 
		# 	self.emxcoe[2] = (0.,0.,0.,0.,0.)

		# self.demxcoe = [(0,)*5,(0,)*5,(0,)*5]
		# self.dfmxcoe = [(0,)*5,(0,0,0,0,0),(0,0,0,0,0)]
		# self.dcmxcoe = [(0,)*5,(0,)*5,(0,)*5]

		# VICTORIA IS BAD @ BASIC ALGEBRA 
		# v_coeff = (128*np.pi*w**9)/(315*(64*D)**2)
		# C1 = np.pi/((64*D)**2)
		# C2 = 128/315
		# C3 = 9
		# lin_1 = C1*C2*w**9*self.Pa**2
		# lin_2 = C1*(C2*w**9*(-2*self.Pa+2*pi))
		# lin_3 = C1*(C3*w**8*self.Pa**2-C3*w**8*2*self.Pa*pi+C3*w**8*pi**2)
		# pi_coeff = (C1*C2*w**9*2*self.Pa - C2*C1*w**9*pi)

		# self.emxcoe = [(0.,0.,0.,0.,0.),	(0.,0.,0.,0.,0.),	(0.,0.,0.,0.,1.)]
		# self.fmxcoe = [(1.,0.,-1.,0.,0.),	(pi_coeff,0.,0.,0.,-1.), (0.,1.,0.,-1.,0.)]	
		# self.cveccoe = [0., -self.volume-v-lin_1-lin_2-lin_3, 0.]

		# if (pi-po >= (64*D/r**3)):
		# 	self.emxcoe[2] = (0.,0.,0.,0.,0.)

		# rho = 1060. #density kg/m^3
		# mu = 0.004 #viscosity kg/m*s 

		# 			  pi qi po qo  V
		# self.emxcoe = [(0.,0.,0.,0.,0.),	(0.,0.,0.,0.,0.),	(0.,0.,0.,0.,1.)]
		# self.fmxcoe = [(1.,0.,-1.,0.,0.),	(v_coeff*pi,0.,0.,0.,-1.),	(0.,1.,0.,-1.,0.)]
		# self.cveccoe = [0., self.volume-v_coeff*self.Pa**2, 0.]	

		# self.emxcoe = [(0.,0.,0.,0.,0.,0.),				(0.,0.,0.,0.,0.,0.), 	(0.,0.,0.,0.,0.,0.), 				(0.,0.,0.,0.,0.,1.)]
		# self.fmxcoe = [(-r**4/(64*D),0.,0.,0.,-1.,0.),	(1.,0.,-1.,0.,0.,0.),	(w**9*pi*v_coeff,0.,0.,0.,v_coeff,-1.),	(0.,1.,0.,-1.,0.,0.)]
		# self.cveccoe = [(self.Pa(t)*r**4)/(64*D), 0., self.volume-w**9*v_coeff*(self.Pa(t)**2), 0.]


# -- Mynard Valve
class MynardValve(LPNBlock):
	def __init__(self,kv,A_ann,connecting_block_list=None,name="NoNameIdealDiode",flow_directions=None):
		
		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
		self.type = "MynardValve"
		self.neq=3
		self.num_block_vars = 1 # State : 0- shunt, 1- interrupt
		self.kv = kv
		self.A_ann = A_ann
		# self.beta = 0.5*1.06/((np.pi*(0.9/2.)**2)**2)
		self.beta = 0.5*1.06/(self.A_ann**2)	# beta = 0.5*rho/(A_ann^2)
		self.eps = 1e-7
		self.f = 1.0

	def check_block_consistency(self): 
		if len(connecting_block_list) != 2:
			raise Exception("MynardValve block can be connected only to two elements")

	def valve_resistance(self,qi,state):
		return (self.beta*qi*abs(qi))/(state**2 + self.eps)

	def d_res_d_q(self,qi,state,eps=1e-3):
		dq = max(abs(qi)*eps,eps**2)
		return (self.valve_resistance(qi+dq,state)-self.valve_resistance(qi-dq,state))/(2.*dq)

	# Use forward differences here only
	def d_res_d_s(self,qi,state,eps=1e-3):
		ds = max(abs(state)*eps,eps**2)
		return (self.valve_resistance(qi,state+ds)-self.valve_resistance(qi,state))/(ds)


	def state_rhs(self,pi,po,state):
		val = self.kv*(1-state)*(pi-po) if pi >= po else self.f*self.kv*(state)*(pi-po) 
		return val

	def local_eq_coe_def(self,args): 
		# Needs to take current solution vector as argument
		t = args['Time']
		curr_y = args['Solution']
		wire_dict = args['Wire dictionary']
		rho = args['rho']
		dt = args['Time step']
		alpha_f = 1.0/(1.0+rho)

		pi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[0]]
		po = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[0]]		
		qi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1]]		
		# Qo = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[1]]				
		# state = max(curr_y[self.LPN_solution_ids[0]],self.eps)
		state = curr_y[self.LPN_solution_ids[0]]

		# Three equations : flow continuity, bernoulli type resistance and state equation
		# (1) --> as is
		# (2) --> pi-po = self.valve_resistance(qi,state)
		# (3) --> dstate/dt = state_rhs(pi,po,state)

		# Constants from linearization of (2)
		c_lin_res = - self.valve_resistance(qi,state) + state*self.d_res_d_s(qi,state) + qi*self.d_res_d_q(qi,state)

		# Constants from linearization of (3)
		rhs_pi = self.kv*(1-state) if pi >= po else self.f*self.kv*state 
		rhs_po = -rhs_pi
		rhs_state = - self.kv*(pi-po) if pi >= po else self.f*self.kv*(pi-po)
		c_lin_state = -self.state_rhs(pi,po,state) + pi*rhs_pi + po*rhs_po + state*rhs_state


		self.emxcoe = [(0,)*5,(0,)*5,(0,0,0,0,1.)]
		self.fmxcoe = [(0,1,0,-1,0),(1,-self.d_res_d_q(qi,state),-1,0,-self.d_res_d_s(qi,state)),(-rhs_pi,0,-rhs_po,0,-rhs_state)]
		self.cveccoe = [0,c_lin_res,c_lin_state]

		self.demxcoe = [(0,)*5,(0,)*5,(0,)*5]
		self.dfmxcoe = [(0,)*5,(0,)*5,(0,)*5]
		self.dcmxcoe = [(0,)*5,(0,)*5,(0,)*5]		



# -- VAD -- assuming V_Berlin = 25 mL; will build in input/if's later 
class PulsatileVAD(LPNBlock):
	def __init__(self, volume, Pcomp, connecting_block_list=None,name="NoNameVAD",flow_directions=None):

		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)

		self.neq = 3
		self.type = "PulsatileVAD"
		self.num_block_vars = 1 # Chamber volume
		self.volume = volume 
		self.Pcomp = Pcomp
		# if (self.volume == 25. or self.volume == 50.): 
		# 	self.d_can = 12.
		# elif (self.volume == 10.): 
		# 	self.d_can = 10.
		self.d_can = 12.*10**(-3) #meters
		self.length_can = 33.*10**(-3) #meters
		#self.BSA = BSA #patient specific 


	def linearize_func(self,t,v,func,eps_dv=1e-4): 

		# Input  : func(t,v), typically t = t_n+1 ie t_n + dt
		# Return : d/dv (func) at (t,v) using central differences

		dv = eps_dv*v 
		return (func(t,v+dv)-func(t,v-dv))/(2.*dv)

	def linearize_func_t(self,t,v,func,dt=1e-3): 

		return (func(t+dt,v)-func(t-dt,v))/(2.*dt)

	def check_block_consistency(self): 
		if len(connecting_block_list) != 2:
			raise Exception("Chamber block can be connected only to two elements")
		return

	def local_eq_coe_def(self,args): 
		t = args['Time']
		curr_y = args['Solution']
		# dt = args['Time step']
		# rho = args['rho']
		# alpha_f = 1/(1+rho)

		wire_dict = args['Wire dictionary']
		# Previous time step flow in and out
		qi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1]]
		qo = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[1]]
		pi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[0]]
		po = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[0]]

		v = curr_y[self.LPN_solution_ids[0]] 


		# print self.name,' has volume: ',v
		# if v<=0: 
		# 	print "Time: ",t
		# 	raise Exception("Zero chamber volume detected in chamber: ",self.name)

		# Misc constants
		rho = 1060. #density kg/m^3
		mu = 0.004 #viscosity kg/m*s 
		R_cannula = (self.length_can*8*np.pi*mu)/((np.pi*(self.d_can/2)**2))**2 #Pouseille

		denom_out = 1
		denom_in = 1

		Psv_coe = -1. 
		Pao_coe = 1.
		Pcompcoe_in = 1.
		Pcompcoe_out = -1.
							       #V 		       #Qin				                 #Qout 
		self.emxcoe = [(0.,0.,0.,0.,1.),     (0.,denom_in,0.,0.,0.),     (0.,0.,0.,denom_out,0.)]
		self.fmxcoe = [(0.,-1.,0.,1.,0.), (Psv_coe,0.,0.,0.,0.,), (0.,0.,Pao_coe,0.,0.)]
		self.cveccoe = [0.,self.Pcomp(t),(-self.Pcomp(t))]

		self.demxcoe = [(0,)*5,(0,)*5,(0,)*5]
		self.dfmxcoe = [(0,)*5,(0,0,0,0,0),(0,0,0,0,0)]
		self.dcmxcoe = [(0,)*5,(0,)*5,(0,)*5]

		# Checks on volume 
		# if (v >= self.volume): 
		# 	self.emxcoe[1]=(0.,)*5
		# 	self.fmxcoe[1] = (0.,1.,0.,0.,0.)
		# 	self.cveccoe[1]=0
		# 	self.dfmxcoe[1]=(0,)*5
		# elif(v <= 0): 
		# 	self.emxcoe[2]=(0.,)*5
		# 	self.fmxcoe[2] = (0.,0.,0.,1.,0.)
		# 	self.cveccoe[2] = 0
		# 	self.dfmxcoe[2]=(0,)*5


		# Confusing stuff including, but not limited to, Rsuc(ks) *badum tssss* 

	# def friction_coe(self, q): 
	# 	mu = 0.0035 #Pa*s
	# 	Re_numerator = q * (self.d_can/10.)
	# 	Re_denominator = mu * (np.pi/4) * (self.d_can/10.)**2
	# 	Re = abs(Re_numerator/Re_denominator)
	# 	if (Re > 0 and Re <= 2300): 
	# 		return 64./Re
	# 	else: 
	# 		return 10

	# def calc_R_suc(self, P_can, P_sv): #P_sv = P_in 
	# 	V_sv0 = 6.18
	# 	if P_sv > V_sv0: 
	# 		return 0 
	# 	elif(P_sv >= P_th):
	# 		return 0.2623 * (0.9787**P_can - 1) * self.BSA**(-0.3492) 

		#mu = 0.0035 #Pa*s
		# kL_in = 0.8 #some cannula variable
		# kL_out = 0.45 #some cannula variable but out 

		# Friction coefficients 
		# f_in = self.friction_coe(qi)
		# f_out = self.friction_coe(qo)
 
		# Cannula resistance 
		# R_cannula_in = ((8.*f_in*rho)/(np.pi**2))*((0.028/0.009**5)+(0.242/(0.012**5)))*7.5*10**(-15)
		# R_cannula_out = ((8.*f_out*0.28*rho)/(np.pi**2*0.012**5))*7.5*10**(-15)

		# Cannula inductance 
		#L_in = rho*((0.028/((np.pi/4.)*0.009**2))+(0.242/((np.pi/4.)*0.012**2)))*7.5*10**(-9)
		#L_out = (rho*0.028)/((np.pi/4.)*0.012**2)*7.5*10**(-9)
		#L_cannula = (rho*self.length_can)/(np.pi*(self.d_can/2)**2) #Pouseille 

		# Denominator terms for coefficient calculations 
		#denom_in = (np.pi/4.) * self.d_can * L_cannula
		#denom_out = (np.pi/4.) * self.d_can * L_cannula

		# R_suc = self.calc_R_suc(P_can, pi):
		# b1 = (0.00375*kL_in*rho*10**(-6))
		# b2 = (0.00375*kL_out*rho*10**(-6))
		# Psv_coe = -1.*(self.volume**(-0.667))
		# Pao_coe = 1.*(self.volume**(-0.667))

		# a1 = R_cannula
		# a2 = R_cannula

		# c1_lin = a1*abs(qi)*qi
		# c2_lin = a2*abs(qo)*qo
		# c1_lin = 0 
		# c2_lin = 0 

		# self.fmxcoe[0] = v eqn; [1] = dQin; [2] = dQout
		# using previous qi, qo as coefficients for the qi^2, qo^2 terms 
		# self.fmxcoe = [(0.,-1.,0.,1.,0.), (Psv_coe,a1*qi+b1,0.,0.,0.),  (0.,0.,Pao_coe,a2*qo+b2,0.)]
		# self.fmxcoe = [(0.,-1.,0.,1.,0.), (Psv_coe,2*a1*abs(qi)+b1,0.,0.,0.),  (0.,0.,Pao_coe,2*a2*abs(qo)+b2,0.)]
		# self.fmxcoe = [(0.,-1.,0.,1.,0.), (Psv_coe,0.,0.,0.,0.,), (0.,0.,Pao_coe,2*a2*abs(qo)+b2,0.)]
		# self.cveccoe = [0.,Pcompcoe_in*self.Pcomp(t)-c1_lin,Pcompcoe_out*self.Pcomp(t)-c2_lin]
		# print "Pcomp at ", t, ": ", self.Pcomp(t)


		# self.fmxcoe = [(0,-1.,0,1.,0), (-1.,0,1.,0,0),  (1.,0,0.,0,-(a_lin+p_lin))]
		# self.cveccoe = [0,0,c_from_lin+self.Pref]


		#Checks on flow: 
		# if (qi < 0): 
		# 	self.fmxcoe[1] = (Psv_coe,0.,0.,0.,0.)
		# if (qo < 0):
		# 	self.fmxcoe[2] = (0.,0.,Pao_coe,0.,0.)

		# Add foing term to eqn, 0.25mL ~ empty 
		# if volume < 0.01 * v_berlin
		# lambda = 100 * (volume_berlin - 0.01)   
		# else lambda = 0 




# -- VAD -- assuming V_Berlin = 25 mL; will build in input/if's later 
class PulsatileVAD2(LPNBlock):
	def __init__(self, volume, Pcomp, connecting_block_list=None,name="NoNameVAD",flow_directions=None):

		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)

		self.neq = 3
		self.type = "PulsatileVAD"
		self.num_block_vars = 1 # Chamber volume
		self.volume = volume 
		self.Pcomp = Pcomp
		# if (self.volume == 25. or self.volume == 50.): 
		# 	self.d_can = 12.
		# elif (self.volume == 10.): 
		# 	self.d_can = 10.
		self.d_can = 12.*10**(-3) #meters
		self.length_can = 33.*10**(-3) #meters
		#self.BSA = BSA #patient specific 


	def linearize_func(self,t,v,func,eps_dv=1e-4): 

		# Input  : func(t,v), typically t = t_n+1 ie t_n + dt
		# Return : d/dv (func) at (t,v) using central differences

		dv = eps_dv*v 
		return (func(t,v+dv)-func(t,v-dv))/(2.*dv)

	def linearize_func_t(self,t,v,func,dt=1e-3): 

		return (func(t+dt,v)-func(t-dt,v))/(2.*dt)

	def check_block_consistency(self): 
		if len(connecting_block_list) != 2:
			raise Exception("Chamber block can be connected only to two elements")
		return

	def local_eq_coe_def(self,args): 
		t = args['Time']
		curr_y = args['Solution']

		wire_dict = args['Wire dictionary']
		# Previous time step flow in and out
		qi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1]]
		qo = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[1]]
		pi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[0]]
		po = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[0]]

		v = curr_y[self.LPN_solution_ids[0]] 

		# Misc constants
		rho = 1060. #density kg/m^3
		mu = 0.004 #viscosity kg/m*s 
		#mu = 0.0035 #Pa*s

		km = 1.0
		area_berlin = 50*10**(-4)
					 #pi  qi po qo v 	      		       			                  
		self.emxcoe = [(0.,0.,0.,0.,1.),     (0.,0,0.,0.,0.),     (0.,0,0.,0,mu)]
		# self.fmxcoe = [(0.,-1.,0.,1.,0.), (1,0.,-1,0.,0.,), (area_berlin**2,0.,0,0.,km)]	

		# Normalized coefficients
		self.fmxcoe = [(0.,-1./self.volume,0.,1./self.volume, 0.), (1.,0.,-1.,0.,0.,), (area_berlin**2,0.,0.,0.,km/self.volume)]


		self.cveccoe = [0.,0,(-self.Pcomp(t))*area_berlin**2]


		# self.fmxcoe = [(0,-1.,0,1.,0), (-1.,0,1.,0,0),  (1.,0,0.,0,-(a_lin+p_lin))]
		# self.cveccoe = [0,0,c_from_lin+self.Pref]

		self.demxcoe = [(0,)*5,(0,)*5,(0,)*5]
		self.dfmxcoe = [(0,)*5,(0,0,0,0,0),(0,0,0,0,0)]
		self.dcmxcoe = [(0,)*5,(0,)*5,(0,)*5]

		# # Checks on volume 
		# if (v >= self.volume): 
		# 	self.emxcoe[1]=(0.,)*5
		# 	self.fmxcoe[1] = (0.,1.,0.,0.,0.)
		# 	self.cveccoe[1]=0
		# 	self.dfmxcoe[1]=(0,)*5
		# elif(v <= 0): 
		# 	self.emxcoe[2]=(0.,)*5
		# 	self.fmxcoe[2] = (0.,0.,0.,1.,0.)
		# 	self.cveccoe[2] = 0
		# 	self.dfmxcoe[2]=(0,)*5

		#Checks on flow: 
		# if (qi < 0): 
		# 	self.fmxcoe[1] = (Psv_coe,0.,0.,0.,0.)
		# if (qo < 0):
		# 	self.fmxcoe[2] = (0.,0.,Pao_coe,0.,0.)

		# Add forcing term to eqn, 0.25mL ~ empty 
		# if volume < 0.01 * v_berlin
		# lambda = 100 * (volume_berlin - 0.01)   
		# else lambda = 0 



class ContinuousVAD(LPNBlock):
	def __init__(self,RPM,coeff=None, connecting_block_list=None,name="NoNameContinuousVAD",flow_directions=None):
		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
		self.type = "TimeDependentCapacitance"
		rangeRPM = [1800 + 200*r for r in range(0, 9)]
		self.RPM = RPM
		if self.RPM >= rangeRPM[0] and self.RPM < rangeRPM[1]: 
			self.coeff = [-0.0037, 0.0186, 43.674]
		elif self.RPM >= rangeRPM[1] and self.RPM < rangeRPM[2]:
			self.coeff = [-0.0042, 0.0558, 53.676]
		elif self.RPM >= rangeRPM[2] and self.RPM < rangeRPM[3]: 
			self.coeff = [-0.0046, 0.0883, 65.357]
		elif self.RPM >= rangeRPM[3] and self.RPM < rangeRPM[4]: 
			self.coeff = [-0.0049, 0.0882, 78.768]
		elif self.RPM >= rangeRPM[4] and self.RPM < rangeRPM[5]: 
			self.coeff = [-0.0048, 0.01017, 91.875]
		elif self.RPM >= rangeRPM[5] and self.RPM < rangeRPM[6]: 
			self.coeff = [-0.0040, 0.0259, 108.02]
		elif self.RPM >= rangeRPM[6] and self.RPM < rangeRPM[7]:
			self.coeff = [-0.0042, 0.0675, 122.83]
		elif self.RPM >= rangeRPM[7] and self.RPM < rangeRPM[8]:
			self.coeff = [-0.0046, 0.1137, 139.02]
		elif self.RPM >= rangeRPM[8]: 
			self.coeff = [-0.0042, 0.0600, 158.29]
		else: 
			raise Exception("Invalid pressure range for continuous VAD")


	def check_block_consistency(self): 
		if len(connecting_block_list) != 2:
			raise Exception("Capacitance block can be connected only to two elements")


	def local_eq_coe_def(self,args): 
		t = args['Time']
		curr_y = args['Solution']
		wire_dict = args['Wire dictionary']
		qi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1]]
		a = self.coeff[0]
		b = self.coeff[1]
		c = self.coeff[2]
		# a = -0.0042
		# b = 0.06
		# c = 158.29

		c_lin = (a*qi*qi+b*qi+c)-(2*a*qi*qi+b*qi)
		# dt = args['Time step']
		# rho = args['rho']
		# alpha_f = 1/(1.0+rho)	
		self.emxcoe = [(0,0,0,0),(0,0,0,0)]
					   #Pin 	#Q    #Pout
		self.fmxcoe = [(-1.,-(2*a*qi+b),1.,0),(0,1.,0,-1.)]
		self.cveccoe = [-c_lin,0]

		self.demxcoe = [(0,)*4]*2
		self.dfmxcoe = [(0,0,0,0),(0,)*4]
		self.dcmxcoe = [(0,)*4]*2

# -- Inductance
class Inductance(LPNBlock):
	def __init__(self,L,connecting_block_list=None,name="NoNameInductance",flow_directions=None):
		
		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
		self.type = "Inductance"
		self.L = L

	def check_block_consistency(self): 
		if len(connecting_block_list) != 2:
			raise Exception("Inductance block can be connected only to two elements")

	def local_eq_coe_def(self,args): 
		self.emxcoe = [(0.0,-1.,0.0,0),(0,0,0,0)]
		self.fmxcoe = [(1./self.L,0.0,-1./self.L,0),(0,1.,0,-1.)]
		self.cveccoe = [0,0]

		self.demxcoe = [(0,)*4]*2
		self.dfmxcoe = [(0,)*4]*2
		self.dcmxcoe = [(0,)*4]*2		


# -- Ideal diode - state variable
class IdealDiode2(LPNBlock):
	def __init__(self,connecting_block_list=None,eps=1e-17,name="NoNameIdealDiode",flow_directions=None):
		
		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
		self.type = "IdealDiode"
		self.neq=3
		self.num_block_vars = 1 # State : 0- shunt, 1- interrupt
		self.eps =eps
		# self.Peps = Peps

	def check_block_consistency(self): 
		if len(connecting_block_list) != 2:
			raise Exception("IdealDiode block can be connected only to two elements")

	def local_eq_coe_def(self,args): 
		# Needs to take current solution vector as argument
		t = args['Time']
		curr_y = args['Solution']
		wire_dict = args['Wire dictionary']
		rho = args['rho']
		dt = args['Time step']
		alpha_f = 1.0/(1.0+rho)
		Qi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1]]		
		Qo = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[1]]				
		Pi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[0]]
		Po = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[0]]		
		state = curr_y[self.LPN_solution_ids[0]]
		self.emxcoe = [(0,)*4]*3

		if state > 0 : 
			# Zero flow
			self.fmxcoe = [(0,1,0,0,0),(0,0,0,1.,0)]
			self.cveccoe = [0.,0.]

		else : 
			# Perfect wire
			self.fmxcoe = [(1.,0,-1.,0,0),(0,1.,0,-1.,0)]
			self.cveccoe = [0,0]

		if Qi < -self.eps or Qo < -self.eps : 
			self.fmxcoe.append((0,0,0,0,1))
			self.cveccoe.append(-alpha_f-(1-alpha_f)*state)
			# self.cveccoe.append(-1)
			# print 'Flow toggle'
			# print state

		elif Pi - Po > -self.eps : 
			self.fmxcoe.append((0,0,0,0,1))
			self.cveccoe.append(-(1-alpha_f)*state)
			# self.cveccoe.append(0)
			# print 'P toggle'
			# print state
		else : 
			self.fmxcoe.append((0,0,0,0,1))
			self.cveccoe.append(-state)			
			# print 'State maintain'
			# print state



# -- Pressure source (two terminal : analogous to a battery)
class PressureSource(LPNBlock):
	def __init__(self,Pfunction,Pref,connecting_block_list=None,name="NoNamePressureSource",flow_directions=None):
		
		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
		self.neq = 2 
		self.type = "PressureSource"
		self.Pfunction=Pfunction
		self.Pref= Pref

	def check_block_consistency(self): 
		if len(connecting_block_list) != 2:
			raise Exception("PressureSource block can be connected only to two elements")

	def local_eq_coe_def(self,args):
		t = args['Time']
		# rho = args['rho']
		# dt = args['Time step']
		# alpha_f = 1.0/(1.0+rho)

		self.emxcoe = [(0,0)]*2
		self.fmxcoe = [(1.,0),(0,0,1.,0)]
		self.cveccoe = [-1.0*self.Pref,-1.0*self.Pfunction(t)]
		# self.cveccoe = [-1.0*self.Pref,-1.0*self.Pfunction(t)]

		self.demxcoe = [(0,)*4]*2
		self.dfmxcoe = [(0,)*4]*2
		self.dcmxcoe = [(0,)*4]*2


# # -- Flow source (single terminal)
# class FlowSource(LPNBlock):
# 	def __init__(self,Qfunction,connecting_block_list=None,name="NoNameFlowSource",flow_directions=None):
		
# 		LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
# 		self.neq = 1 
# 		self.type = "FlowSource"
# 		self.Qfunction=Qfunction

# 	def check_block_consistency(self): 
# 		if len(connecting_block_list) != 1:
# 			raise Exception("FlowSource block can be connected only to one elements")

# 	def local_eq_coe_def(self,args):
# 		t = args['Time'] 

# 		# rho = args['rho']
# 		# dt = args['Time step']		
# 		# alpha_f = 1.0/(1.0+rho)

# 		self.emxcoe = [(0,0)]
# 		self.fmxcoe = [(0,1.)]
# 		# self.cveccoe = [-1.0*self.Qfunction(t+alpha_f*dt)]
# 		self.cveccoe = [-1.0*self.Qfunction(t)]

# 		self.demxcoe = [(0,)*4]*2
# 		self.dfmxcoe = [(0,)*4]*2
# 		self.dcmxcoe = [(0,)*4]*2

#========================================================================
# 			Connection functions ---~---> 
#========================================================================

def check_block_pair_flow_consistency(bA,bB): 
	if bB.name not in bA.connecting_block_list: 
		raise Exception('Block '+bB.name+' not in connecting list of '+bA.name)
	else: 
		id_bB = bA.connecting_block_list.index(bB.name)

	if bA.name not in bB.connecting_block_list: 
		raise Exception('Block '+bA.name+' not in connecting list of '+bB.name)	
	else: 
		id_bA = bB.connecting_block_list.index(bA.name)

	if bA.flow_directions[id_bB]*bB.flow_directions[id_bA] != -1 : 
		print 'Flow direction of '+bB.name+' :',bB.flow_directions[id_bA]
		print 'Flow direction of '+bA.name+' :',bA.flow_directions[id_bB]
		raise Exception('Flow directions of '+bA.name+' do not conform to that of '+bB.name)

def connect_blocks_by_inblock_list(block_list): 

	connectivity = []

	wire_dict = {}

	bnames = [ _.name for _ in block_list]

	# Check if connection definition is consistent
	for bA in block_list: 
		for bBnm in bA.connecting_block_list:
			bB = block_list[bnames.index(bBnm)]
			check_block_pair_flow_consistency(bA,bB)  

	# If you reached here, it means each block has a consistent (connecting_block_list) and (flow_directions)
	for bA in block_list :
		i = -1  
		id_bA = block_list.index(bA)
		for bBnm in bA.connecting_block_list:
			id_bB = bnames.index(bBnm)
			bB = block_list[id_bB]
			i+=1 
			if bA.flow_directions[i]==+1 and (id_bA,id_bB) not in connectivity : 
				name_wire = bA.name+'_'+bB.name
				connecting_elements = (block_list[id_bA],block_list[id_bB])
				# wire_dict[name_wire] = wire(connecting_elements,name=name_wire)
				connectivity.append((id_bA,id_bB))
			elif bA.flow_directions[i]==-1 : 
				name_wire = bB.name+'_'+bA.name
				connecting_elements = (block_list[id_bB],block_list[id_bA])
			# 	block_list[id_bA].add_connecting_wire(name_wire)
			# 	block_list[id_bB].add_connecting_wire(name_wire)
			else :
				continue 
			wire_dict[name_wire] = wire(connecting_elements,name=name_wire)
			block_list[id_bA].add_connecting_wire(name_wire)

	return connectivity,wire_dict

def connect_blocks_by_connectivity_list(block_list,connectivity): 

	wire_dict = {}

	for e in connectivity: 
		e1,e2 = e
		e1name = block_list[e1].name
		e2name = block_list[e2].name

		connecting_elements = (block_list[e1],block_list[e2])
		name_wire = e1name+'_'+e2name

		wire_dict[name_wire] = wire(connecting_elements,name=name_wire)

		if e2name not in block_list[e1].connecting_block_list: 
			block_list[e1].add_connecting_wire(name_wire)
			block_list[e1].add_connecting_block(e2name,+1)

		if e1name not in block_list[e2].connecting_block_list: 
			block_list[e2].add_connecting_wire(name_wire)
			block_list[e2].add_connecting_block(e1name,-1)

		# print name_wire
		# print block_list[e1].name, block_list[e1].flow_directions
		# print block_list[e2].name, block_list[e2].flow_directions

	#print wire_dict

	return wire_dict	


def check_block_connection(block): 


	if len(block.flow_directions) != block.num_connections : 

		print "Block name: "+block.name
		print "Block number of flows: ",len(block.flow_directions)
		print "Block number of eqs: ",block.num_connections

		raise Exception("Number of connections donot match the number of inflows+outflows for this block")


	# print block.connecting_wires_list
	reorder_inblock_connectivity(block)


# Reorder blocks to have connecting_block_list and connecting_wires_list arranged in ascending flow_directions 
# This will give robustness to initial ordering during setup

def reorder_inblock_connectivity(block):

	indx = sorted(range(len(block.flow_directions)),key=lambda k: block.flow_directions[k])

	block.flow_directions = [ block.flow_directions[_] for _ in indx ]
	block.connecting_wires_list = [ block.connecting_wires_list[_] for _ in indx ]
	block.connecting_block_list = [ block.connecting_block_list[_] for _ in indx ]	

# Function to compute number of equations from blocks and wires
def compute_neq(block_list,wire_dict): 

	neq = 0 
	block_vars = 0
	for b in block_list:
		neq += b.neq
		block_vars += b.num_block_vars

	if 2*len(wire_dict.values()) + block_vars != neq : 
		print "Expected number of variables : ", 2*len(wire_dict) + block_vars 
		print "Instead you have: ", neq
		raise Exception('Mismatch between number of variables and equations')

	return neq


def initialize_solution_structures(neq): 
	# Return y,ydot
	return np.zeros(neq),np.zeros(neq)

def initialize_solution_matrices(neq): 
	# Return E,F,C,dE,dF,dC
	return np.zeros((neq,neq)),np.zeros((neq,neq)),np.zeros(neq),np.zeros((neq,neq)),np.zeros((neq,neq)),np.zeros((neq,neq))

def assign_global_ids(block_list,wire_dict): 

	# Ordering of solution variables : 
	# P0,Q0,P1,Q1,...,Pn,Qn, V1,V2,..,Vm

	i = 0

	var_name_list = []

	for w in wire_dict.values(): 
		w.LPN_solution_ids = [i,i+1]
		var_name_list.append('P_'+w.name)
		var_name_list.append('Q_'+w.name)		
		i+=2 

	for b in block_list : 
		for j in range(b.num_block_vars): 
			b.LPN_solution_ids.append(i)
			var_name_list.append('var_'+str(j)+'_'+b.name)		
			i+=1

	for b in block_list: 
		for local_id in range(b.num_block_vars+2*len(b.connecting_block_list)): 
			b.eq_id_list.append(b.eqids(wire_dict,local_id))


	# print var_name_list

	return var_name_list

def assemble_structures(E,F,C,dE,dF,dC,args,block_list): 
	ieq = 0 
	wire_dict = args['Wire dictionary']

	for b in block_list : 
		# print "Setting up E equations for "+b.name
		b.local_eq_coe_def(args)
		for tpl in b.emxcoe :
			ientry = 0 
			for t in tpl : 
				E[ieq,b.eq_id_list[ientry]] = t
				ientry += 1
			# if len(tpl)!= 0 : 
				# print E[ieq,:] 
			ieq+=1


	ieq = 0
	for b in block_list : 
		# print "Setting up F equations for "+b.name		
		for tpl in b.fmxcoe :
			ientry = 0 
			for t in tpl : 
				F[ieq,b.eq_id_list[ientry]] = t
				ientry += 1
			# if len(tpl)!= 0: 
			# 	print F[ieq,:] 
			ieq+=1	

	ieq = 0 
	for b in block_list : 
		# print "Setting up C equations for "+b.name		
		for t in b.cveccoe :
			# print b.cveccoe, t
			C[ieq] = t
			# if len(tpl)!= 0: 			
			# 	print C[ieq] 
			ieq+=1


	ieq = 0
	for b in block_list : 
		# print "Setting up dE equations for "+b.name		
		for tpl in b.demxcoe :
			ientry = 0 
			for t in tpl : 
				dE[ieq,b.eq_id_list[ientry]] = t
				ientry += 1
			# if len(tpl)!= 0: 
			# 	print F[ieq,:] 
			ieq+=1	

	ieq = 0
	for b in block_list : 
		# print "Setting up dF equations for "+b.name		
		for tpl in b.dfmxcoe :
			ientry = 0 
			for t in tpl : 
				dF[ieq,b.eq_id_list[ientry]] = t
				ientry += 1
			# if len(tpl)!= 0: 
			# 	print F[ieq,:] 
			ieq+=1	

	ieq = 0
	for b in block_list : 
		# print "Setting up dF equations for "+b.name		
		for tpl in b.dcmxcoe :
			ientry = 0 
			for t in tpl : 
				dC[ieq,b.eq_id_list[ientry]] = t
				ientry += 1
			# if len(tpl)!= 0: 
			# 	print F[ieq,:] 
			ieq+=1


# Mid code Taoisms: 
# Implicit codes have 3 use cases for PDE - 1) Over resolved meshes; 2) Skipping initial transients for steady soln ; 
# 3) Handling instabilities due to quasi-steadiness. DAEs form the 4th use-case. 

#=================================================================================
# Generalized alpha integration code - with Newton Raphson
#=================================================================================

# Generalized alpha matrix solve for constant coe E,F,C 
# Mid-steps re-added for ya_f and ydota_m

def form_matrix_NR(E,F,dE,dF,dC,alpha_f,alpha_m,gamma,dt): 

	ydot_update = E*alpha_m/(alpha_f*gamma*dt)

	return (F+(dE+dF+dC+ydot_update)) 

def form_rhs_NR(E,F,C,y,ydot):

	v1 = np.dot(E,ydot)
	v2 = np.dot(F,y)
	rhs = -v1-v2-C
	# print "Residual: ",np.linalg.norm(rhs)
	return rhs


def min_ydot_least_sq_init(neq,eps_min,yinit,block_list,args,dt,rho,eps_init=10.0,eps_factor=5.0):
	# System : min (over y) ||Fy+C||^2 + eps||(Ay-yinit)/yinit||^2
	# Inversion equation: 
	# y <-- inv(F'F+eps*D^2) (-F'C+eps*D*yinit)
	# yinit : Desired set of initial conditions
	# D : diag([(i in yinit)/(abs(yinit)) for i in range(neq) ] )
	# If yinit == zeros, set D = I 

	# ydot is solved in an approximate sense: ydot <-- np.linalg.lstsq(E,-Fy-C)

	# Solve as a sequence of problems : eps ---> eps_min 
	eps = eps_init
	iit = 0
	args['Time'] = 0	
	# y0 = np.zeros(neq)

	y0 = yinit

	E,F,C,dE,dF,dC = initialize_solution_matrices(neq)


	if np.linalg.norm(yinit) == 0.: 
		D = np.eye(neq)
		D2 = D
	else : 
		D = np.diag([ 1.0/abs(_) if _!=0 else 0. for _ in yinit ])
		D2 = np.diag([ 1.0/(_*_) if _!=0 else 0. for _ in yinit ])

	print "Approximate consistent initialization : \n\n"

	while eps > eps_min : 
		iit +=1
		args['Solution'] = y0
		assemble_structures(E,F,C,dE,dF,dC,args,block_list)

		M = np.dot(F.transpose(),F)+eps*D2
		v = -np.dot(F.transpose(),C) + eps*np.dot(D,yinit)
		y0,_,_,_ = np.linalg.lstsq(M,v)

		ydot0,_,_,_ = np.linalg.lstsq(E,-np.dot(F,y0)-C)

		print "Iteration ",iit,", Initializing residual: ",np.linalg.norm(form_rhs_NR(E,F,C,y0,ydot0))
		eps = eps/eps_factor
	return y0,ydot0


# Equation: E*ydot + F*y + C = 0 
def gen_alpha_dae_integrator_NR(y,ydot,t,block_list,args,dt,rho,nit=16): 

	# Constants for generalized alpha
	alpha_m = 0.5*(3.0-rho)/(1.0+rho)
	alpha_f = 1.0/(1.0+rho)
	gamma = 0.5 + alpha_m - alpha_f
	n = y.shape[0]

	damping_step = 1.5

	# # Initial guess for n+1-th step -- steady state type guess
	# curr_y = np.copy(y)
	# curr_ydot = np.copy(ydot)*((gamma-1.0)/gamma)


	# Initial guess for n+1-th step -- explicit euler type guess
	# curr_y = y+dt*ydot
	# curr_ydot = np.copy(ydot)

	# Initial guess for n+1-th step -- explicit euler type guess, half step
	curr_y = y+0.5*dt*ydot
	curr_ydot = np.copy(ydot)*((gamma-0.5)/gamma)


	# Substep level quantities
	yaf = y + alpha_f*(curr_y-y)
	ydotam = ydot + alpha_m*(curr_ydot-ydot)

	# print t
	args['Time'] = t+alpha_f*dt	

	iit = 0

	args['Solution'] = yaf 

	E,F,C,dE,dF,dC = initialize_solution_matrices(n)
	assemble_structures(E,F,C,dE,dF,dC,args,block_list)

	res0 = form_rhs_NR(E,F,C,yaf,ydotam)
	res = res0
	while max(abs(res0)) > 1e-9 and iit < nit: 

		damping = 1.

		M = form_matrix_NR(E,F,dE,dF,dC,alpha_f,alpha_m,gamma,dt)
		res = form_rhs_NR(E,F,C,yaf,ydotam)
		dy = np.linalg.solve(M,res)
		while np.linalg.norm(res) >= np.linalg.norm(res0)  and damping > 1e-7 : 
			yaf2 = yaf +  damping*dy 
			ydotam2 = ydotam +  damping*alpha_m*dy/(alpha_f*gamma*dt)
			damping /= damping_step
			res = form_rhs_NR(E,F,C,yaf2,ydotam2)

		yaf = yaf +  damping_step*damping*dy 
		ydotam = ydotam +  damping_step*damping*alpha_m*dy/(alpha_f*gamma*dt)


		res0 = res

		# Check this equation up
		# ydotam = (1-alpha_m/gamma)*ydot + (alpha_m/(gamma*dt*alpha_f))*(yaf-y)

		args['Solution'] = yaf 

		E,F,C,dE,dF,dC = initialize_solution_matrices(n)
		assemble_structures(E,F,C,dE,dF,dC,args,block_list)
		iit+=1


	if iit>= nit : 
		print "Max NR iterations reached at time: ", t, " , max error: ", max(abs(res0))
		# print "Condition number of F ", np.linalg.cond(F)
		# print "Condition number of NR matrix: ",np.linalg.cond(M)
		# print M

	curr_y = y + (yaf-y)/alpha_f
	curr_ydot = ydot + (ydotam- ydot)/alpha_m

	args['Time'] = t+dt

	return curr_y,curr_ydot
