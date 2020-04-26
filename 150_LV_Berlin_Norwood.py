#-- Berlin Heart using thin plate model for membrane ------#
#--- with native left ventricle contribution in a ----------#
#--- a simplified open loop of systemic circulation. -------#
#--- Norwood procedure w/ values from Boston Children's ----#

import sys
import numpy as np
from network_util_NR import * 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#matplotlib.use('Agg')
import os 
import time

def createPressureLambdaFunction(time, pressure):
    return lambda t: np.interp(t, time, pressure, period=0.5)

Lap = 4. #4.15 #mmHg 
# 60-80 mL/kg 
# Normalization constants
T_CONV = 1.0
P_CONV = 1333.32 # mmHg to dyn/cm^2 
Q_CONV = 1.0 # ml/sec to cm3/sec
R_CONV = P_CONV/Q_CONV 
C_CONV = Q_CONV*T_CONV/P_CONV # Compliance in cgs 
L_CONV = P_CONV*T_CONV/Q_CONV # Inertance from Wood to cgs 
V_CONV = Q_CONV*T_CONV

H_CONV = 1.0
E_CONV = P_CONV


# Time variables 
HR = 150.

# DFR = 0.6 #60./100.
# t_sys = 0.2*t_VAD+(1-DFR)*t_VAD # Function to track systolic time 

tinit = 0
nts= 250000
dt = 0.0001*T_CONV #0.0001*T_CONV
tlist = np.array([tinit+ _*dt for _ in range(1,nts+1)])
COI = [18, 19, 20, 21, 22, 23, 24]

########--------------- LEFT VENT THINGS ----------------########

# Mynard valves for heart; Areas from Boston Children's Zscore 
A_mitral = 1.12 # Boston Children's Z-score
A_aortic = 0.65 # Boston Children's Z-score
rxn = 7. # 5. 2. 
v_mitral = MynardValve(kv=rxn, A_ann=A_mitral, name="vMitral")
v_aortic = MynardValve(kv=rxn, A_ann=A_aortic, name="vAortic")

# Psource = UnsteadyPressureRef(Pfunc= lambda t: (Lap*P_CONV)+P_CONV*100.0*np.max(np.cos(2*np.pi*t/t_VAD),0), name="Lap")
Psource = PressureRef(Pref=(Lap*P_CONV), name="Lap")

LVESV = 5.98 #4.48 #3.485 #1 STD #13.2 
LVEDV = 14.43 #10.76 #14.285 #1 STD #23.8

# print "\nLVESV = ", LVESV, "\nLVEDV = ", LVEDV 
vu = LVESV

a_factor = 7.15 #6.5, 8.
def active_vent(t,v,Vu=vu*V_CONV): #prev Vu=5.13
	t_cyc = 60./HR
	E = 8.5*P_CONV*a_factor
	E2 = -0.042*P_CONV*a_factor
	# E2 = -.1*P_CONV
	# t_off = 0.02*t_cyc
	# t_sys = 0.3*t_cyc+0.16
	t_off = 0. 
	t_sys = (1./3.)*t_cyc

	tn = t - np.floor(t/t_cyc)*t_cyc

	activationfunc = (0.5*(np.sin(np.pi*(tn-t_off)/t_sys)+abs(np.sin(np.pi*(tn-t_off)/t_sys))))  if (tn <= t_sys+t_off and tn >= t_off) else 0.0 
	Pa = activationfunc*(E+E2*(v-Vu))*(v-Vu) 
	return Pa

p_factor = 9.35 #6.
def passive_vent(t,v,Vu=vu*V_CONV): #prev 11.7 
	P0 = 0.9*P_CONV
	Ke = 0.062*p_factor 
	Pp = P0*Ke*(v-Vu) # P0*(np.exp(Ke*(v-Vu))-1.)  
	return Pp

LV = ChamberModel(Activefunc=active_vent, Passivefunc=passive_vent, Pref=0.0, name="LV")

########-------------- Activation funcs ----------------########
def Pa_LV(t):
	t_cyc_Pa = 60./HR
	t_off_Pa = 0.
	t_sys_Pa = (1./3.)*t_cyc_Pa 
	tn_Pa = t - np.floor(t_cyc_Pa)*t_cyc_Pa
	active = (0.5*(np.sin(np.pi*(tn_Pa-t_off_Pa)/t_sys_Pa)+abs(np.sin(np.pi*(tn_Pa-t_off_Pa)/t_sys_Pa)))) if (tn_Pa <= t_sys_Pa+t_off_Pa and tn_Pa >= t_off_Pa) else 0.0
	E_Pa = 8.5*P_CONV*10.
	E2_Pa = -0.024*P_CONV*10.
	return active*(E_Pa+E2_Pa*7.)*7. 

########-------------- BERLIN DECLARATION ---------------########

A = 120. # [140.-->241, 135.-->231, 130.-->221]
C = Lap*-5. #3.25, 2.

#const = 0.
#phase_shift = const*np.pi #(np.pi)/(-2.) 

# print "\nBerlin pressure params: ", A, C, "\n\n"
#print "\nBerlin phase shift: ", phase_shift/np.pi
# Pa_Berlin = lambda t: P_CONV*((A*np.sin(2*np.pi*t/t_VAD)+abs(A*np.sin(2*np.pi*t/t_VAD))) + C) #T = 1.5 
#Pa_Berlin = lambda t: P_CONV*((A*np.sin(2*np.pi*t/t_VAD + phase_shift)+abs(A*np.sin(2*np.pi*t/t_VAD + phase_shift))) + C) #T = 1.5

# phase_shift = 0.5*2 #[0.625, 0.1875, 0.3125, 0.4375]*2. # [0.5, 0.375. 0.25, 0.125, 0.]
# Pa_Berlin = lambda t: P_CONV*((A*np.sin(2*np.pi*t/t_VAD)+abs(A*np.sin(2*np.pi*t/t_VAD))) + C) #T = 1.5 
# def Pa_Berlin(t, c):
# 	t_VAD_off = 0. 
# 	sys_time = 0.4
# 	t_VAD_sys = sys_time*t_VAD      #(1./3.)*t_VAD
# 	c = 0.5*2 #-(0.75/2.)*2 # (+/-)[0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375]
# 	# c = -(6./16.)/2.
# 	# phase = c*np.pi
# 	# tn_VAD = t - np.floor(t/t_VAD)*t_VAD
# 	tn_VAD = t+c*t_VAD  - np.floor((t+c*t_VAD)/t_VAD)*t_VAD  #tracking current time
# 	active_Berlin = P_CONV*((A*np.sin(np.pi*((tn_VAD-t_VAD_off)/t_VAD_sys)) \
# 	+abs((A*np.sin(np.pi*((tn_VAD-t_VAD_off)/t_VAD_sys)))))+C) \
# 	if ( (tn_VAD <= t_VAD_sys)  ) and (tn_VAD >= 0) else P_CONV*C	
# 	return active_Berlin


# Explicit declaration of Berlin 
mu = 0.5
E = 1.*10**6*E_CONV
h = 0.057*H_CONV #0.075<->30mL 0.041<->10mL
volume = 10.
# BERLIN TUBING 
r_mu = 0.04
# can_radius = 0.9/2. # cm
can_diameter = 0.66
can_radius = can_diameter/2. # 10mL 
can_area = np.pi*can_radius**2 

# in_length = 220*10**-1 # mm to cm for 10mL
# out_length = 240*10**-1 # mm to cm for 10mL

in_length = 250*10**-1
out_length = 250*10**-1
Rpouseille_in = (8*np.pi*r_mu*in_length)/(can_area**2)
Rpouseille_out = (8*np.pi*r_mu*out_length)/(can_area**2)
Rin = Resistance(R=Rpouseille_in, name="Rin")
Rout = Resistance(R=Rpouseille_out, name="Rout")

# MYNARD VALVES - 0 chiuso, 1 aperto
rxn_in = 0.5
rxn_out = rxn_in
f = 1.0 
A_Bvalve = (np.pi*(can_radius)**2)
vIn = MynardValve(kv=rxn_in, A_ann=A_Bvalve, name="vIn")
vOut = MynardValve(kv=rxn_out, A_ann=A_Bvalve, name="vOut")

########-------------- LPN BLOCK THINGS ----------------########
def main():
	all_DR = [60., 70., 75., 80., 90., 100., 110., 120.]
	for DR in all_DR: 
		print "HR/DR: ", HR, DR
		t_VAD = 60./DR
		phase = [-0.5, -0.25, 0.125, 0., 0.125, 0.25, 0.5]
		for c in phase:
			print "\nPHASE: ", c 
			def Pa_Berlin(t, c=c):
				t_VAD_off = 0. 
				sys_time = 0.4
				t_VAD_sys = sys_time*t_VAD      #(1./3.)*t_VAD
				# c = 0.5*2 #-(0.75/2.)*2 # (+/-)[0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375]
				# c = -(6./16.)/2.
				# phase = c*np.pi
				# tn_VAD = t - np.floor(t/t_VAD)*t_VAD
				c = c*2. 
				tn_VAD = t+c*t_VAD  - np.floor((t+c*t_VAD)/t_VAD)*t_VAD  #tracking current time
				active_Berlin = P_CONV*((A*np.sin(np.pi*((tn_VAD-t_VAD_off)/t_VAD_sys)) \
				+abs((A*np.sin(np.pi*((tn_VAD-t_VAD_off)/t_VAD_sys)))))+C) \
				if ( (tn_VAD <= t_VAD_sys)  ) and (tn_VAD >= 0) else P_CONV*C	
				return active_Berlin
				
			Berlin = pVADsphere(volume=volume*V_CONV, Pa=Pa_Berlin, h=h, Young=E, Poisson=mu, threshold=0.5, name="Berlin")

			# Simple heart model blocks
			L_la = 5.*10**-6
			# R_c = 90.*0.25*1.e-2
			# Upper body physiology 
			#C_uba = 0.271*10**-2 # Migliavacca, et al. 
			#L_ub = 0.05		 # Migliavacca, et al. 

			# C_d = smush capacitors from UB, LB circulations 
			V_d = 30.*V_CONV
			#C_d = 3*C_uba*C_CONV
			#Q_0 = V_d*V_CONV/C_d
			# Heart
			Lla = Inductance(L=L_la*L_CONV, name="Lla")

			# Rc = Resistance(R=R_c*R_CONV, name="Rc")
			# SVR = 23.8596 
			# R_uba = 34. # Woods 
			# R_c = 533.

			R_CONV_Woods = 80. # Convert from Woods to cgs 
			#R_uba = 4.7662 #[Woods] prev = 8.08
			#R_lba = R_uba
			#C_uba = 21.96 #2.94088 #2.81276  #2.8764
			# R_c = 135.756 # [cgs]
			
			BSA = 0.33
			SVR = 21.6
			R_uba = (SVR/BSA)*2.
			R_lba = R_uba
			MAP = 67. 
			wt = 6.51
			TBV= 85*6.51
			C_ao = 0.05*((TBV-LVEDV)/MAP)
			C_uba = 0.475*((TBV-LVEDV)/MAP)
			R_ao = 0.05*((1./R_uba + 1./R_lba)**-1)
			
			r_factor = 2. 
			R_ub_dist = R_uba/r_factor
			R_ub_prox = R_uba - R_ub_dist
			R_lb_dist = R_ub_dist
			R_lb_prox = R_ub_prox

			C_d = 3*C_uba*C_CONV
			Q_0 = V_d*V_CONV/C_d

			l_factor = 1.  #0.1 
			L = 22.083*l_factor #L_ub*L_CONV
			# L = 0.05
			megaL = Inductance(L=L, name="megaL")
			RCub = RCRBlock(Rp=R_ub_prox*R_CONV_Woods, C=C_uba*C_CONV*0.5, Rd=R_ub_dist*R_CONV_Woods, name="RCub")
			RClb = RCRBlock(Rp=R_lb_prox*R_CONV_Woods, C=C_uba*C_CONV*0.5, Rd=R_lb_dist*R_CONV, name="RClb")
			# r_factor = 1.
			#R_ao = 0.05*((1./R_uba + 1./R_lba)**-1)
			#C_ao = 0.29952 #0.29608 #0.309904 #0.3028
			RCao = RCBlock(R=R_ao*R_CONV_Woods, C=C_ao*C_CONV*0.5, name="RCao")

			print "Unstressed Volume [mL]: ", vu
			print "Aortic Resistance [cgs]: ", RCao.R
			print "Vascular Resistance [cgs]: ", RClb.Rp, RClb.Rd
			print "Systemic Inductance [cgs]: ", megaL.L
			print "Sum of capacitances [cgs]: ", (RCao.C + RCub.C + RClb.C)*P_CONV
			print "Up/Low Capacitances [cgs]: ", RCub.C*P_CONV

			# Grounds 
			GND0 = PressureRef(Pref=8.*P_CONV, name="GND0")
			GND1 = PressureRef(Pref=8.*P_CONV, name="GND1")
			GND2 = PressureRef(Pref=0, name="GND2")
			GND3 = PressureRef(Pref=0, name="GND3")
			GND4 = PressureRef(Pref=0, name="GND4")
			# Jxns 
			J0 = Junction(name="J0")
			J1 = Junction(name="J1")
			J2 = Junction(name="J2")
			J3 = Junction(name="J3")

			#				0		1		2	  3			4	   5	 6	  7.   8	  9	   10	  
			# block_list = [Psource, v_mitral, LV, v_aortic, RCao, megaL, J0, RCub, GND0, RClb, GND1]
			# connect_list = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (6,9), (7,8), (9,10)]

			# 				0		1	  2		 3		4		5	 6 		7	  8		9	 10	  11.  	12	  13   14	 15	   16	17	
			# block_list = [Psource, v_mitral, J0, LV, v_aortic, Rin, vIn, Berlin, vOut, Rout, J1, RCao, megaL, J2, RCub, GND0, RClb, GND1]
			# connect_list = [(0,1), (1,2), (2,3), (3,4), (4,10), (2,5), (5,6), (6,7), (7,8), (8,9), (9,10), (10,11), (11,12), (12,13), (13,14), (14,15), (13,16), (16,17)]

			# 				0		1	  2		 3		4		5	 6 		7	  8		9	 10	  11.  	12	  13   14	 15	   16	17	
			block_list = [Psource, v_mitral, LV, J0, v_aortic, Rin, vIn, Berlin, vOut, Rout, J1, RCao, megaL, J2, RCub, GND0, RClb, GND1]
			connect_list = [(0,1), (1,2), (2,3), (3,4), (4,10), (3,5), (5,6), (6,7), (7,8), (8,9), (9,10), (10,11), (11,12), (12,13), (13,14), (14,15), (13,16), (16,17)]
			wdict = connect_blocks_by_connectivity_list(block_list,connect_list)

			# print [w for w in wdict.keys()]

			# Compute number of equations
			neq = compute_neq(block_list,wdict)

			# Check block level consistency
			for b in block_list: 
			    check_block_connection(b)

			# Assign an ID number to wire and block level variables
			var_name_list = assign_global_ids(block_list, wdict)

			# Print list of solution variables (debug)
			print "\nVariables: ", var_name_list

			# Initialize solution structures
			y0,ydot0 = initialize_solution_structures(neq) 
			ind_p_down = var_name_list.index('P_J1_RCao')
			ind_p_ubody = var_name_list.index('P_J2_RCub')
			ind_p_aorta = var_name_list.index('P_RCao_megaL')
			ind_pin = var_name_list.index('P_LV_J0')
			ind_pout = var_name_list.index('P_J0_vAortic')
			ind_qin = var_name_list.index('Q_LV_J0')
			ind_qout = var_name_list.index('Q_J0_vAortic')
			ind_vAortic = var_name_list.index('var_0_vAortic')
			ind_vMitral = var_name_list.index('var_0_vMitral')
			ind_Qd = var_name_list.index('Q_megaL_J2')
			ind_LV = var_name_list.index('var_0_LV')
			# Berlin parameters of interest 
			ind_pin_Berlin = var_name_list.index('P_vIn_Berlin')
			ind_pout_Berlin = var_name_list.index('P_Berlin_vOut')
			ind_qin_Berlin = var_name_list.index('Q_vIn_Berlin')
			ind_qout_Berlin = var_name_list.index('Q_Berlin_vOut')
			ind_Berlin = var_name_list.index('var_0_Berlin')
			ind_vIn = var_name_list.index('var_0_vIn')
			ind_vOut = var_name_list.index('var_0_vOut')
			ind_ascAo = var_name_list.index('Q_J1_RCao')

			P_aorta = 60.*P_CONV
			for v in var_name_list: 
				if 'P_' in v and 'J0' in v: 
					index = var_name_list.index(v)
					y0[index] = P_aorta
					print v
			y0[ind_Qd] = Q_0
			y0[ind_vAortic] = 0. 
			y0[ind_vMitral] = 1.
			y0[var_name_list.index('P_J1_RCao')] = P_aorta 
			y0[ind_LV] =  LVEDV*V_CONV 
			y0[ind_Berlin] = 0.7*Berlin.volume
			y0[var_name_list.index('var_0_vIn')] = 0.5
			y0[var_name_list.index('var_0_vOut')] = 0.5
			print "Intial conditions: ", y0
			curr_y = y0
			curr_ydot = ydot0

			# rho is the spectral radius for numerical noise dissipation used by Generalized-alpha. Standard value is set to 0.5. 
			# Set rho to zero if backward Euler is needed. 
			rho = 0.0

			# Args passes arguments needed by blocks to define the equations through a dictionary
			args = {}
			args['Time step'] = dt
			args['rho']=rho
			args['Wire dictionary']=wdict
			args['Displacement']=[]


			# raise Exception('Stop now')
			# curr_y,curr_ydot = min_ydot_least_sq_init(neq,1e0,y0,block_list,args,dt,rho,eps_factor=2.5)

			print curr_y
			ylist = [curr_y]

			# Time-stepping loop----->
			for t in tlist:
				if np.mod(t, 0.5) == 0:
					print t
				args['Time'] = t 
				args['Solution'] = curr_y 
				curr_y,curr_ydot = gen_alpha_dae_integrator_NR(curr_y,curr_ydot,t,block_list,args,dt,rho,nit=90) # Change nit to 75 if necessary 
				ylist.append(curr_y)
				args['Displacement'].append(Berlin.displacement_for_plots(args))

			ylist = np.array(ylist)

			# Save arrays 
			filepath = os.getcwd()
			#print filepath
			# phase_label = str(phase_shift).replace('.', '_')
			# h_label = str(Berlin.h).replace('.', '_')
			ylist_file = 'ylist_' + str(HR) + '_' + str(DR) + '.npy'
			args_file  = 'args_' + str(HR) + '_' + str(DR) + '.npy'
			vars_file = 'vars_' + str(HR) + '_' + str(DR) + '.npy'
			phase_label = str(c*2.).replace('.', '_')
			dir_name  = str(Berlin.volume)[:2] +'_' + str(HR) + '_' + str(DR) + '_' + phase_label
			#big_dir = 'HR_' + str(HR) + '_' + str(DR)
			#dir_name = big_dir + '/' + dir_name
			#big_directory = os.path.join(filepath, big_dir)
			#if not os.path.exists(big_directory):
			#	os.mkdir(big_directory)
			directory = os.path.join(filepath, dir_name)
			print directory
			if not os.path.exists(directory):
				os.mkdir(directory)
				print "Directory made"
			np.save(directory + '/' + ylist_file, ylist, allow_pickle=True)
			np.save(directory + '/' + args_file, args['Displacement'], allow_pickle=True)
			np.save(directory + '/' + vars_file, var_name_list, allow_pickle=True)

			# fig=plt.figure()
			# Berlin_Pa = [Berlin.Pa(t) for t in tlist]
			# plt.plot(tlist[0:19999], np.array(Berlin_Pa)[0:19999]/P_CONV)
			# LV_Pa = [Pa_LV(t) for t in tlist]
			# plt.plot(tlist[0:19999], np.array(LV_Pa)[0:19999]/P_CONV)
			# plt.xlabel('Time (s)') 
			# plt.ylabel('Activation Pressure (mmHg)')
			# plt.title('Activation functions of BH and LV')
			# plt.legend(['Berlin', 'LV'])
			# fig.savefig('ActFuncs.png')
			# plt.close(fig)

			#os.system('sbatch pulsatile.sh')
			
			time.sleep(60)

			print ("Sleep done. Next sim.")

if __name__ == "__main__":
	main()
















