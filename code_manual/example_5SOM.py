import numpy as np
import os
import subprocess

pi = np.pi
clight = 299792458
E0Au = 196.9665687*931.5e6
E0Elec = 0.51099895000e6
E0P = 938.27208816e6
qe = 1.60217662e-19

working_folder = '20200907_ion_591_store_A0.8_SOM5_fine_mesh'
home = os.getcwd()
cwd = working_folder
try:
    os.mkdir(cwd)
except OSError:
    print ("Creation of the directory %s failed" % cwd)
    
else:
    print ("Successfully created the directory %s" % cwd)
tempinput = {}

        
type_of_particle =  0 # 0 means pronton, 1 means electron, 
csv_out = 0
dynamicOn = 0
n_per_show = 10
turn_start = 0

mainRF = 0
main_detune = 0
detune_slow_factor = 1.0 
R = 610.1754 
GMTSQ = 23.45**2#736.58
Gamma = 293.09#19600.0 

nBeam = 1
beam_shift = 94

beta = np.sqrt(1-1/Gamma**2)
T0 = 2*np.pi*R/(clight*beta)
f0 = 1/T0

#----------------------------#
#important inputs
IbDC = 1
n_turns = 30000
n_dynamicOn = 3000
n_bunches = 630 
n_fill = 1000
n_q_ramp = 2000
n_detune_start = 1000
n_detune_ramp = 3000
n_detune_ramp_tot = 3000
n_I_ramp_start = 1000
n_I_ramp_end = 3000
step_store = 1000
Prad = 1e1#9.72e6
t_rad_long=  1e11#0.03555
Ek_damp = 1e11
Ek_damp_always = 1e11  # to keep the beam stable if needed.
Npar = 1440*10
NperBunch = IbDC/n_bunches/f0/qe
N_bins = 333
fill_step = 12
siglong = 11.368 
A = 0.8

nRF = 5
nRF1 = 5
nRFc = 0
nRF2 = 0
nHOM = 0
nCav = np.array([2,2,2,2,2])
h = np.array([7560,7371,7422,7486,7539])
RoQ = np.array([251,9.47e-3/2,4.53e-3/2,3.18e-2/2,1.14e-2/2])*nCav
delay_time = np.array([1e-6,1e-6,1e-6,1e-6,1e-6]) # in unit of second
delay = [(delay_time[i]*f0*h[0])*N_bins for i in range(len(delay_time))] # number of data points corresponding to the delay time.

n_fb_on = np.array([3000,3000,3000,3000,3000])
gII = np.array([0,0,0,0,0])
gQQ = np.array([0,0,0,0,0])
gIQ = np.array([0.0,0,0,0,0])
gQI = np.array([0.0,0,0,0,0])
gIIi = np.array([0.0,0,0,0,0])
gQQi = np.array([0.0,0,0,0,0])
gIQi = np.array([0.0,0,0,0,0])
gQIi = np.array([0.0,0,0,0,0])

epsilon_comb = np.array([1e-2,1e-2,1e-2,1e-2,1e-2])
g_comb = np.array([0,0,0,0,0])

PA_cap = np.array([1.0,1,1,1,1])

#----------------------------#
# the following parameters need to be derived from the input parameters
# the numbers here are just some place holders.
QL = np.array([1,1,1,1,1])
Vref_I = np.array([1,1,1,1,1])
Vref_Q = np.array([1,1,1,1,1])
Iref_I = np.array([1,1,1,1,1])
Iref_Q = np.array([1,1,1,1,1])
I_I_ref_ini = np.array([1,1,1,1,1])
I_I_ref_final = np.array([1,1,1,1,1])
I_Q_ref_ini = np.array([1,1,1,1,1])
I_Q_ref_final = np.array([1,1,1,1,1])
detune = np.array([0.0,0,0,0,0])
detune_ini = np.array([0.0,0,0,0,0])
detune_mid = np.array([0,0,0,0,0])
detune_final = np.array([0,0,0,0,0])

n_ini_CBI = np.array([1])
mu = np.array([0])
CBI_ini_amp = np.array([0])


beta = np.sqrt(1-1/Gamma**2)
T0 = 2*np.pi*R/(clight*beta)
f0 = 1/T0

if int(type_of_particle==2):
    atomicZ = 79
    Ek = Gamma*E0Au
else:
    atomicZ =1
if int(type_of_particle==1):  
    Ek = Gamma*E0Elec
if int(type_of_particle==0):  
    Ek = Gamma*E0P

eta = 1/GMTSQ-1/Gamma**2
if nRF == 1:
    Qs = np.sqrt(h[int(mainRF)]*atomicZ*np.abs(Vref_I[int(mainRF)])*eta/(2*np.pi*Ek))
elif nRF != 1 :
    Qs = np.sqrt(h[int(mainRF)]*atomicZ*np.abs(Vref_I[int(mainRF)]+Vref_I[1])*eta/(2*np.pi*Ek))
bucket_height_need = 1.2e-2
Qs_need = 0.017#h[0]*eta*bucket_height_need/2

print("Qs_need : ",Qs_need)
omegarf = 2*np.pi*(np.array(h)*f0)
omegac = 2*np.pi*(np.array(h)*f0+detune_final)
Trf = 2*np.pi/omegarf
Rsh = [RoQ[i]*QL[i] for i in range(int(nRF))]

Th = 2*np.pi/omegarf[0]
dthat =Th/N_bins
bucket_height = 2*Qs/(h[mainRF]*eta)*Gamma

print(bucket_height)
print(Ek)
print(Qs)

# setup the parameters 

print("Generating input parameters...")
nPar = NperBunch
NC = nCav[0] #+nCav[1]
NF = 0
ND = 0
if nRF == 1:
	NC = nCav[0] #+nCav[1]
	NF = nCav[0]
	ND = 0
elif nRF ==2 or nRF==3 :
	NC = nCav[0]+nCav[1]
	NF = nCav[0]
	ND = nCav[1]

thetaL = np.zeros(nRF)
Vs = np.zeros(nRF)
Vq = np.zeros(nRF)
Phis = np.zeros(nRF)
PhisPhasor = np.zeros(nRF)
PhiIQ = np.zeros(nRF)
PhiIQIg = np.zeros(nRF)
VrefI = np.zeros(nRF)
VrefQ = np.zeros(nRF)
Vreftot = np.zeros(nRF)
IrefI = np.zeros(nRF)
IrefQ = np.zeros(nRF)
QL = np.zeros(nRF)
Rsh = np.zeros(nRF)

Vs[0] = 1
Vq[0] = 28e6
Vs[1] = 1e-10#Vsynch_need
Vq[1] = 1e-10#-Vquard_need/(NF-ND)
Vs[2] = 1e-10#Vsynch_need
Vq[2] = 1e-10#-Vquard_need/(NF-ND)
Vs[3] = 1e-10#Vsynch_need
Vq[3] = 1e-10#-Vquard_need/(NF-ND)
Vs[4] = 1e-10#Vsynch_need
Vq[4] = 1e-10#-Vquard_need/(NF-ND)

print("Vs,Vq: ",Vs,Vq)
print("Qs = ",Qs)
PhisPhasor = np.arctan(Vq/Vs)
Pbeam0 = Prad0/NC # beam power per fundamental cavity, 
Pbeam = Prad0*nPar/nPar0/NC # beam power per fundamental cavity

IbDC = nPar*f0*qe*n_bunches
Pbeam = IbDC*Vs
print("Beam power per cavity: ",Pbeam)


f = h*f0
# convert RoQ from total to per cavity
RoQ = RoQ/nCav
RoQacc = RoQ*2
print("RoQ per cavity: ", RoQ)
print("Number of cavity : ",nCav)

Vreftot = np.sqrt(Vs**2+Vq**2) 
Qbeam0 = Vreftot**2/(RoQacc*np.abs(Pbeam))
Qbeam = Qbeam0
QL[:] = Qbeam[:]
QL[0] = 3e6#Qbeam
QL[1] = 3e6#Qbeam
QL[2] = 3e6#Qbeam
QL[3] = 3e6#Qbeam
QL[4] = 3e6#Qbeam
if nRF == 3:
	PhisPhasor[2] = PhisPhasor[2]+pi
	#QL[2] = 1e8
	#QL[0] = 1e5
	#QL[1] = 1e5
	
Rsh = RoQ*QL

# Now calculate the inputs
thetaL[0] = (ThetaL_min[0]+dThetaL[0]*thetaL_factor)/180.0*pi  # angle between Ig and Vc
if nRF >=2:
	thetaL[1] = (ThetaL_min[1]+dThetaL[1]*thetaL_factor2)/180.0*pi
if nRF ==3:
	thetaL[2] = (ThetaL_min[2]+dThetaL[2]*thetaL_factor3)/180.0*pi 
Vbr = 2*IbDC*Rsh
print("Vbr = ",Vbr)
Vgr = Vreftot/np.cos(thetaL)*(1+Vbr/Vreftot*np.cos(PhisPhasor))

tgPhi = -(Vbr*np.sin(PhisPhasor)/Vreftot+(1+Vbr*np.cos(PhisPhasor)/Vreftot)*np.tan(thetaL))
tgPhi_ini = -np.tan(thetaL)
delta_f_ini = f*(tgPhi_ini/2/QL+np.sqrt((tgPhi_ini/2/QL)**2+1))-f
delta_f = f*(tgPhi/2/QL+np.sqrt((tgPhi/2/QL)**2+1))-f

delta_f_ini[1] = 45452
delta_f[1] = 45452
delta_f_ini[2] = 18195
delta_f[2] = 18195
delta_f_ini[3] = 10731
delta_f[3] = 10731
delta_f_ini[4] = -2157
delta_f[4] = -2157

VrefI = Vreftot*np.sin(PhisPhasor)
VrefQ = -Vreftot*np.cos(PhisPhasor)
VrefI[1] = 0
VrefQ[1] = 0
VrefI[2] = 0
VrefQ[2] = 0
VrefI[3] = 0
VrefQ[3] = 0
VrefI[4] = 0
VrefQ[4] = 0

I_I = Vgr/Rsh*np.sin(PhisPhasor+thetaL) 
I_Q = -Vgr/Rsh*np.cos(PhisPhasor+thetaL)
I_I[1] = 0
I_Q[1] = 0
I_I[2] = 0
I_Q[2] = 0
I_I[3] = 0
I_Q[3] = 0
I_I[4] = 0
I_Q[4] = 0


I_I_ini = Vreftot/(Rsh)/np.cos(thetaL)*np.sin(PhisPhasor+thetaL)
I_Q_ini = -Vreftot/(Rsh)/np.cos(thetaL)*np.cos(PhisPhasor+thetaL)
I_I_ini[1] = 0
I_Q_ini[1] = 0
I_I_ini[2] = 0
I_Q_ini[2] = 0
I_I_ini[3] = 0
I_Q_ini[3] = 0
I_I_ini[4] = 0
I_Q_ini[4] = 0


print("Vnew : ",Vreftot)
print("QL : ",QL)
print("ThetaL : ",thetaL/pi*180, " [degree]")

print("Tan(PhisPhasor) : ",np.tan(PhisPhasor))
print("PhisPhasor : ",PhisPhasor/pi*180)

print("detune tan : ", tgPhi)
print("detune angle : ", np.arctan(tgPhi)/pi*180, " [degree]")
print("delta_f : ",delta_f, " [Hz]")
print("VrefI : ",VrefI)
print("VrefQ : ",VrefQ)

print("II : ",I_I)
print("IQ : ",I_Q)
print("II_ini : ",I_I_ini)
print("IQ_ini : ",I_Q_ini)

print("VrefTot : ",Vreftot)
print("IbDC : ", IbDC)

# now write to the input file
tempinput['type'] = np.array([type_of_particle])
tempinput['csv_out'] = np.array([csv_out])
tempinput['dynamicOn'] = np.array([dynamicOn])
tempinput['n_per_show'] = np.array([n_per_show])
tempinput['turn_start'] = np.array([turn_start])
tempinput['mainRF'] = np.array([mainRF])
tempinput['main_detune'] = np.array([main_detune])
tempinput['detune_slow_factor'] = np.array([detune_slow_factor])
tempinput['R'] = np.array([R])
tempinput['GMTSQ'] = np.array([GMTSQ])
tempinput['Gamma'] = np.array([Gamma])
tempinput['nBeam'] = np.array([nBeam])
tempinput['beam_shift'] = np.array([beam_shift])

tempinput['n_turns'] = np.array([n_turns])
tempinput['n_dynamicOn'] = np.array([n_dynamicOn])
tempinput['n_bunches'] = np.array([n_bunches])
tempinput['n_fill'] = np.array([n_fill])
tempinput['n_q_ramp'] = np.array([n_q_ramp])
tempinput['n_detune_start'] = np.array([n_detune_start])
tempinput['n_detune_ramp'] = np.array([n_detune_ramp])
tempinput['n_detune_ramp_tot'] = np.array([n_detune_ramp_tot])
tempinput['n_I_ramp_start'] = np.array([n_I_ramp_start])
tempinput['n_I_ramp_end'] = np.array([n_I_ramp_end])
tempinput['step_store'] = np.array([step_store])
tempinput['Prad'] = np.array([Prad])
tempinput['t_rad_long'] = np.array([t_rad_long])
tempinput['Ek_damp'] = np.array([Ek_damp])
tempinput['Ek_damp_always'] = np.array([Ek_damp_always])
tempinput['Npar'] = np.array([Npar])
tempinput['NperBunch'] = np.array([nPar])
tempinput['N_bins'] = np.array([N_bins])
tempinput['fill_step'] = np.array([fill_step])
tempinput['siglong'] = np.array([siglong])
tempinput['A'] = np.array([A])

tempinput['nRF'] = np.array([nRF])
tempinput['nRF1'] = np.array([nRF1])
tempinput['nRFc'] = np.array([nRFc])
tempinput['nRF2'] = np.array([nRF2])
tempinput['nHOM'] = np.array([nHOM])
tempinput['nCav'] = np.array(nCav)
tempinput['h'] = np.array(h)
tempinput['I_Q_ref_ini'] = np.array(I_Q_ref_ini)
tempinput['I_Q_ref_final'] = np.array(I_Q_ref_final)
tempinput['RoQ'] = np.array(RoQ)*nCav
tempinput['delay'] = np.array(delay)
tempinput['n_fb_on'] = np.array(n_fb_on)
tempinput['gII'] = np.array(gII)
tempinput['gQQ'] = np.array(gQQ)
tempinput['gIQ'] = np.array(gIQ)
tempinput['gQI'] = np.array(gQI)
tempinput['gIIi'] = np.array(gIIi)
tempinput['gQQi'] = np.array(gQQi)
tempinput['gIQi'] = np.array(gIQi)
tempinput['gQIi'] = np.array(gQIi)
tempinput['PA_cap'] = np.array(PA_cap)
tempinput['epsilon_comb'] = np.array(epsilon_comb)
tempinput['g_comb'] = np.array(g_comb)

tempinput['QL'] = np.array(QL)
tempinput['Vref_I'] = np.array(VrefI)*nCav
tempinput['Vref_Q'] = np.array(VrefQ)*nCav
tempinput['Iref_I'] = np.array(I_I)
tempinput['Iref_Q'] = np.array(I_Q)
tempinput['I_I_ref_ini'] = np.array(I_I_ini)
tempinput['I_I_ref_final'] = np.array(I_I)
tempinput['I_Q_ref_ini'] = np.array(I_Q_ini)
tempinput['I_Q_ref_final'] = np.array(I_Q)
tempinput['detune'] = np.array(delta_f_ini)
tempinput['detune_ini'] = np.array(delta_f_ini)
tempinput['detune_mid'] = np.array(delta_f_ini+(delta_f-delta_f_ini)/2)
tempinput['detune_final'] = np.array(delta_f)

tempinput['n_ini_CBI'] = np.array(n_ini_CBI)
tempinput['mu'] = np.array(mu)
tempinput['CBI_ini_amp'] = np.array(CBI_ini_amp)



fn1 = 'input.txt'
inputfile1 = os.path.join(cwd,fn1)
with open(inputfile1,'w') as wrt_to_input:
	for i in tempinput:
		wrt_to_input.write(str(i)+' ')
		#print(i)
		for j in range(len(tempinput[i])):
			wrt_to_input.write(str(tempinput[i][j])+' ')
			#print(tempinput[i][j])
		wrt_to_input.write('\n')
print("Generated input file.")
# please leave the one that is going be used uncommented, and comment out the other two.

#args = ("../APES")
#args = ("../APESAVX2")
#args = ("../APESGCC")
args = ("~/Dropbox/code/Cpp/APES_pack/APES8.41/run512.sh")
#args = ("~/Dropbox/code/Cpp/APES_pack/APES8.41/runavx2.sh")
#args = ("")
print(cwd)

popen = subprocess.Popen(args, shell = True, stdout=subprocess.PIPE,cwd=cwd)
print("Simulation started...")
err = popen.wait()
output = popen.stdout.read()
print(output.decode("utf-8"))

if "Beam Lost!" in output.decode("utf-8"):
	path = os.path.join(cwd,"Lost_{0:03d}".format(charge_factor)+"{0:03d}".format(thetaL_factor)+"{0:03d}".format(delay_factor)+"{0:03d}".format(gain_factor)+"nmacro{0:.0f}".format(Npar)+"_nBin{0:.0f}".format(N_bins)+"_Idc{0:.2f}A".format(n_bunches*nPar*f0*qe)+"_ThetaL1_{0:.2f}degree".format(180/pi*thetaL[0])+"_3rd_{0:.3f}us".format(Vq3rd_factor)+"_delay{0:.3f}us".format(delay_time[0]/1e-6)+"_gain{0:.2f}_{1:.2f}".format(gII[0],gQQ[0])+"_epsilon_{0:.2e}".format(epsilon_comb[-1]))
	try:
		os.mkdir(path)
	except OSError:
		print ("Creation of the directory %s failed" % path)
	else:
		print ("Successfully created the directory %s" % path)
		print ("===========================================================\n")
	files = os.listdir(cwd)
	result_fn = [i for i in files if i[-3:]=='bin' and i!='par.bin' and i != 'data.bin' and i != 'init.bin']
	for i in result_fn:
		path_result_fn = os.path.join(cwd,i)
		subprocess.call(["mv",path_result_fn,path])
	path_out = os.path.join(cwd,"out")
	subprocess.call(["cp",path_out,path])

	path_in = os.path.join(cwd,"input.txt")
	subprocess.call(["cp",path_in,path])
	# convert the RoQ back to per cavity, otherwise it's wrong
	RoQ = RoQ*nCav
else: 
	path = os.path.join(cwd,"{0:03d}".format(charge_factor)+"{0:03d}".format(thetaL_factor)+"{0:03d}".format(delay_factor)+"{0:03d}".format(gain_factor)+"nmacro{0:.0f}".format(Npar)+"_nBin{0:.0f}".format(N_bins)+"_Idc{0:.2f}A".format(n_bunches*nPar*f0*qe)+"_ThetaL1_{0:.2f}degree".format(180/pi*thetaL[0])+"_3rd_{0:.3f}".format(Vq3rd_factor)+"_delay{0:.3f}us".format(delay_time[0]/1e-6)+"_gain{0:.2f}_{1:.2f}".format(gII[0],gQQ[0])+"_epsilon_{0:.2e}".format(epsilon_comb[-1]))
	try:
		os.mkdir(path)
	except OSError:
		print ("Creation of the directory %s failed" % path)
	else:
		print ("Successfully created the directory %s" % path)
		print ("===========================================================\n")
	files = os.listdir(cwd)
	result_fn = [i for i in files if i[-3:]=='bin' and i!='par.bin' and i!='init.bin']
	for i in result_fn:
		path_result_fn = os.path.join(cwd,i)
		subprocess.call(["mv",path_result_fn,path])
	path_out = os.path.join(cwd,"out")
	subprocess.call(["mv",path_out,path])

	path_in = os.path.join(cwd,"input.txt")
	subprocess.call(["mv",path_in,path])
	# convert the RoQ back to per cavity, otherwise it's wrong
	RoQ = RoQ*nCav
os.chdir(home)
