#module mdcall
#using PyPlot
@everywhere using kin_init
@everywhere using MD_calc
#const plt = PyPlot
#plt.pygui(true)
#plt.close()
#reload("MD_calc")

Ni = 2#2#8; #number of atoms per side in original cubic configuration 
Ts = 1.0; #desired (and initial) temperature in LJ units (Ts = kB*T/epsilon)
ns = 0.05; #LJ number density ns = n*sigma^3

rc = 3.0; #cut-off radius, in LJ units
STEPS = 500#1e4; #number of timesteps
STEPS_thermostat = round(Int64,(STEPS/10)) #3*10^3; #number of steps to leave thermostat on until
STEPS_equilib = round(Int64,(STEPS/2)) #4*10^3; #number of steps before starting to average
alpha = 0.1; #thermostat relaxation parameter 
dt = 0.005; #time step, in LJ units

N = Ni^3; #number of atoms 
Vs = N/ns; #volume of domain, in LJ units, Vs = V/sigma^3
Ls = Vs^(1/3); #length of the domain in each direction, in LJ units
#T_set = Ts*epsilon/kB; #set temperature, in Kelvin
P_LRC = 32/9*pi*ns^2*rc^-9. - 16/3*pi*ns^2*rc^-3.; #long-range P correction
E_LRC = 8/9*pi*ns*rc^-9. - 8/3*pi*ns*rc^-3.; #long-range U correction (per particle)
rc2 = rc^2; #square of cut-off radius

r = zeros(3,N); 
v = zeros(3,N);  
F = zeros(3,N); 
U = zeros(1,N); 

Res = zeros(STEPS,2); 
#positions = zeros(N,3,STEPS);

r,v=kin_init.initialize(Ls,Ni,Ts);
c = [sum(v[1,:]) sum(v[2,:]) sum(v[3,:])]/N;
v = broadcast(-,v,c');
T_inst =2/(3.0*N)*0.5*sum(v.^2)
t=1
println(T_inst)   
#addprocs(4)
F = convert(SharedArray,zeros(3,N))
	U = convert(SharedArray,zeros(1,N))
 	r = convert(SharedArray,r)
	v = convert(SharedArray,v)

@time for t=1:STEPS
	
		
    MD_calc.force_calculation_parallel(N,r,Ls,rc2,F,U);
  	v,r,T_inst = MD_calc.take_a_step_parallel(Ts,dt,N,STEPS_thermostat,alpha,t,F,v,r,Ls)
	T_inst =2/(3.0*N)*0.5*sum(v.^2)
	#println(T_inst)
	if(T_inst > 10.0)
			break
	end
	
	#plt.scatter3D(r[1,:],r[2,:],r[3,:])
	#plt.hold(true)
	
end

#end
