@everywhere module MD_calc

#include("particle_manipulations.jl")
#include("potential.jl")
#################################
#Initialization
function set_MD_params(n::Int64=27,
		dt::Float64=0.005,
		tsteps::Int64=100,
		steps::Int64=1000,
		alpha::Float64=0.1,
		rc::Float64=2.5,
		Ts::Float64=1.0)

	N = n
	DT = dt
	STEPS_THERMOSTAT= tsteps
	TS = Ts
	α = alpha
	RC = rc
	STEPS=steps

	return (DT,N,STEPS_THERMOSTAT,TS,α,RC)

end


#################################
#Velocity Verlet
function take_a_step(TS::Float64,
		DT::Float64,
		N::Int64,
		STEPS_THERMOSTAT::Int64,
		α::Float64,
		t::Int64,
		F::Array{Float64,2},
		V::Array{Float64,2},
		POS::Array{Float64,2},LS::Float64)



VCM = mean(V,2);
V_THERM_SQ =0.e0
#for i = 1:3*size(POS,1)
	  V = V +DT/2.0*F
	  V = broadcast(-,V,VCM)
	  #V[i] = V[i] + DT/2.0*F[i] - VCM[ceil(Int64,i/N)]
	 
for i = 1:3*size(POS,2)	  
	 
	  if( t < STEPS_THERMOSTAT)	
		V_THERM_SQ += (V[i])^2.0
	  end
end

T_inst = 2/(3*N)*0.5*V_THERM_SQ;
c = sqrt( 1+α*(TS/T_inst-1) )
if( t < STEPS_THERMOSTAT)	
V = V*c
end

POS = POS + DT*V
[POS[i] = POS[i] >=LS ? POS[i]-LS: POS[i] for i = 1:length(POS)]
[POS[i] = POS[i] <0 ? POS[i]+LS: POS[i] for i = 1:length(POS)]


return V,POS,T_inst
end

############################
#Parallel take_a_step

#Velocity Verlet
function take_a_step_parallel(TS::Float64,
		DT::Float64,
		N::Int64,
		STEPS_THERMOSTAT::Int64,
		α::Float64,
		t::Int64,
		F::SharedArray{Float64,2},
		V::SharedArray{Float64,2},
		POS::SharedArray{Float64,2},LS::Float64)



VCM = mean(V,2);
V_THERM_SQ =0.e0
#for i = 1:3*size(POS,1)
	  V = V +DT/2.0*F
	  V = broadcast(-,V,VCM)
	  #V[i] = V[i] + DT/2.0*F[i] - VCM[ceil(Int64,i/N)]
	 
for i = 1:3*size(POS,2)	  
	 
	  if( t < STEPS_THERMOSTAT)	
		V_THERM_SQ += (V[i])^2.0
	  end
end

T_inst = 2/(3*N)*0.5*V_THERM_SQ;
c = sqrt( 1+α*(TS/T_inst-1) )
if( t < STEPS_THERMOSTAT)	
V = V*c
end


@sync @parallel for i = 1:length(POS)
		POS[i] = POS[i] + DT*V[i]
		POS[i] = POS[i] >=LS ? POS[i]-LS: POS[i] 
		POS[i] = POS[i] <0 ? POS[i]+LS: POS[i] 
end


return V,POS,T_inst
end


#################################

function force_calculation(N::Int64,
						   pos::Array{Float64,2},
						   Ls::Float64,
						   rc2::Float64,
						   F::Array{Float64,2},
						   U::Array{Float64,2}) 
	 

				   
			for i = 1:3:3*(N-1) 
				for j = i+3:3:3*N
					xij = pos[i]-pos[j]
					yij = pos[i+1]-pos[j+1]
					zij = pos[i+2] - pos[j+2]
					if(xij > Ls/2.0)
							xij-=Ls
					end
					if(xij < -Ls/2.0)
							xij+=Ls
					end
					if(yij > Ls/2.0)
							yij-=Ls
					end
					if(yij < -Ls/2.0)
							yij+=Ls
					end
					if(zij > Ls/2.0)
							zij-=Ls
					end
					if(zij < -Ls/2.0)
							zij+=Ls
					end


					rij2 = xij^2.0 + yij^2.0 + zij^2.0
					rij = sqrt(rij2)
					if(rij2 < rc2^2)
						p_i= ceil(Int64,i/3)
						p_j = ceil(Int64,j/3)
						rij6 = 1./(rij2^3.0)
						rij12 = 1./(rij2^6.0)
						c =4.e0*(rij12-rij6)
						d =24.e0/rij2*(2.0*rij12-rij6)
						U[p_i] += c
						U[p_j] += c

						F[i] += d*xij
						F[i+1] += d*yij
						F[i+2] += d*zij
						F[j] += -d*xij
						F[j+1] += -d*yij
						F[j+2] += -d*zij
					end

				end
		end




end

#### Parallel Force Calculation

function force_calculation_parallel(N::Int64,
pos::SharedArray{Float64,2},
Ls::Float64,
rc2::Float64,
F::SharedArray{Float64,2},
U::SharedArray{Float64,2})



@sync @parallel for i = 1:3:3*(N-1)
@sync @parallel for j = i+3:3:3*N
xij = pos[i]-pos[j]
yij = pos[i+1]-pos[j+1]
zij = pos[i+2] - pos[j+2]
if(xij > Ls/2.0)
xij-=Ls
end
if(xij < -Ls/2.0)
xij+=Ls
end
if(yij > Ls/2.0)
yij-=Ls
end
if(yij < -Ls/2.0)
yij+=Ls
end
if(zij > Ls/2.0)
zij-=Ls
end
if(zij < -Ls/2.0)
zij+=Ls
end


rij2 = xij^2.0 + yij^2.0 + zij^2.0
rij = sqrt(rij2)
if(rij2 < rc2^2)
p_i= ceil(Int64,i/3)
p_j = ceil(Int64,j/3)
rij6 = 1./(rij2^3.0)
rij12 = 1./(rij2^6.0)
c =4.e0*(rij12-rij6)
d =24.e0/rij2*(2.0*rij12-rij6)
U[p_i] += c
U[p_j] += c

F[i] += d*xij
F[i+1] += d*yij
F[i+2] += d*zij
F[j] += -d*xij
F[j+1] += -d*yij
F[j+2] += -d*zij
end

end
end
end

end
