#I've loaded all the packages you need
using OrdinaryDiffEq, OCV, CSV, DataFrames, Turing, ForwardDiff, NeuralECM, PythonPlot, DiffEqFlux, ComponentArrays, Parameters, DataInterpolations


# Change this function to change the current
function current_function(t)
    return 3
end

p = ComponentArray(
    R_0 = 0.02,
    R_1 = 0.2, 
    C_1 = 10000.0
)

# This is the model function
function EasyECM(du, u, p, t)
    #Get parameters and states
    @unpack R_0, R_1, C_1 = p
    discharge_moles, I_R1 = u

    #current from the current function!
    I = current_function(t)

    #discharge_moles is the DOD, but in moles of "cycled" Li.
    discharge_moles = u[1]
    du[1] = I/NeuralECM.F # (faraday's constant)

    #from plett
    dIR_1_dt = (-I_R1 + I)/(R_1*C_1)
    du[2] = dIR_1_dt
end

# This converts the solution into voltage
function get_voltage(solution, p, I)
    V = similar(sol.t)
    @unpack R_0, R_1 = p
    V = NeuralECM.easyOCV.(solution[1,:]) .- (sol[2,:].*R_1 .+ I.*R_0)
    return V
end


# Start at 0 depth of discharge
DOD₀ = 0.0
I_R1 = 0.0

u0 = [DOD₀, I_R1]


prob = ODEProblem(EasyECM, u0, (0.0, 3600.0), p)

sol = solve(prob, Tsit5(), dtmax=1.0)
V = get_voltage(sol, p, 3.0)

figure(1)
clf()
plot(sol.t, V)

