module NeuralECM
    using OCV
    const F = 96485.


    # Build OCV
    A_p = [-123466.90208905171,
            -65215.58255076725,
        -77739.69488015345,
        -132750.8136541972,
        -109183.62667589705]

    A_n = [-9.04729951400691e6,
        -8.911003506471587e6,
        -9.04657963355355e6,
        -8.904669509837592e6,
        -9.0363250622869e6,
        -9.05878665345401e6,
        -9.606335422964232e6,
        -8.023042975317075e6,
        -2.3190474522951595e6,
         1.4303914546788693e6]

    U_0_p = 3.4271972387993173      # U_0_p
    U_0_n = -46.09780594535385    # U_0_n


    c_s_max⁻ = 50000.0
    c_s_max⁺ = 43478.26086956522
    c_s_min⁺ = 0.0
    c_s_min⁻ = 0.0

    cathode_ocv = OCV.RKPolynomial(A_p,U_0_p,c_s_max⁺,c_s_min⁺,length(A_p))
    anode_ocv   = OCV.RKPolynomial(A_n,U_0_n,c_s_max⁻,c_s_min⁻,length(A_n))


    # Geometry
    T⁻ = 86.7e-6
    T⁺ = 66.2e-6
    W⁻ = 2*61.5e-2
    W⁺ = 2*61.5e-2
    L⁺ = 5.8e-2
    L⁻ = 5.8e-2
    Vₛ⁻ = W⁻*L⁻*T⁻
    Vₛ⁺ = W⁺*L⁺*T⁺
    # OCV Initialization
    Vᵢ = 4.2
    x⁻₀ = 0.6
    V = 4.2
    V⁺₀ = V + calcocv(anode_ocv,x⁻₀,297.0)
    x⁺₀ = OCV.get_x_from_voltage(cathode_ocv,V⁺₀,297.0)
    c_n_init = x⁻₀*(anode_ocv.c_s_max-anode_ocv.c_s_min)
    c_p_init = x⁺₀*(cathode_ocv.c_s_max-cathode_ocv.c_s_min)

    function easyOCV(discharge_moles)
        θ⁺ = (c_p_init + (discharge_moles/(Vₛ⁺*0.8)))/c_s_max⁺
        θ⁻ = (c_n_init - (discharge_moles/(Vₛ⁻*0.8)))/c_s_max⁻
        U⁺ = cathode_ocv(θ⁺)
        U⁻ = anode_ocv(θ⁻)
        return U⁺ - U⁻
    end
    function easyOCV(discharge_moles::T) where {T<:AbstractArray}
        θ⁺ = (c_p_init .+ (discharge_moles./(Vₛ⁺.*0.8)))./c_s_max⁺
        θ⁻ = (c_n_init .- (discharge_moles./(Vₛ⁻.*0.8)))./c_s_max⁻
        U⁺ = cathode_ocv.(θ⁺)
        U⁻ = anode_ocv.(θ⁻)
        return U⁺ .- U⁻
    end
end # module NeuralECM
