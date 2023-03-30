"""
    rooth_smooth(u, p, t; kwargs...)
Rooth's 3-box model of the interhemispheric Thermohaline Circulation, smooth version.

The model is slightly modified by including a salinity restoration term and replacing
the absolute value of the flow strength `q` with a smooth hyperbolic tangent approximation.

## State vector
`u = [T1, T2, T3, S1, S3]`

## Parameters 
`p = [F1, F3, lambd_T, lambd_S, tau1, tau2, tau3, sig1, sig2, sig3, smooth]`

## Keyword arguments
* `k = 1.5e-6`: hydraulic constant
* `alpha = 1.5e-4`: thermal expansion coefficient
* `beta = 8e-4`: haline expansion coefficient
* `V = 2`: volume ratio between equatorial and polar boxes
"""
function rooth_smooth(u, p, t;
    k = 1.5e-6,
    alpha = 1.5e-4,
    beta = 8e-4,
    V = 2,
    S_tot = 140)

    T1, T2, T3, S1, S3 = u
    F1, F3, lambd_T, lambd_S, tau1, tau2, tau3, sig1, sig2, sig3, smooth = p[1]
    tau = [tau1, tau2, tau3]
    sigma = [sig1, sig2, sig3]

    q = k*(alpha*(T3-T1) - beta*(S3-S1))
    S2 = (S_tot - (S1+S3))/V
    c_pos = 1 + tanh(smooth*q)
    c_neg = 1 - tanh(smooth*q)

    dT1 = q/2*(c_pos*(T2-T1) + c_neg*(T1-T3)) + lambd_T*(tau[1] - T1)
    dT2 = q/(2V)*(c_pos*(T3-T2) + c_neg*(T2-T1)) + lambd_T*(tau[2] - T2)
    dT3 = q/2*(c_pos*(T1-T3) + c_neg*(T3-T2)) + lambd_T*(tau[3] - T3)

    dS1 = q/2*(c_pos*(S2-S1) + c_neg*(S1-S3)) + lambd_S*(sigma[1] - S1) - F1
    dS3 = q/2*(c_pos*(S1-S3) + c_neg*(S3-S2)) + lambd_S*(sigma[3] - S3) - F3

    SA[dT1, dT2, dT3, dS1, dS3]
end