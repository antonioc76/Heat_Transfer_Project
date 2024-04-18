function T_m_i_plus_1 = conduction_radiation(T_m_minus_1_i, T_m_i, T_m_plus_1_i, k, A1, epsilon, A2, rho, C_p, delta_x, delta_t)
    A_avg = (A1 + A2) / 2;
    sigma = 5.67 * 10^(-8);
    conduction_const = k * A1 / delta_x;
    radiation_const = sigma * epsilon * A2;
    tau = rho * C_p * A_avg * delta_x / delta_t;

    T_m_i_plus_1 = (conduction_const * T_m_minus_1_i + ... 
        -1 * conduction_const * T_m_i - radiation_const * T_m_i^4 + tau * T_m_i + ...
        radiation_const * T_m_plus_1_i^4) / tau;
end

