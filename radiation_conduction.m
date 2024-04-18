function T_m_i_plus_1 = radiation_conduction(T_m_minus_1_i, T_m_i, T_m_plus_1_i, epsilon, A1, k, A2, rho, C_p, delta_x, delta_t)
    A_avg = (A1 + A2) / 2;
    sigma = 5.67 * 10^(-8);
    radiation_const = sigma * epsilon * A1;
    conduction_const = k * A2 / delta_x;
    tau = rho * C_p * A_avg * delta_x / delta_t;

    T_m_i_plus_1 = (radiation_const * T_m_minus_1_i^4 - radiation_const * T_m_i^4 - conduction_const * T_m_i + ...
        conduction_const * T_m_plus_1_i + tau * T_m_i ) / tau;
end

