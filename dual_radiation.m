function T_m_i_plus_1 = dual_radiation(T_m_minus_1_i, T_m_i, T_m_plus_1_i, epsilon1, A1, epsilon2, A2, rho, C_p, delta_x, delta_t)
    A_avg = (A1 + A2) / 2;
    sigma = 5.67 * 10^(-8);
    radiation_const_1 = sigma * epsilon1 * A1;
    radiation_const_2 = sigma * epsilon2 * A2;
    tau = rho * C_p * A_avg * delta_x / delta_t;

    T_m_i_plus_1 = (radiation_const_1 * T_m_minus_1_i^4 - (radiation_const_1 + radiation_const_2) * T_m_i^4 + ...
        radiation_const_2 * T_m_plus_1_i^4 + tau * T_m_i) / tau;
end

