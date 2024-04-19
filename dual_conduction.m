function T_m_i_plus_1 = dual_conduction(T_m_minus_1_i, T_m_i, T_m_plus_1_i, k1, A1, k2, A2, rho1, rho2, C_p_1, C_p_2, delta_x, delta_t)
    A_avg = (A1 + A2) / 2;
    rho_avg = (rho1 + rho2) / 2;
    C_p_avg = (C_p_1 + C_p_2) / 2;
    conduction_const_1 = k1 * A1 / delta_x;
    conduction_const_2 = k2 * A2 / delta_x;
    tau = rho_avg * C_p_avg * A_avg * delta_x / delta_t;
    middle_term = conduction_const_1 + conduction_const_2 - tau;

    T_m_i_plus_1 = (conduction_const_1 * T_m_minus_1_i - middle_term * T_m_i + conduction_const_2 * T_m_plus_1_i ) / tau;
end

