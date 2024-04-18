function T_m_i_plus_1 = backward_conduction(T_m_i, T_m_plus_1_i, k, A1, rho, C_p, delta_x, delta_t)
    A_avg = A1;
    conduction_const = k * A1 / delta_x;
    tau = rho * C_p * A_avg * delta_x / delta_t;

    T_m_i_plus_1 = (conduction_const * (T_m_plus_1_i - T_m_i) + tau * T_m_i ) / tau;
end

