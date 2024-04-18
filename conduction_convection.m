function T_m_i_plus_1 = conduction_convection(T_m_minus_1_i, T_m_i, T_m_plus_1_i, k, A1, h, As, rho, C_p, delta_x, delta_t)
    conduction_const = k * A1 / delta_x;
    convection_const = h * As;
    tau = rho * C_p * A1 * delta_x / delta_t;
    middle_term = conduction_const + convection_const - tau;

    T_m_i_plus_1 = (conduction_const * T_m_minus_1_i - middle_term * T_m_i + convection_const * T_m_plus_1_i) / tau;
end

