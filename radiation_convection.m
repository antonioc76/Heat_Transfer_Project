function T_m_i_plus_1 = radiation_convection(T_m_minus_1_i, T_m_i, T_m_plus_1_i, epsilon, A1, h, As, rho, C_p, delta_x, delta_t)
    sigma_male = 5.67 * 10^(-8);
    radiation_const = sigma_male * epsilon * A1;
    convection_const = h * As;
    tau = rho * C_p * A1 * delta_x / delta_t;

    T_m_i_plus_1 = (radiation_const * T_m_minus_1_i^4 - radiation_const * T_m_i^4 - convection_const * T_m_i + ...
        convection_const * T_m_plus_1_i + tau * T_m_i ) / tau;
end

