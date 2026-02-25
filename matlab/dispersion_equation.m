function L = dispersion_equation(T, h)
    g = 9.81;
    tolerance = 0.0001;
    err = tolerance + 1;
    L = g*T^2/(2*pi);
    
    while err >= tolerance 
        f = L - ( (g / (2 * pi)) * T^2 * tanh((2 * pi / L) * h ) );
        df = 1 + (g / (2 * pi)) * T^2 * (2 * pi * h / L^2) * sech((2 * pi / L) * h)^2;
        
        L_new = L - (f/df);
        err = abs((L_new-L)/L_new);
        L = L_new;
    end
    
end