

classdef Atmosphere < handle
    properties
        T0_arr = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65];
        p0_arr = [101325, 22632.1, 5474.89, 868.02, 110.91, 66.94, 3.96];
        B_arr = [-0.0065, 0, 0.001, 0.0028, 0.0, -0.0028, -0.002];
        z_arr = [0, 11000, 20000, 32000, 47000, 51000, 71000];
        g0 = 9.81;
        R0 = 6370;
        pMin = 0.000000001;
        Tmin = 3;
        R = 287;
    end
    
    methods
        function [T0, p0, B, z0] = Dynamics(o, z)
            idx = max(find(z >= o.z_arr));
            T0 = o.T0_arr(idx);
            p0 = o.p0_arr(idx);
            B = o.B_arr(idx);
            z0 = o.z_arr(idx);
        end
        
        function [p, rho, T, g] = ParametersAtAltitude(o, z)
            % Get the dynamics variables
            [T0, p0, B, z0] = o.Dynamics(z);

            % Now calculate gravity
            g = o.g0 * (o.R0 / (o.R0 + z/1000)).^2;
            
            % Above critical alt?
            if z > 180000
                p = o.pMin;
                T = o.Tmin;
            else
                % First calculate temperature
                T = T0 + B * (z - z0);
                if T < o.Tmin || ~isreal(T)
                    T = o.Tmin;
                end
                
                % Now calculate pressure
                if B == 0
                    p = p0 * exp(-g * (z - z0) / (o.R * T0));
                else
                    p = p0 * (T0 / T).^(g / (o.R * B));
                end

                % Keep to a min
                if p < o.pMin || ~isreal(p)
                    p = o.pMin;
                end
            end

            % Then calculate density
            rho = p / (o.R * T);
        end
        
    end
    
end
