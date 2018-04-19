classdef SimulationData < handle
    properties
        dt = 0.01;
        t = 0:0.01:500;

        Fd = [];
        Ft = [];
        m = [];
        theta = [];
        z = [];
        vz = [];
        az = [];
        x = [];
        vx = [];
        ax = [];
        
        I = 2;
    end
    
    methods
        function o = SimulationData(dt, tMax)
            o.dt = dt;
            o.t = 0:dt:tMax;

            o.Fd = zeros(length(o.t), 1);
            o.Ft = zeros(length(o.t), 1);
            o.m = zeros(length(o.t), 1);
            o.theta = zeros(length(o.t), 1);
            o.z = zeros(length(o.t), 1);
            o.vz = zeros(length(o.t), 1);
            o.az = zeros(length(o.t), 1);
            o.x = zeros(length(o.t), 1);
            o.vx = zeros(length(o.t), 1);
            o.ax = zeros(length(o.t), 1);
        end
        
        function SetInitialConditions(o, mFull)
            o.z(1) = 0;
            o.vz(1) = 0;
            o.az(1) = -9.81;
            o.x(1) = 0;
            o.vx(1) = 0;
            o.ax(1) = 0;
            o.m(1) = mFull; % Start full
            o.Fd(1) = 0;
            o.Ft(1) = 0;
        end
        
        function ResetFrame(o)
            o.I = 2;
        end
        
        function Thrust(o, mdot, vjet, Pe, Pa, Ae, N)
            o.m(o.I) = o.m(o.I-1) - mdot * o.dt * N;
            o.Ft(o.I) = (mdot * vjet + (Pe - Pa) * Ae) * N;
        end
        
        function NoThrust(o)
            o.m(o.I) = o.m(o.I-1);
            o.Ft(o.I) = 0;
        end
        
        function relWind = RelativeWind(o)
            zpart = abs(o.vz(o.I-1));
            xpart = abs(o.vx(o.I-1));
            relWind = abs(cos(o.theta(o.I-1)) * zpart) ...
                        + abs(sin(o.theta(o.I-1)) * xpart);
        end
        
        function UpdateVZ(o, dv_axial, g)
            o.vz(o.I) = o.vz(o.I-1) + dv_axial * cos(o.theta(o.I-1)) - g*o.dt;
        end
        
        function UpdateAngle(o, dv_axial, thetaWindow, thetaAngle)
            % Migrate to the angle over first few seconds
            if o.t(o.I) < thetaWindow
                o.theta(o.I) = thetaAngle * (o.t(o.I) / thetaWindow);
                o.vx(o.I) = tan(o.theta(o.I-1)) * o.vz(o.I-1);
            else
                o.vx(o.I) = o.vx(o.I-1) + dv_axial * sin(o.theta(o.I-1));
                o.theta(o.I) = atan(o.vx(o.I)/o.vz(o.I));
            end
        end
        
        function CalculateNewXandZ(o)
            % Net acceleration
            dvx = o.vx(o.I) - o.vx(o.I-1);
            o.ax(o.I) = dvx / o.dt; 
            o.x(o.I) = o.x(o.I-1) + 0.5 * dvx * o.dt + o.vx(o.I-1) * o.dt;
            dvz = o.vz(o.I) - o.vz(o.I-1);
            o.az(o.I) = dvz / o.dt;
            o.z(o.I) = o.z(o.I-1) + 0.5 * dvz * o.dt + o.vz(o.I-1) * o.dt;
        end
        
        function NextFrame(o)
            o.I = o.I + 1;
        end
        
        function mass = RemainingMass(o, m0)
            mass = o.m(o.I) - m0;
        end
        
        function v = velocity(o, offset)
            v = (o.vz(o.I + offset).^2 + o.vx(o.I + offset).^2).^0.5;
        end
        
    end
end