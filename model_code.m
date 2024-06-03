%% objective function for least squares fitting
function out = objfun(lat,pars,xx,yy)

% transient
Y0 = [0 0];
dt = 1;
tspan = 0:dt:24*365;
opts = odeset('RelTol',1e-10);
[T,Y] = ode45(@(t,y) tpo(t,y,lat,pars),tspan,Y0,opts);

Y0 = Y(end,:);
[T,Y] = ode45(@(t,y) tpo(t,y,lat,pars),tspan,Y0,opts);

% phase gap, converted to hours
pg = mod(abs(Y(:,1)-Y(:,2))*12/pi,24);
pg = min(pg,24-pg);

sim_y = pg(xx*24+1);
out = yy/8.7427e03 - sim_y/4.3545;

end


%% get solutions given latitude and parameters
function [T,pg] = solcurve(lat,pars)

% transient
Y0 = [0 0];
dt = 1;
tspan = 0:dt:24*365;
opts = odeset('RelTol',1e-10);
[T,Y] = ode45(@(t,y) tpo(t,y,lat,pars),tspan,Y0,opts);

Y0 = Y(end,:);
[T,Y] = ode45(@(t,y) tpo(t,y,lat,pars),tspan,Y0,opts);

% phase gap, converted to hours
pg = mod(abs(Y(:,1)-Y(:,2))*12/pi,24);
pg = min(pg,24-pg);

end


%% ODEs
function dydt = tpo(t,y,lat,pars)

ae = pars(1);
am = pars(2);
be = pars(3);
bm = pars(4);
Ae = pars(5);
Am = pars(6);
ip = pars(7);

dydt = zeros(2,1);

day = t/24; % convert t to days
% day_length (source: https://www.mathworks.com/matlabcentral/fileexchange/20390-day-length)
L = day_length(day,lat); 

% light intensity is proportional to L
eli = ae*L; % E
mli = am*L; % M

% intrinsic period depends on light intensity
De = 1/(1+exp(-(ip*(eli-1)))); % nominal model
Dm = 1/(1+exp(-(ip*(mli-1))));

we = 2*pi/(24+De);
wm = 2*pi/(24-Dm);

% coupling to light oscillator depends on light intensity
Aze = be*eli;
Azm = bm*mli;

lp = light_phase(t,L);

dydt(1) = we + Ae*sin(y(2)-y(1)) + Aze*sin(lp-y(1));
dydt(2) = wm + Am*sin(y(1)-y(2)) + Azm*sin(lp-y(2));

end


%% light phase given time and photoperiod
function out = light_phase(T,L)

out = zeros(length(T),1);

if L>0 % if L=0 out will just be 0s
    for i=1:length(T)
    
        if mod(T(i),24)<L
            out(i) = pi/L * mod(T(i),24);
        else
            out(i) = pi/(24-L) * (mod(T(i),24)-L) + pi;
        end
    end
end
end