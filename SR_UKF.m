%	Author: GhasemHeydari 
%	Email: 	ghasem.heydari2@gmail.com
% 	SR-UKF: Square Root Unscented Kalman Filter for nonlinear dynamic systems
% [x, S] =  ukf(f,x,S,h,z,Qs,Rs) returns state estimate, x and state covariance, P 
% for nonlinear dynamic system (for simplicity, noises are assumed as additive):
%           x_k+1 = f(x_k) + w_k
%           z_k   = h(x_k) + v_k
% where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
%       v ~ N(0,R) meaning v is gaussian noise with covariance R
% Inputs:   f: function handle for f(x)
%           x: "a priori" state estimate
%           S: "a priori" estimated the square root of state covariance
%           h: fanction handle for h(x)
%           z: current measurement
%           Qs: process noise standard deviation
%           Rs: measurement noise standard deviation
% Output:   x: "a posteriori" state estimate
%           S: "a posteriori" square root of state covariance

function [X, S] = SR_UKF(X, S, Z, Q, R, alpha, ki, beta, f_b, Wb_nb, flagEarth, lat0, long0, h0, dt, Throttle, ds, dr, db, flagUpdate)

L       = numel(X);                                 % numer of states
m       = numel(Z);                                 % numer of measurements
lambda  = ( (alpha^2)*(L + ki) ) - L;               % scaling factor
c       = L+lambda;                                 % scaling factor
Wm      = [lambda/c  (0.5/c)+zeros(1,2*L)];         % weights for means
Wc      = Wm;
Wc(1)   = Wc(1) + (1 - alpha^2 + beta);             %weights for covariance
gamma   = sqrt(c);

%%%% PREDICT %%%%
SigmaPoints = CalcSigmaPoints(X, S, gamma);         %sigma points around x
Ls          = size(SigmaPoints, 2);
Xminus      = zeros(L, 1);
Xx          = zeros(L, Ls);
for i = 1:1:numel(SigmaPoints(1,:))
   Point    = SigmaPoints(:,i);
   Xx(:,i)  = f(Point, f_b, Wb_nb, flagEarth, lat0, long0, h0, dt);
   Xminus   = Xminus + Wm(i) * Xx(:,i); 
   
end
Y1      = Xx - Xminus(:,ones(1,Ls));
ResPer  = sqrt(Wc(2)) * Y1;

RbnF    = Body2NED([Xminus(7) Xminus(8) Xminus(9)]);
GAMMA   = [zeros(3,3) zeros(3,3); RbnF zeros(3,3); zeros(3,3) RbnF];

ResPer  = Y1 * diag(sqrt(abs(Wc)));
[~,Sm] = qr([ResPer(:,2:2*L+1) chol(GAMMA * (Q) * GAMMA' + 1e-70*eye(9))]',0);

if Wc(1)<0
    Sminus = cholupdate(Sm,ResPer(:,1),'-');
else
    Sminus = cholupdate(Sm,ResPer(:,1),'+');
end

%%%% UPDATE  %%%%
if flagUpdate == 1
    Zplus   = zeros(m,1);
    Zm      = zeros(m,(2*L+1) );
    for k = 1:(2*L+1)
        Zm(:,k) = h(Xx(:,k), Wb_nb, lat0, long0, h0, dt, Throttle, ds, dr, db, flagEarth);
        Zplus   = Zplus + Wm(k) * Zm(:,k);
    end
    Y2      = Zm - Zplus(:,ones(1,(2*L+1)));
    % ResUp   = Y2 * diag(sqrt(abs(Wc)));
    ResUp   = sqrt(Wc(2)) * Y2;
    [~,Sp] = qr([ResUp(:,2:2*L+1) chol(R)]',0);
    if Wc(1)<0
        Splus = cholupdate(Sp, ResUp(:,1),'-');
    else
        Splus = cholupdate(Sp, ResUp(:,1),'+');
    end
    
    P12     = Y1 * diag(Wc) * Y2';
    K       = P12 / Splus / Splus';
    X       = Xminus + K * (Z - Zplus);	 %state update
    
    U       = K * Splus';
    for i = 1:m
        Sminus = cholupdate(Sminus, U(:,i), '-');
    end
    S       = Sminus;
    
else
    X       = Xminus;
    S       = Sminus;
end

end

function SigmaPoints = CalcSigmaPoints(X,S,gamma)
% Sigma points around reference point
% Inputs:
%       x: reference point
%       P: covariance
%       gamma: coefficient
% Output:
%       X: Sigma points

A           = gamma * S';
Y           = X(:, ones(1, numel(X)));
SigmaPoints = [X Y+A Y-A] ;

end

function X = f(X0, dt)

t       = 0;
h       = dt;
k1      = Dynamic( t,         X0);  
k2      = Dynamic( t+0.5*h,   X0+0.5*h*k1);
k3      = Dynamic((t+0.5*h), (X0+0.5*h*k2));
k4      = Dynamic((t+h),     (X0+k3*h));
X       = X0 + (1/6)*(k1+2*k2+2*k3+k4)*h;  % main equation

end

function Z = h(X, Wb_nb, lat0, long0, h0, dt, Throttle, ds, dr, db, flagEarth)
global   param1 param2

state 	= X;  % This line is for example. you must set your formulations
Z       = [state(1), state(2), state(3), state(1) state(2) state(3)];

end




