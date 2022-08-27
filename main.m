clc;clear;close;
%% model information
% the actual process model -- ss
A = [1.582 -0.5916; 1 0];
B = [1; 0];
C = [1.69 1.419];
D = 0;
% the predictive model
[Am, Bm, Cm, ~] = tf2ss(7.31, [1 -0.92]);

%% configuration
% total operating time
T = 100;
% the number of iteration
k = 10;
% the setpoint: piecewise function
c = ones(1, T);
c(1, 1 : 49) = 15;
c(1, 50 : 100) = 30;
% the prediction horizon
P = 3;
% the number of base function
N = 2;
% the smoothing factor
w = 0.8;
% the weighting matrix Gamma
base_vector = ones(1, P);
% gamma = 1;
% Gamma = diag(base_vector * gamma);

gamma = diag(2,1);
Gamma = GetDiagMarix(gamma, P);
% the weighting matrix Lambda
lambda = 5;
Lambda = diag(base_vector * lambda);
% the weighting matrix Beta
beta = 0.01;
Beta = diag(base_vector * beta);
% the parameter tao
tao = 1;

%% extended state-space model's coefficient matrix
Ae = [Am, 0; Cm * Am, 1];
Be = [Bm; Cm * Bm];
Ce = [0; -1];

%% calculate predicton model's coefficient matrix
[row_Ae, col_Ae] = size(Ae);

baseMat = eye(row_Ae * P, col_Ae);
for ii = 2 : P
    baseMat((ii - 1) * row_Ae + 1: ii * row_Ae, :) ...
        = baseMat((ii - 2) * row_Ae + 1 : (ii - 1) * row_Ae, :) * Ae;
end
% matrix theta
theta = baseMat * Ae;

% matrix Phi
[row_Be, col_Be] = size(Be);
Phi = zeros(row_Be * P, col_Be * P);
for ii = 1 : P
    Phi((ii - 1) * row_Be + 1 : row_Be * P, (ii - 1) * col_Be + 1 : ii * col_Be) ...
        = baseMat(1 : (P - ii + 1) * row_Be, :) * Be;
end
% matrix kec
[row_Ce, col_Ce] = size(Ce);
kec = zeros(row_Ce * P, col_Ce * P);
for ii = 1 : P
    kec((ii - 1) * row_Ce + 1 : end, (ii - 1) * col_Ce + 1 : ii * col_Ce) ...
        = baseMat(1 : (P - ii + 1) * row_Ce, :) * Ce;
end

% matrix theta_u
theta_u = Phi(:, 1 : col_Be);

% matrix H
H = zeros(P, N);
for ii = 0 : (P - 1)
    H(ii + 1, :) = baseFunc(ii, N);
end

% matrix G
G = zeros(P, N);
G(1, :) = baseFunc(0, N);
for ii = 1 : (P - 1)
    G(ii + 1, :) = baseFunc(ii, N) - baseFunc(ii - 1, N);
end

% X matrix
X = Phi * G;

%% start simulation
% set initial value
% initial value of predictive state/output variable
xm_0 = 0;
ym_0 = 0;
% initial value of real process state/output variable
xp_0 = [0; 0];
yp_0 = 0;
% initial input u
u_0 = 0.2;

% storage vector
Xm = zeros(T, k);
Xp = repmat({xp_0}, T, k);
Ym = zeros(T, k);
Yp = zeros(T, k);
U = zeros(T, k);
% coefficient matrix
W = repmat({zeros(N,1)}, T, k);
%% k = 1, t = 1
% predictive model
Xm(1, 1) = Am * xm_0 + Bm * u_0;
Ym(1, 1) = Cm * Xm(1, 1);

% real process model
Xp{1, 1} = A * xp_0 + B * u_0;
Yp(1, 1) = C * Xp{1, 1};

delta_xm = Xm(1, 1) - xm_0;
et = Ym(1, 1) - Ym(1, 1);
% extended state vector
z = [delta_xm; et];
Z = repmat({z}, T, k);
Z{1, 1} = z;

deltaYr_P = GetDeltaYr(P, w, Yp(1, 1), c(1));

ev = [zeros(size(delta_xm)); et];
Ek = repmat(ev, P, 1);

Ut = repmat(u_0, P, 1);

mat1 = -inv((X.') * Gamma * X + (G.') * Lambda * G);
mat2 = (X.') * Gamma * (theta * Z{1, 1} + kec * deltaYr_P + Ek - theta_u * u_0) ...
    - (G.') * Lambda * Ut;
W{1, 1} = mat1 * mat2;

baseVec = baseFunc(0, N);
U(1, 1) = baseVec * W{1, 1};

%% k = 1, t = 2 : T
for t = 2 : T
    % predictive model
    Xm(t, 1) = Am * Xm(t - 1, 1) + Bm * U(t - 1, 1);
    Ym(t, 1) = Cm * Xm(t, 1);

    % real process model
    Xp{t, 1} = A * Xp{t - 1, 1} + B * U(t - 1, 1);
    Yp(t, 1) = C * Xp{t, 1};

    delta_xm = Xm(t, 1) - Xm(t - 1, 1);
    et = Ym(t, 1) - Ym(t, 1);
    % extended state vector
    z = [delta_xm; et];
    Z{t, 1} = z;

    deltaYr_P = GetDeltaYr(P, w, Yp(t, 1), c(t));

    ev = [zeros(size(delta_xm)); et];
    Ek = repmat(ev, P, 1);

    Ut = repmat(U(t - 1, 1), P, 1);

    mat2 = (X.') * Gamma * (theta * Z{t, 1} + kec * deltaYr_P + Ek - theta_u * U(t - 1, 1)) ...
        - (G.') * Lambda * Ut;
    W{t, 1} = mat1 * mat2;
    U(t, 1) = baseVec * W{t, 1};
end

plot(1 : T, c, 'r');
hold on;
plot(1 : T, Yp(:,1),'b-');
hold on;
plot(1 : T, Ym(:,1), 'g+')