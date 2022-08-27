clc;clear;close;
%% PS
% 第一次迭代中，体现出跟踪能力
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
iterK = 10;
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
[theta, Phi, Kec, theta_u, H, G, X] = GetCoefMat(Ae, Be, Ce, P, N);
%% start simulation
% set initial value, t = -1
% initial value of predictive state/output variable
Xm_neg = 0;
Ym_neg = 0;
% initial value of real process state/output variable
Xp_neg = [0; 0];
Yp_neg = 0;
% initial input u
u_neg = 0;
% create storage vector
Xm = zeros(T, iterK);
Xp = repmat({Xp_neg}, T, iterK);
Ym = zeros(T, iterK);
Yp = zeros(T, iterK);
U = zeros(T, iterK);
% coefficient matrix
W = repmat({zeros(N,1)}, T, iterK);
% store predictive control input
predU = repmat({zeros(P, 1)}, T, iterK);
%% K = 1, t = 0
% predictive model
Xm_0 = Am * Xm_neg + Bm * u_neg;
Ym_0 = Cm * Xm_0;
% real process model
Xp_0 = A * Xp_neg + B * u_neg;
Yp_0 = C * Xp_0;
delta_xm = Xm_0 - Xm_neg;
et = Ym(1, 1) - Ym(1, 1);
% extended state vector
z = [delta_xm; et];
Z_0 = z;
deltaYr_P = GetDeltaYr(P, w, Yp_0, 15);
mat1 = -inv((X.') * Gamma * X + (G.') * Lambda * G);
mat2 = (X.') * Gamma * (theta * Z_0 + Kec * deltaYr_P - theta_u * u_neg);
W{1, 1} = mat1 * mat2;
u_0 = W{1, 1}(1);
%% k = 1, t = 1
% predictive model
Xm(1, 1) = Am * Xm_0 + Bm * u_0;
Ym(1, 1) = Cm * Xm(1, 1);
% real process model
Xp{1, 1} = A * Xp_0 + B * u_0;
Yp(1, 1) = C * Xp{1, 1};
delta_xm = Xm(1, 1) - Xm_0;
et = Ym(1, 1) - Ym(1, 1);
% extended state vector
z = [delta_xm; et];
Z = repmat({z}, T, iterK);
Z{1, 1} = z;

% YrVec = GetYr_P(1, w, Yp(1, 1), c(1), P);
% Yr_P = repmat({YrVec}, T, iterK);
% deltaYr_P = YrVec - GetYr_P(0, w, Yp(1, 1), c(1), P);
deltaYr_P = GetDeltaYr(P, w, Yp(1, 1), c(1));
% Ut = repmat(u_0, P, 1);
mat1 = -inv((X.') * Gamma * X + (G.') * Lambda * G);
mat2 = (X.') * Gamma * (theta * Z{1, 1} + Kec * deltaYr_P - theta_u * u_0);
W{1, 1} = mat1 * mat2;

predU{1, 1} = H * W{1, 1}; 
U(1, 1) = predU{1, 1}(1, 1);

% YpVec = GetYp_P(Xp(1, 1), Yp(1, 1), predU(1, 1), P, A, B, C);
% Yp_P = repmat({YpVec}, T, iterK);

%% K = 1, t = 2 : T
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

%     YrVec = GetYr_P(1, w, Yp(t, 1), c(t), P);
%     Yr_P{t, 1} = YrVec;
%     deltaYr_P = YrVec - GetYr_P(0, w, Yp(t, 1), c(t), P);
    deltaYr_P = GetDeltaYr(P, w, Yp(t, 1), c(t));
%     Ut = repmat(U(t - 1, 1), P, 1);
    mat2 = (X.') * Gamma * (theta * Z{t, 1} + Kec * deltaYr_P - theta_u * U(t - 1, 1));
    W{t, 1} = mat1 * mat2;
    predU{t, 1} = H  * W{t, 1};
    U(t, 1) = predU{t, 1}(1, 1);
%     YpVec = GetYp_P(Xp(t, 1), Yp(t, 1), predU(t, 1), P, A, B, C);
%     Yp_P{t, 1} = YpVec;
end
%% assessment of performance
IAE = zeros(iterK, 1);
for i = 1 : T
    IAE(1, 1) = IAE(1, 1) + abs(c(i) - Yp(i, 1));
end
%% plot
plot(1 : T, c, 'r');
hold on;
plot(1 : T, Yp(:,1),'b-');
hold on;
plot(1 : T, Ym(:,1), 'g+')