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
% w = 0.8;
w = 0.65;
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
et = 0;
% extended state vector
z = [delta_xm; et];
Z_0 = z;

deltaYr_P = GetDeltaYr(P, w, Yp_0, 15);

Ut = [0; zeros(P - 1, 1)];

mat1 = -inv((X.') * Gamma * X + (G.') * Lambda * G);
mat2 = (X.') * Gamma * (theta * Z_0 + Kec * deltaYr_P - theta_u * u_neg);
mat3 = (G.') * Lambda * Ut;

W_0 = mat1 * (mat2 - mat3);
predU_0 = repmat({H * W_0}, iterK, 1);
u_0 = W_0(1);

%% k = 1, t = 1
% predictive model
Xm(1, 1) = Am * Xm_0 + Bm * u_0;
Ym(1, 1) = Cm * Xm(1, 1);

% real process model
Xp{1, 1} = A * Xp_0 + B * u_0;
Yp(1, 1) = C * Xp{1, 1};

delta_xm = Xm(1, 1) - Xm_0;
% extended state vector
z = [delta_xm; et];
Z = repmat({z}, T, iterK);
Z{1, 1} = z;

deltaYr_P = GetDeltaYr(P, w, Yp(1, 1), c(1));

% Ut(1) = u_0;
Ut = [u_0; zeros(P - 1, 1)];

mat1 = -inv((X.') * Gamma * X + (G.') * Lambda * G);
mat2 = (X.') * Gamma * (theta * Z{1, 1} + Kec * deltaYr_P - theta_u * u_0);
mat3 = (G.') * Lambda * Ut;
W{1, 1} = mat1 * (mat2 - mat3);

predU{1, 1} = H * W{1, 1};
U(1, 1) = predU{1, 1}(1, 1);
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

    deltaYr_P = GetDeltaYr(P, w, Yp(t, 1), c(t));
    %Ut(1) = U(t - 1, 1);
    Ut = [U(t - 1, 1); zeros(P - 1, 1)];

    mat2 = (X.') * Gamma * (theta * Z{t, 1} + Kec * deltaYr_P - theta_u * U(t - 1, 1));
    mat3 = (G.') * Lambda * Ut;

    W{t, 1} = mat1 * (mat2 - mat3);
    predU{t, 1} = H  * W{t, 1};
    U(t, 1) = predU{t, 1}(1, 1);
end
%% assessment of performance: the first circle
IAE = zeros(iterK, 1);
for i = 1 : T
    IAE(1, 1) = IAE(1, 1) + abs(c(i) - Yp(i, 1));
end
IAE(1)

IACS = zeros(iterK, 1);
for i = 1 : T
    IACS(1, 1) = IACS(1, 1) + abs(U(i, 1));
end
IACS(1)

%% k = 2 : iterK
for k = 2 : iterK
    % t = 0
    deltaYr_P = GetDeltaYr(P, w, Yp_0, 15);
    Ut = [0; zeros(P - 1, 1)];
    Uk = predU_0{k - 1};
    mat2 = (X.') * Gamma * (theta * Z_0 + Kec * deltaYr_P - theta_u * u_neg);
    mat3 = (G.') * Lambda * Ut;
    mat4 = (H.') * Beta * Uk;

    W_0 = mat1 * (mat2 - mat3- mat4);
    predU_0{k} = H * W_0;
    u_0 = predU_0{k}(1);

    % t = 1
    % predictive model
    Xm(1, k) = Am * Xm_0 + Bm * u_0;
    Ym(1, k) = Cm * Xm(1, k);
    % real process model
    Xp{1, k} = A * Xp_0 + B * u_0;
    Yp(1, k) = C * Xp{1, k};
    delta_xm = Xm(1, k) - Xm_0;
    et = 0;
    % extended state vector
    z = [delta_xm; et];
    Z{1, k} = z;
    deltaYr_P = GetDeltaYr(P, w, Yp(1, k), c(1));
    Ut = [u_0; zeros(P - 1, 1)];
    Uk = predU{1, k - 1};

    mat2 = (X.') * Gamma * (theta * Z{1, k} + Kec * deltaYr_P - theta_u * u_0);
    mat3 = (G.') * Lambda * Ut;
    mat4 = (H.') * Beta * Uk;
    W{1, k} = mat1 * (mat2 -mat3 - mat4);

    predU{1, k} = H * W{1, k};
    U(1, k) = predU{1, k}(1);

    for t = 2 : T
        % predictive model
        Xm(t, k) = Am * Xm(t - 1, k) + Bm * U(t - 1, k);
        Ym(t, k) = Cm * Xm(t, k);

        % real process model
        Xp{t, k} = A * Xp{t - 1, k} + B * U(t - 1, k);
        Yp(t, k) = C * Xp{t, k};

        delta_xm = Xm(t, k) - Xm(t - 1, k);
        et = Ym(t, 1) - Ym(t, 1);
        % extended state vector
        z = [delta_xm; et];
        Z{t, k} = z;

        deltaYr_P = GetDeltaYr(P, w, Yp(t, k), c(t));

        Ut = [U(t - 1, k); zeros(P - 1, 1)];
        Uk = predU{t, k - 1};

        mat2 = (X.') * Gamma * (theta * Z{t, k} + Kec * deltaYr_P - theta_u * U(t - 1, k));
        mat3 = (G.') * Lambda * Ut;
        mat4 = (H.') * Beta * Uk;
        W{t, k} = mat1 * (mat2 -mat3 - mat4);
        predU{t, k} = H  * W{t, k};
        U(t, k) = predU{t, k}(1);
    end
end

%% assessment of performance: the last circle
IAE = zeros(iterK, 1);
for i = 1 : T
    IAE(1, 1) = IAE(1, 1) + abs(c(i) - Yp(i, iterK));
end
IAE(1)

IACS = zeros(iterK, 1);
for i = 1 : T
    IACS(1, 1) = IACS(1, 1) + abs(U(i, iterK));
end
IACS(1)

%% plot
plot(1 : T, c, 'r');
hold on;
plot(1 : T, Yp(:,iterK),'b-');
hold on;
plot(1 : T, Ym(:,iterK), 'g+')
