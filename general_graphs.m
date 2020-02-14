%% tracking parameters
%{
load('tracking_vars','x','xd','zx','y','yd','zy');
ref_x = x;
ref_y = y;
clear x xd zx y yd zy
%}

%% General parameters
delta = 0.01;       % time step
w = 0.08/2;         % half foot size = 0.04

%{
x(1) = 0.0;         % xc
xd(1) = 0.0;        % xc_dot
zx(1) = 0.0;
% dfx(1) = 0.0;

y(1) = 0;           % yc
yd(1) = 0;          % yc_dot
zy(1) = 0;
% dfy(1) = 0.0;
%}

% %{
shifting_x = 100;
shifting_y = 100;
x(1) = shifting_x;         % xc
xd(1) = 0.0;               % xc_dot
zx(1) = shifting_x;
x_int(1) = 0;

y(1) = shifting_y;         % yc
yd(1) = 0;                 % yc_dot
zy(1) = shifting_y;
y_int(1) = 0;
%}

warning('off','all');

%% Compute constraints
f1_y = -foot_distance_y;
f2_y = foot_distance_y;

% additionalFirstStepDuration = 10;           % why wait for a while at first???
additionalFirstStepDuration = 0;

% one step duration                         % why zeros at first step???
fs_sequence_x = zeros(S+D+additionalFirstStepDuration,1);
fs_sequence_y = zeros(S+D+additionalFirstStepDuration,1);

% fs_matrix is the preassigned footsteps matrix
for i = 1:12
    f1_x = fs_matrix(i,1);
    f2_x = fs_matrix(i+1,1);
    
    f1_y = fs_matrix(i,2);
    f2_y = fs_matrix(i+1,2);
    
    % [old sequence; single support has same value; 
    %  double support has an interpolated value (linear)]
    fs_sequence_x = [fs_sequence_x; ones(S,1) * f1_x; f1_x + (1:D)'*(f2_x-f1_x)/D];
    fs_sequence_y = [fs_sequence_y; ones(S,1) * f1_y; f1_y + (1:D)'*(f2_y-f1_y)/D];
end

% why this deletion???
fs_sequence_x(1) = [];
fs_sequence_y(1) = [];

% ZMP constraints are the boundaries of the support polygon at time step
%{
zx_min = fs_sequence_x - w;
zx_max = fs_sequence_x + w;

zy_min = fs_sequence_y - w;
zy_max = fs_sequence_y + w;
%}

% %{
zx_min = fs_sequence_x - w + shifting_x;
zx_max = fs_sequence_x + w + shifting_x;

zy_min = fs_sequence_y - w + shifting_y;
zy_max = fs_sequence_y + w + shifting_y;
%}

% make support polygon of stand still contain both feet
zy_min(1:S+D+additionalFirstStepDuration-1) = zy_min(1:S+D+additionalFirstStepDuration-1) - foot_distance_y;
zy_max(1:S+D+additionalFirstStepDuration-1) = zy_max(1:S+D+additionalFirstStepDuration-1) + foot_distance_y;


%% Compute matrices
%{
% p = ones(N,1);
% P = delta*tril(ones(N,N));
% A = [P;-P];

% Vu = [];
% Vs = [];
% 
% for i = 1:N
% 
%     A_power = eye(3);
%     ch = cosh(omega*delta);
%     sh = sinh(omega*delta);
%     A_upd_j = [ch      , sh/omega, 1-ch; 
%                omega*sh, ch      , -omega*sh; 
%                0       , 0       , 1];
%     for j = 1:i
%         A_power = A_power * A_upd_j;
%     end
%     
%     Ps_newline = A_power(1,:);
%     Vs_newline = A_power(2,:);
%     
%     Pu_newline = zeros(1,N);
%     Vu_newline = zeros(1,N);
%     
%     for j = 1:i
%         A_power = eye(3);
%         A_upd_n = A_upd_j;
%         for n = j+1:i
%             A_power = A_power * A_upd_n;
%         end
%         
%         B_upd_j = [delta-sh/omega; 1-ch; delta];
%         coeff = A_power*B_upd_j;
%         
%         Pu_newline(1,j) = coeff(1);
%         Vu_newline(1,j) = coeff(2);
%     end
%     
%     Vu = [Vu; Vu_newline];
%     Vs = [Vs; Vs_newline];
% end
%}

%% Compute stability constraint
%{
% (1/omega) * (1-lambda)/(1-lambda^N) * [1 lambda ... lambda^(N-1)]
% Aeq = (1-exp(-omega*delta))/omega * exp(-omega*delta*(0:N-1)) / (1-exp(-omega*delta*N));
% phi_x = eye(3);
% phi_y = eye(3);
% zd_x(1) = 0;
% zd_y(1) = 0;
%}

%% compute vector ...
%{
% ch = cosh(omega*delta);     % ( e^x + e^-x )/2
% sh = sinh(omega*delta);     % ( e^x - e^-x )/2
% A_upd = [ch      , sh/omega, 1-ch;
%          omega*sh, ch      , -omega*sh;
%          0       , 0       , 1];
% 
% dist_effect_xz_vector = zeros(N,1);
% dist_effect = 0;
% D = [0; 0; 0.02];
% for i=1:N
%     for j=0:i
%         dist_effect = dist_effect + A_upd^j*D;
%     end
%     dist_effect_xz = dist_effect(3);
%     dist_effect_xz_vector(i) = dist_effect_xz;
% end
% disp(dist_effect_xz_vector)
%}

%% repeated computations

% system discretization
ch = cosh(omega*delta);     % ( e^x + e^-x )/2
sh = sinh(omega*delta);     % ( e^x - e^-x )/2

A_upd = [ch      , sh/omega, 1-ch;
         omega*sh, ch      , -omega*sh;
         0       , 0       , 1];
B_upd = [delta-sh/omega;
         1-ch;
         delta];

Ts=0.01;

% %{
% A_matrix = [0      , 1, 0;
%             omega^2, 0, -omega^2;
%             0      , 0, 0];
% B_matrix = [0;
%             0;
%             1];

% k_lqr = dlqr(A_upd, B_upd, diag([100,100,100]), 1);     % not working
% k_lqr = dlqr(A_upd, B_upd, diag([1,1,1]), 1);     % not working
% k_lqr = dlqr(A_upd, B_upd, diag([0.1,0.1,0.1]), 1);     % working
% k_lqr = dlqr(A_upd, B_upd, diag([0,10^-10,10^-10]), 10^-10);    % better
k_lqr = dlqr(A_upd, B_upd, diag([0,0.1,0]), 1); 

phi = A_upd - B_upd * k_lqr;
% disp(eig(phi))
%}

%{
C_lqi = [1 0 0];

sys_d_lqi = ss(A_upd, B_upd, C_lqi, 0, Ts);
K_lqi = lqi(sys_d_lqi,diag([100,100,100,100]), 1);

integr = c2d(ss(tf([1],[1 0])),Ts);
%}
         
%{
A_upd_d = [ch      , sh/omega , 1-ch     , 0;
           omega*sh, ch       , -omega*sh, 0;
           0       , 0        , 1        , 0;
           1-ch    , -sh/omega, ch-1     , 1];
B_upd_d = [delta-sh/omega;
           1-ch;
           delta;
           sh/omega];
       
% k_lqr_d = dlqr(A_upd_d, B_upd_d, diag([0,0,0,100]), 10);
k_lqr_d = dlqr(A_upd_d, B_upd_d, diag([0,0,0,10000]), 10000);
% k_lqr_d = [-1.001   -1    1   0];
phi_d = A_upd_d - B_upd_d * k_lqr_d;
%minimize y component only so that it can move
%}

% %{
% disturbance

D1 = 0*[sh/omega;ch-1;0];                      % D = [1;0;0]
D2 =[ch/(omega^2)-1/(omega^2);sh/omega;0];    % D = [0;1;0]
D3 = [delta-sh/omega; 1-ch; delta];          % D = [0;0;1]

% %{
% D3 = 1.5*D3;
% D2 = 0.2*D2;
% D1 = 0.06*D1*0;
% D3 = 0.5*D3;
% D2 = 0.1*D2;
% D1 = 0.06*D1;
D3 = 0.12*D3;
D2 = 0.12*D2;
D1 = 0.03*D1*0.3;
d = 1;
%}

%{
d = 0;
D3 = 0*D3;
D2 = 0*D2;
D1 = 0*D1;
%}

D = 0*D1 + D2 + D3;
d_vector = ones(N,1) * d;

% my H,Q,S,T,D
% my H,Q,S,T,D
R = 10;
R_bar = diag(ones(N,1)*R);
num_states = size(A_upd, 1);
% Q = eye(num_states);
Q = eye(num_states)*0;
Q_bar = [];
for j = 1:N
    Q_bar = blkdiag(Q_bar, Q);
end

S_bar = zeros(N*num_states, N);
for j = 1:N
    for k = 1:j
        % num_states*(j-1)+1 = j*num_states-2
        S_bar(num_states*(j-1)+1:j*num_states,k) = A_upd^(j-k)*B_upd;
    end
end

T_bar = zeros(N*num_states, num_states);
for j = 1:N
    T_bar(num_states*(j-1)+1:j*num_states, :) = A_upd^j;
end

D_bar = zeros(N*num_states, N);
for j = 1:N
    for k = 1:j
        D_bar(num_states*(j-1)+1:j*num_states,k) = A_upd^(j-k)*D;
    end
end

H = 2*(R_bar + S_bar'*Q_bar*S_bar);
F = 2*S_bar'*Q_bar*T_bar;
L = 2*S_bar'*Q_bar*D_bar;

C = [0 0 1];
S_bar_y = zeros(N, N);
for j = 1:N
    for k = 1:j
        S_bar_y(j,k) = C*A_upd^(j-k)*B_upd;
    end
end

% %{
F_bar = zeros(N, N);
for j = 1:N
    for k = 1:j
        if j==k
            F_bar(j,k) = 1;
        else
            F_bar(j,k) = -k_lqr*phi^(j-k-1)*B_upd;
        end
    end
end
%}

T_bar_y = zeros(N, num_states);
for j = 1:N
    T_bar_y(j, :) = C*A_upd^j;
end

% %{
T_F_bar = zeros(N, num_states);
for j = 1:N
    T_F_bar(j, :) = -k_lqr*phi^(j-1);
end
%}

% K = [1 1 1];
K = [0 0 1];
D_bar_y = zeros(N, N);
D_bar_y_abs = zeros(N, N);

%{
for j = 1:N-1
    for k = 1:j
%         D_bar_y_abs(j+1,k) = abs(K*A_upd^(j-k)*D);
        D_bar_y(j+1,k) = K*A_upd^(j-k)*D;
    end
end
%}

% %{
for j = 1:N-1
    for k = 1:j
        D_bar_y_abs(j+1,k) = abs(K*phi^(j-k)*D);
%         D_bar_y(j+1,k) = K*phi^(j-k)*D;
    end
end
%}

% constraints_restrict = D_bar_y*d_vector;
constraints_restrict_abs = D_bar_y_abs*d_vector;
constraints_restrict = constraints_restrict_abs;
% disp([constraints_restrict constraints_restrict_abs])

%}

%{
d = 0;
D = [0; 0; 0; 0];
d_vector = ones(N,1) * d;

% my H,Q,S,T,D
R = 10;
R_bar = diag(ones(N,1)*R);
num_states = size(A_upd_d, 1);
Q = eye(num_states)*0;
Q_bar = [];
for j = 1:N
    Q_bar = blkdiag(Q_bar, Q);
end

S_bar = zeros(N*num_states, N);
for j = 1:N
    for k = 1:j
        S_bar(num_states*(j-1)+1:j*num_states,k) = A_upd_d^(j-k)*B_upd_d;
    end
end

T_bar = zeros(N*num_states, num_states);
for j = 1:N
    T_bar(num_states*(j-1)+1:j*num_states, :) = A_upd_d^j;
end

D_bar = zeros(N*num_states, N);
for j = 1:N
    for k = 1:j
        D_bar(num_states*(j-1)+1:j*num_states,k) = A_upd_d^(j-k)*D;
    end
end

H = 2*(R_bar + S_bar'*Q_bar*S_bar);
F = 2*S_bar'*Q_bar*T_bar;
L = 2*S_bar'*Q_bar*D_bar;

C_d = [0 0 1 0];
S_bar_y = zeros(N, N);
for j = 1:N
    for k = 1:j
        S_bar_y(j,k) = C_d*A_upd_d^(j-k)*B_upd_d;
    end
end

F_bar = zeros(N, N);
for j = 1:N
    for k = 1:j
        if j==k
            F_bar(j,k) = 1;
        else
            F_bar(j,k) = -k_lqr_d*phi_d^(j-k-1)*B_upd_d;
        end
    end
end

T_bar_y = zeros(N, num_states);
for j = 1:N
    T_bar_y(j, :) = C_d*A_upd_d^j;
end

T_F_bar = zeros(N, num_states);
for j = 1:N
    T_F_bar(j, :) = -k_lqr_d*A_upd_d^(j-1);
%     T_F_bar(j, :) = k_lqr*A_upd^(j-1);
end

% K = [1 1 1];
K_d = [0 0 1 0];
D_bar_y = zeros(N, N);
% D_bar_y_abs = zeros(N, N);

for j = 1:N-1
    for k = 1:j
%         D_bar_y_abs(j+1,k) = abs(K_d*A_upd_d^(j-k)*D);
        D_bar_y(j+1,k) = K_d*A_upd_d^(j-k)*D;
    end
end

constraints_restrict = D_bar_y*d_vector;
%}

%% Solve
for i = 1:300

%     Q = 0;
    P = delta*tril(ones(N,N));      % tril: Lower triangular part of matrix
%     A = [P;-P];

    % from eqn(7) in the paper: X_z_k+1 = x_z_k + P*X_z_dot_k
    % X_z_k+1 as [x1 ... xN]', x_z_k as x0, P is S_bar since X_z_dot_k is
    % the input sequence
    % he use Q instead of Q_bar!!!
%     H = Q*(P'*P) + eye(N);
    
    % f_x = x0'*F', F'=T_bar'*Q_bar*S_bar
    % T_bar' = [1 ... 1], P is S_bar, he use Q instead of Q_bar!!!
    % (- fs_sequence_x(i:i+N-1) + zx(i)) should be x0'!!!
%     f_x = Q*P'*(- fs_sequence_x(i:i+N-1) + zx(i));
%     f_y = Q*P'*(- fs_sequence_y(i:i+N-1) + zy(i));
    
%     disp([- fs_sequence_x(i:i+N-1), ones(N,1)*zx(i)])

% %{
    f_x = F*[x(i); xd(i); zx(i)] + L*d_vector;
    f_y = F*[y(i); yd(i); zy(i)] + L*d_vector;
%}
%{
    f_x = F*[x(i); xd(i); zx(i); dfx(i)] + L*d_vector;
    f_y = F*[y(i); yd(i); zy(i); dfx(i)] + L*d_vector;
%}    
    
%     b_x = [ zx_max(i:i+N-1) - T_bar_y*[x(i); xd(i); zx(i)] - constraints_restrict;...
%            -zx_min(i:i+N-1) + T_bar_y*[x(i); xd(i); zx(i)] - constraints_restrict];
%     b_y = [ zy_max(i:i+N-1) - T_bar_y*[y(i); yd(i); zy(i)] - constraints_restrict;...
%            -zy_min(i:i+N-1) + T_bar_y*[y(i); yd(i); zy(i)] - constraints_restrict];

%{
    Zx = [zx_max(i:i+N-1) - T_bar_y*[x(i); xd(i); zx(i)], ...
          zx_min(i:i+N-1) - T_bar_y*[x(i); xd(i); zx(i)]]';
    Zy = [zy_max(i:i+N-1) - T_bar_y*[y(i); yd(i); zy(i)], ...
          zy_min(i:i+N-1) - T_bar_y*[y(i); yd(i); zy(i)]]';
    CR = [constraints_restrict -constraints_restrict]';
    
%     Polyhedron('lb', [-0.04;-0.22;-inf], 'ub', [0.04;0.22;inf],'Ae',[0 0 1],'be',5)
%     Polyhedron('lb', [c(2);v(2);-inf], 'ub', [c(1);v(1);inf],'Ae',[0 0 1],'be',j)
    
    for j = 1:size(Zx,2)
        cx(1) = max(Zx(:, j));
        cx(2) = min(Zx(:, j));
        cy(1) = max(Zy(:, j));
        cy(2) = min(Zy(:, j));
        v(1) = max(CR(:, j));
        v(2) = min(CR(:, j));
%         ZZ(j) = Polyhedron([j cx(1);j cx(2)]);
%         CRR(j) = Polyhedron([j v(1);j v(2)]);
        CRR3(j) = Polyhedron('lb', [v(2);v(2);-inf], 'ub', [v(1);v(1);inf],'Ae',[0 0 1],'be',j);
        square(j) = Polyhedron('lb', [cx(2);cy(2);-inf], 'ub', [cx(1);cy(1);inf],'Ae',[0 0 1],'be',j);
%         pont(j) = ZZ(j) - CRR(j);
        pont3(j) = square(j) - CRR3(j);
    end
    
%     U1 = PolyUnion(ZZ);
%     U2 = PolyUnion(CRR);
%     Uscr = PolyUnion(CRR3);
%     U3 = PolyUnion(pont);
%     Us12 = PolyUnion(square);
    Us12cr = PolyUnion(pont3);
%     plot(Us12); figure;
    plot(Us12cr)
%     view(0, 90)
    axis([-inf inf -inf inf 0 100])
    break;
    
%     b_x = [ zx_max(i:i+N-1) - T_bar_y*[x(i); xd(i); zx(i)] ...
%            -zx_min(i:i+N-1) + T_bar_y*[x(i); xd(i); zx(i)]]
%     b_x = [ max(unique(Us12.Set(100).V(:,1))) ; ...
%            -min(unique(Us12.Set(100).V(:,1))) ]

    U_num = Us12cr.Num;
    for j = 1:U_num
        b_x(j,1) = max(unique(Us12cr.Set(j).V(:,1)));
    end
    for j = 1:U_num
        b_x(U_num+j,1) = -min(unique(Us12cr.Set(j).V(:,1)));
    end
%     disp([b_x(1:U_num) b_x(U_num+1:2*U_num)]);

    for j = 1:U_num
        b_y(j,1) = max(unique(Us12cr.Set(j).V(:,2)));
    end
    for j = 1:U_num
        b_y(U_num+j,1) = -min(unique(Us12cr.Set(j).V(:,2)));
    end
%}
    
%{    
    % T3
    T_restrictions_no_cut = [zx_max(i:i+N-1) - T_bar_y*[x(i); xd(i); zx(i)] - constraints_restrict, ...
        zx_min(i:i+N-1) - T_bar_y*[x(i); xd(i); zx(i)] + constraints_restrict];
    
    for j = 1:length(constraints_restrict)
        if constraints_restrict(j) >= 0.03
            constraints_restrict(j) = 0.03;
        end
    end
    
%     disp([zx_max(i:i+N-1) - T_bar_y*[x(i); xd(i); zx(i)] ...
%         zx_max(i:i+N-1) - T_bar_y*[x(i); xd(i); zx(i)] - constraints_restrict ... 
%         constraints_restrict ...
%         -zx_min(i:i+N-1) + T_bar_y*[x(i); xd(i); zx(i)] ...
%         -zx_min(i:i+N-1) + T_bar_y*[x(i); xd(i); zx(i)] - constraints_restrict])
%     waitforbuttonpress

    % T
    T_no_restrictions = [zx_max(i:i+N-1) - T_bar_y*[x(i); xd(i); zx(i)], ...
        zx_min(i:i+N-1) - T_bar_y*[x(i); xd(i); zx(i)]];
    % T2
    T_restrictions_cut = [zx_max(i:i+N-1) - T_bar_y*[x(i); xd(i); zx(i)] - constraints_restrict, ...
        zx_min(i:i+N-1) - T_bar_y*[x(i); xd(i); zx(i)] + constraints_restrict];

%     plot(0:99,T(:,1), color, 0:99,T(:,2), color, ...
%         0:99,T2(:,1), color, 0:99,T2(:,2), color,...
%         0:99,T3(:,1), color, 0:99,T3(:,2), color)
%     disp([T T2])
%     waitforbuttonpress
%     plot(0:99,T(:,1), 0:99,T(:,2), ...
%         0:99,T2(:,1), 0:99,T2(:,2), ...
%         0:99,T3(:,1), 0:99,T3(:,2))

%     p1 = plot(0:5,T_no_restrictions(1:6,1), 'm', 'LineWidth',2);
%     hold on
%     p2 = plot(0:5,T_no_restrictions(1:6,2), 'm', 'LineWidth',2);
%     p3 = plot(0:5,T_restrictions_no_cut(1:6,1), ':r', 'LineWidth',2.5);
%     p4 = plot(0:5,T_restrictions_no_cut(1:6,2), ':r', 'LineWidth',2.5);
%     p5 = plot(0:5,T_restrictions_cut(1:6,1), '-.b', 'LineWidth',2);
%     p6 = plot(0:5,T_restrictions_cut(1:6,2), '-.b', 'LineWidth',2);
%     axis([-inf inf -0.05 0.05]);
%     legend([p1 p3 p5],{'No Restrictions','Restrictions/Not limited','Restrictions/limited'});
    
    p1 = plot(0:99,T_no_restrictions(:,1), 'm', 'LineWidth',2);
    hold on
    p2 = plot(0:99,T_no_restrictions(:,2), 'm', 'LineWidth',2);
    p3 = plot(0:99,T_restrictions_no_cut(:,1), ':r', 'LineWidth',2.5);
    p4 = plot(0:99,T_restrictions_no_cut(:,2), ':r', 'LineWidth',2.5);
    p5 = plot(0:99,T_restrictions_cut(:,1), '-.b', 'LineWidth',2);
    p6 = plot(0:99,T_restrictions_cut(:,2), '-.b', 'LineWidth',2);
    axis([-inf inf -0.05 0.05]);
    legend([p1 p3 p5],{'No Restrictions','Restrictions/Not limited','Restrictions/limited'});

    break;

    A = [S_bar_y; -S_bar_y];
    b_x = [ zx_max(i:i+N-1) - T_bar_y * [x(i); xd(i); zx(i)] - constraints_restrict;...
           -zx_min(i:i+N-1) + T_bar_y * [x(i); xd(i); zx(i)] - constraints_restrict];
    b_y = [ zy_max(i:i+N-1) - T_bar_y * [y(i); yd(i); zy(i)] - constraints_restrict;...
           -zy_min(i:i+N-1) + T_bar_y * [y(i); yd(i); zy(i)] - constraints_restrict];
%}

% %{
    on = 1;
    Zx = [zx_max(i:i+N-1) - (T_bar_y + S_bar_y * T_F_bar * on)*[x(i); xd(i); zx(i)], ...
          zx_min(i:i+N-1) - (T_bar_y + S_bar_y * T_F_bar * on)*[x(i); xd(i); zx(i)]]';
    Zy = [zy_max(i:i+N-1) - (T_bar_y + S_bar_y * T_F_bar * on)*[y(i); yd(i); zy(i)], ...
          zy_min(i:i+N-1) - (T_bar_y + S_bar_y * T_F_bar * on)*[y(i); yd(i); zy(i)]]';
    CR = [constraints_restrict -constraints_restrict]';
    if i==1
    for j = 1:size(Zx,2)
        cx(1) = max(Zx(:, j));
        cx(2) = min(Zx(:, j));
        cy(1) = max(Zy(:, j));
        cy(2) = min(Zy(:, j));
        v(1) = max(CR(:, j));
        v(2) = min(CR(:, j));
        CRR3(j) = Polyhedron('lb', [v(2);v(2);-inf], 'ub', [v(1);v(1);inf],'Ae',[0 0 1],'be',j);
        square(j) = Polyhedron('lb', [cx(2);cy(2);-inf], 'ub', [cx(1);cy(1);inf],'Ae',[0 0 1],'be',j);
        pont3(j) = square(j) - CRR3(j);
    end
   
    Us12 = PolyUnion(square);
    Us12cr = PolyUnion(pont3);
    figure;
    plot(Us12cr)
    title('ZMP Constraints restriction over the prediction horizon')
    ylabel('y[m]','FontSize',16')
    xlabel('x[m]','FontSize',16,'FontName','courier')
    zlabel('Time[iterations]','FontSize',16,'FontName','courier')
    axis([-inf inf -inf inf 0 100])
    break;
   end
%{
    U_num = Us12cr.Num;
    epsilon_over_2 = 0.0016;
    for j = 1:N
        if j <= U_num
            b_x(j,1) = max(unique(Us12cr.Set(j).V(:,1)));
        else
            Us12_j_unique = unique(Us12.Set(j).V(:,1));
            b_x(j,1) = (max(Us12_j_unique)+min(Us12_j_unique))/2 + epsilon_over_2;
        end
    end
    for j = 1:N
        if j <= U_num
            b_x(N+j,1) = -min(unique(Us12cr.Set(j).V(:,1)));
        else
            Us12_j_unique = unique(Us12.Set(j).V(:,1));
            b_x(N+j,1) = -(max(Us12_j_unique)+min(Us12_j_unique))/2 + epsilon_over_2;
        end
    end
%     disp([b_x(1:U_num) b_x(U_num+1:2*U_num)]);
    
    for j = 1:N
        if j <= U_num
            b_y(j,1) = max(unique(Us12cr.Set(j).V(:,2)));
        else
            Us12_j_unique = unique(Us12.Set(j).V(:,2));
            b_y(j,1) = (max(Us12_j_unique)+min(Us12_j_unique))/2 + epsilon_over_2;
        end
    end
    for j = 1:N
        if j <= U_num
            b_y(N+j,1) = -min(unique(Us12cr.Set(j).V(:,2)));
        else
            Us12_j_unique = unique(Us12.Set(j).V(:,2));
            b_y(N+j,1) = -(max(Us12_j_unique)+min(Us12_j_unique))/2 + epsilon_over_2;
        end
    end
    
%     S_bar_y_F_bar = S_bar_y*F_bar;
%     A = [S_bar_y_F_bar(1:U_num,:); -S_bar_y_F_bar(1:U_num,:)];
    A = [S_bar_y*F_bar; -S_bar_y*F_bar];
%}    
%}

% %{
    for j = 1:length(constraints_restrict)
        if constraints_restrict(j) >= 0.03
            constraints_restrict(j) = 0.03;
        end
    end

    A = [S_bar_y*F_bar; -S_bar_y*F_bar];
%     A = [S_bar_y; -S_bar_y];
    b_x = [ zx_max(i:i+N-1) - (T_bar_y + S_bar_y * T_F_bar)*[x(i); xd(i); zx(i)] - constraints_restrict;...
           -zx_min(i:i+N-1) + (T_bar_y + S_bar_y * T_F_bar)*[x(i); xd(i); zx(i)] - constraints_restrict];
    b_y = [ zy_max(i:i+N-1) - (T_bar_y + S_bar_y * T_F_bar)*[y(i); yd(i); zy(i)] - constraints_restrict;...
           -zy_min(i:i+N-1) + (T_bar_y + S_bar_y * T_F_bar)*[y(i); yd(i); zy(i)] - constraints_restrict];
%}

%{
    on = 1;
    A = [S_bar_y*F_bar; -S_bar_y*F_bar];
%     A = [S_bar_y; -S_bar_y];
    b_x = [ zx_max(i:i+N-1) - (T_bar_y + S_bar_y * T_F_bar * on)*[x(i); xd(i); zx(i); dfx(i)] - constraints_restrict;...
           -zx_min(i:i+N-1) + (T_bar_y + S_bar_y * T_F_bar * on)*[x(i); xd(i); zx(i); dfx(i)] - constraints_restrict];
    b_y = [ zy_max(i:i+N-1) - (T_bar_y + S_bar_y * T_F_bar * on)*[y(i); yd(i); zy(i); dfy(i)] - constraints_restrict;...
           -zy_min(i:i+N-1) + (T_bar_y + S_bar_y * T_F_bar * on)*[y(i); yd(i); zy(i); dfy(i)] - constraints_restrict];
       
    square = [zx_max(i:i+N-1)-(T_bar_y + S_bar_y * T_F_bar * on)*[x(i); xd(i); zx(i); dfx(i)] ...
              zx_min(i:i+N-1)-(T_bar_y + S_bar_y * T_F_bar * on)*[x(i); xd(i); zx(i); dfx(i)] ...
              zy_max(i:i+N-1)-(T_bar_y + S_bar_y * T_F_bar * on)*[y(i); yd(i); zy(i); dfy(i)] ...
              zy_min(i:i+N-1)-(T_bar_y + S_bar_y * T_F_bar * on)*[y(i); yd(i); zy(i); dfy(i)]];
          
%     square_y = [zy_max(i:i+N-1)-(T_bar_y + S_bar_y * T_F_bar * on)*[y(i); yd(i); zy(i)]*0 ...
%                 zy_min(i:i+N-1)-(T_bar_y + S_bar_y * T_F_bar * on)*[y(i); yd(i); zy(i)]*0 ...
%                 zy_max(i:i+N-1)-(T_bar_y + S_bar_y * T_F_bar * on)*[y(i); yd(i); zy(i)] ...
%                 zy_min(i:i+N-1)-(T_bar_y + S_bar_y * T_F_bar * on)*[y(i); yd(i); zy(i)]];
%}

%{
    for j = 1:length(constraints_restrict)
        if constraints_restrict(j) >= 0.03
            constraints_restrict(j) = 0.03;
        end
    end

    A = [S_bar_y; -S_bar_y];
    b_x = [ zx_max(i:i+N-1) - T_bar_y*[x(i); xd(i); zx(i)] - constraints_restrict;...
           -zx_min(i:i+N-1) + T_bar_y*[x(i); xd(i); zx(i)] - constraints_restrict];
    b_y = [ zy_max(i:i+N-1) - T_bar_y*[y(i); yd(i); zy(i)] - constraints_restrict;...
           -zy_min(i:i+N-1) + T_bar_y*[y(i); yd(i); zy(i)] - constraints_restrict];

    Aeq = (1-exp(-omega*delta))/omega * exp(-omega*delta*(0:N-1)) / (1-exp(-omega*delta*N));

    % x(i) is xc, xd(i) is xc_dot
    beq_x = x(i) + xd(i)/omega - zx(i);
    beq_y = y(i) + yd(i)/omega - zy(i);
    
    zd_x = quadprog(H,f_x,A,b_x,Aeq,beq_x);
    zd_y = quadprog(H,f_y,A,b_y,Aeq,beq_y);

    x_updated = A_upd*[x(i); xd(i); zx(i)] + B_upd*zd_x(1) + D*d;
    y_updated = A_upd*[y(i); yd(i); zy(i)] + B_upd*zd_y(1) + D*d;
    
    % eqn (7) at the paper
    z_pred_x = P*zd_x + zx(i);
    z_pred_y = P*zd_y + zy(i);
    
    x(i+1) = x_updated(1);
    xd(i+1) = x_updated(2);
    zx(i+1) = x_updated(3);
    
    y(i+1) = y_updated(1);
    yd(i+1) = y_updated(2);
    zy(i+1) = y_updated(3);
%}

%{
    u_lqi_x = -K_lqi*[x(i); xd(i); zx(i); x_int(i)];
    u_lqi_y = -K_lqi*[y(i); yd(i); zy(i); y_int(i)];
    
    x_int(i+1) = integr.a * x_int(i) + integr.b * (ref_x(i)-x(i));
    y_int(i+1) = integr.a * y_int(i) + integr.b * (ref_y(i)-y(i));

    zd_x = quadprog(H,f_x,A,b_x);
    zd_y = quadprog(H,f_y,A,b_y);

    x_updated = A_upd*[x(i); xd(i); zx(i)] + B_upd*(u_lqi_x*0+zd_x(1)) + D*d;
    y_updated = A_upd*[y(i); yd(i); zy(i)] + B_upd*(u_lqi_y*0+zd_y(1)) + D*d;
    
    % eqn (7) at the paper
    z_pred_x = P*zd_x + zx(i);
    z_pred_y = P*zd_y + zy(i);
    
    x(i+1) = x_updated(1);
    xd(i+1) = x_updated(2);
    zx(i+1) = x_updated(3);
    
    y(i+1) = y_updated(1);
    yd(i+1) = y_updated(2);
    zy(i+1) = y_updated(3);
%}

%{
    Aeq = (1-exp(-omega*delta))/omega * exp(-omega*delta*(0:N-1)) / (1-exp(-omega*delta*N));

    beq_x_F = (x(i) + xd(i)/omega - zx(i)) - Aeq*T_F_bar*[x(i); xd(i); zx(i)];
    beq_y_F = (y(i) + yd(i)/omega - zy(i)) - Aeq*T_F_bar*[y(i); yd(i); zy(i)];
    
    c_u_x = quadprog(H,f_x,A,b_x,Aeq*F_bar,beq_x_F);
    c_u_y = quadprog(H,f_y,A,b_y,Aeq*F_bar,beq_y_F);
    
    x_updated = phi*[x(i); xd(i); zx(i)] + B_upd*c_u_x(1);
    y_updated = phi*[y(i); yd(i); zy(i)] + B_upd*c_u_y(1);

    zd_x = F_bar*c_u_x + T_F_bar*[x(i); xd(i); zx(i)];
    zd_y = F_bar*c_u_y + T_F_bar*[y(i); yd(i); zy(i)];

    % eqn (7) at the paper
    z_pred_x = P*zd_x + zx(i);
    z_pred_y = P*zd_y + zy(i);
    
    x(i+1) = x_updated(1);
    xd(i+1) = x_updated(2);
    zx(i+1) = x_updated(3);
    
    y(i+1) = y_updated(1);
    yd(i+1) = y_updated(2);
    zy(i+1) = y_updated(3);
%}

% %{
    c_u_x = quadprog(H,f_x,A,b_x);
    c_u_y = quadprog(H,f_y,A,b_y);
    
    x_updated = phi*[x(i); xd(i); zx(i)] + B_upd*c_u_x(1) + 0*D*d;
    y_updated = phi*[y(i); yd(i); zy(i)] + B_upd*c_u_y(1) + 0*D*d;

    zd_x = F_bar*c_u_x + T_F_bar*[x(i); xd(i); zx(i)];
    zd_y = F_bar*c_u_y + T_F_bar*[y(i); yd(i); zy(i)];

    % eqn (7) at the paper
    z_pred_x = P*zd_x + zx(i);
    z_pred_y = P*zd_y + zy(i);
    
    x(i+1) = x_updated(1);
    xd(i+1) = x_updated(2);
    zx(i+1) = x_updated(3);
    
    y(i+1) = y_updated(1);
    yd(i+1) = y_updated(2);
    zy(i+1) = y_updated(3);
%}

%{
    c_u_x_d = quadprog(H,f_x,A,b_x);
    c_u_y_d = quadprog(H,f_y,A,b_y);

%     if abs(zy(i)) > 0
%         disp([square(1,:); [zx(i); zy(i); c_u_x(1); c_u_y(1)]']);
%     end
    
    x_updated_d = phi_d*[x(i); xd(i); zx(i); dfx(i)] + B_upd_d*c_u_x_d(1);
    y_updated_d = phi_d*[y(i); yd(i); zy(i); dfy(i)] + B_upd_d*c_u_y_d(1);
    
    zd_x = F_bar*c_u_x_d + T_F_bar*[x(i); xd(i); zx(i); dfx(i)];
    zd_y = F_bar*c_u_y_d + T_F_bar*[y(i); yd(i); zy(i); dfy(i)];

    % eqn (7) at the paper
    z_pred_x = P*zd_x + zx(i);
    z_pred_y = P*zd_y + zy(i);
    
    x(i+1) = x_updated_d(1);
    xd(i+1) = x_updated_d(2);
    zx(i+1) = x_updated_d(3);
    dfx(i+1) = x_updated_d(4);
    
    y(i+1) = y_updated_d(1);
    yd(i+1) = y_updated_d(2);
    zy(i+1) = y_updated_d(3);
    dfy(i+1) = y_updated_d(4);
%}

    clf
    hold on

    rect_x = [w,w,-w,-w,w]+shifting_x;
    rect_y = [foot_distance_y+w,-foot_distance_y-w,-foot_distance_y-w,foot_distance_y+w,foot_distance_y+w]+shifting_y;
    plot(rect_x,rect_y,'m','lineWidth',2);
    
    rect_x = [w,w,-w,-w,w]+shifting_x;
    rect_y = [w,-w,-w,w,w]+shifting_y;
    
    nPlottedFootsteps = 5;
    
    for j = 1:nPlottedFootsteps
        rect_x = [w,w,-w,-w,w]+shifting_x;
        rect_y = [w,-w,-w,w,w]+shifting_y;
        h1 = plot(fs_matrix(j,1)+rect_x,fs_matrix(j,2)+rect_y,'m','lineWidth',2);
    end
    
    h2 = plot(x,y,'r','lineWidth',2);
    h3 = plot(zx,zy,'b','lineWidth',2);
    h4 = plot(z_pred_x(1:end),z_pred_y(1:end)','g','lineWidth',2);
    axis equal
    axis([-0.2 1 -0.5 0.5]+shifting_x)
    
    xlabel('x [m]')
    ylabel('y [m]')
    
    grid on
    
    drawnow
%     break;
end