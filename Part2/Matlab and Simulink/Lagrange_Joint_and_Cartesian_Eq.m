clc
clear
close all

syms t1 t2 t3 t4 t5 t6 

d1=0.128; 
a2=0.6129;
a3=0.5716; 
d4=0.16389; 
d5=0.1157; 
d6=0.0922;

DOF=6;
a = [0 0 a2 a3 0 0];      % Link lengths
alpha = [0 pi/2 0 0 -pi/2 pi/2];        % Link twists
d = [d1 0 0 d4 d5 d6];          % Link offsets
theta = [t1 t2 t3 t4 t5 t6];      % Joint angles

% Mass properties
P_c1=[0; -0.01; -0.01];
P_c2=[0.25; 0; 0.17];
P_c3=[0.26; 0; 0.05];
P_c4=[0; 9.74; -7.6]*1e-3;
P_c5=[0; -9.74; -7.6]*1e-3;
P_c6=[0; -0.95; -17.46]*1e-3;
P_c=[P_c1 P_c2 P_c3 P_c4 P_c5 P_c6];

I1=0.03*eye(3);
I2=[0.05 0 0.01; 0 1.23 0; 0.01 0 1.22];
I3=[0.02 0 0; 0 0.54 0; 0 0 0.54];
I4=[3 0 0; 0 2.4 0.25; 0 0.25 2.8]*1e-3;
I5=[3 0 0; 0 2.4 -0.25; 0 -0.25 2.8]*1e-3;
I6=[2.2 0 0; 0 2.4 0; 0 0 4]*1e-4;
I=[I1 I2 I3 I4 I5 I6];

%% Transformation Matrix
T = eye(4);
z = [theta; theta; theta];
O_i = z;
O_ci = z;
R=sym('x',[3,18]); % Rotation Matrix [i to base]

for i=1:DOF

    T_i = [cos(theta(i))  -sin(theta(i))  0  a(i);
           sin(theta(i))*round(cos(alpha(i)))  cos(theta(i))*round(cos(alpha(i)))  round(-sin(alpha(i)))  round(-sin(alpha(i)))*d(i);
           sin(theta(i))*round(sin(alpha(i)))  cos(theta(i))*round(sin(alpha(i)))  round(cos(alpha(i)))    round(cos(alpha(i)))*d(i);
           0 0 0 1];
    
    T=T*T_i;
    T=simplify(T, 'IgnoreAnalyticConstraints', true);
    R(1:3,3*(i-1)+1:3*i)=T(1:3,1:3);

    z(:,i)=T(1:3,3);
    O_i(:,i)=T(1:3,4);
    O_ci(:,i)=T(1:3,4)+T(1:3,1:3)*P_c(:,i);

end

%% jacobian
jv_c=zeros(3)*sym('x',[3,36]); % Linear Jacobian Matrix [Center of i to base]
jw_c=zeros(3)*sym('x',[3,36]); % Rotational Jacobian Matrix [Center of i to base]

for q=1:DOF
    for i=1:q

        vec1=z(:,i);
        vec2=O_ci(:,q)-O_i(:,i);
        jv_c(1:3,6*(q-1)+i)=cross(vec1,vec2);
        jw_c(1:3,6*(q-1)+i)=z(:,i);
        jv_c=simplify(jv_c, 'IgnoreAnalyticConstraints', true);
        jw_c=simplify(jw_c, 'IgnoreAnalyticConstraints', true);

    end
end

%% Inertia Matrix
m=[8.3 23.52 12.56 1.967 1.967 0.43]; %[Kg]
M=zeros(6)*sym('x',6);

for i=1:DOF

    M = M + m(i)*transpose(jv_c(1:3,6*(i-1)+1:6*i))*jv_c(1:3,6*(i-1)+1:6*i) + transpose(jw_c(1:3,6*(i-1)+1:6*i))*(R(1:3,3*(i-1)+1:3*i)*I(1:3,3*(i-1)+1:3*i)*transpose(R(1:3,3*(i-1)+1:3*i)))*jw_c(1:3,6*(i-1)+1:6*i);
%     M=simplify(M, 'IgnoreAnalyticConstraints', true);
end

% M=subs(M,[t1 t2 t3 t4 t5 t6 d1 a2 a3 d4 d5 d6],[0 0 0 0 0 0 0.128 0.6129 0.5716 0.16389 0.1157 0.0922]);
% M=vpa(M,4);
% disp('M(0)')
% disp(M)

%% Gravity vector
z_ci=sym('x',[6,1]);
syms g 
Pot=0;
G=sym('x',[6,1]);
theta=[t1; t2; t3; t4; t5; t6]; % Angular Displacement

for i=1:DOF
    z_ci(i)=O_ci(3,i);
    Pot=Pot + m(i)*g*z_ci(i);
end

for i=1:DOF
    G(i)=diff(Pot, theta(i));
end

G=simplify(G, 'IgnoreAnalyticConstraints', true);
% G=subs(G,[t1 t2 t3 t4 t5 t6 d1 a2 a3 d4 d5 d6],[0 0 0 0 0 0 0.128 0.6129 0.5716 0.16389 0.1157 0.0922]);
% G=vpa(G,4);
% disp('G(0)')
% disp(G)

%% Coriolis and centrifugal vector
syms td1 td2 td3 td4 td5 td6
theta_dot=[td1; td2; td3; td4; td5; td6]; % Angular Velocities
C=zeros(6)*sym('x',6);

for i=1:DOF
    for n=1:DOF
        for m=1:DOF

        C(m,n) = C(m,n) + 1/2 * (diff(M(m,n),theta(i)) + diff(M(m,i),theta(n)) - diff(M(i,n),theta(m))) * theta_dot(i);
%         C(m,n) = simplify(C(m,n), 'IgnoreAnalyticConstraints', true);
        
        end
    end
end
V=C*theta_dot;
%%
% V=subs(V,[t1 t2 t3 t4 t5 t6 td1 td2 td3 td4 td5 td6 d1 a2 a3 d4 d5 d6],[0 0 0 0 0 0 1 0 0 0 0 0 0.128 0.6129 0.5716 0.16389 0.1157 0.0922]);
% V=vpa(V,4);
% disp('V')
% disp(V)

%% Cartesian State_Space Torque Equations
syms a2 a3 d1 d4 d5 d6 t1(time) t2(time) t3(time) t4(time) t5(time) t6(time)

DOF=6;
a = [0 0 a2 a3 0 0];      % Link lengths
alpha = [0 pi/2 0 0 -pi/2 pi/2];        % Link twists
d = [d1 0 0 d4 d5 d6];          % Link offsets
theta = [t1(time) t2(time) t3(time) t4(time) t5(time) t6(time)];      % Joint angles by time

T = eye(4);
z = [theta; theta; theta];
O_i= z;

for i=1:DOF
    T_i = [cos(theta(i))  -sin(theta(i))  0  a(i);
           sin(theta(i))*round(cos(alpha(i)))  cos(theta(i))*round(cos(alpha(i)))  round(-sin(alpha(i)))  round(-sin(alpha(i)))*d(i);
           sin(theta(i))*round(sin(alpha(i)))  cos(theta(i))*round(sin(alpha(i)))  round(cos(alpha(i)))    round(cos(alpha(i)))*d(i);
           0 0 0 1];

    T=T*T_i;
    T=simplify(T, 'IgnoreAnalyticConstraints', true);
    z(:,i)=T(1:3,3);
    O_i(:,i)=T(1:3,4);
end

%%
O_e=T(1:3,4);
j=sym('x',DOF);
for i=1:DOF
        vec1=z(:,i);
        vec2=O_e-O_i(:,i);
        j(1:3,i)=cross(vec1,vec2);
        j(4:6,i)=z(:,i);
        j=simplify(j, 'IgnoreAnalyticConstraints', true);
%         disp(['J', num2str(i)])
%         disp(j(1:6,i))
end

%%
M_x=M*j^(-1); % [Mass Matrix in Cartesian Space]

% tdi is the Angular Velocity of i_th joint
j_dot=diff(j, time);
j_dot=subs(j_dot, [diff(t1(time),time) diff(t2(time),time) diff(t3(time),time) diff(t4(time),time) diff(t5(time),time) diff(t6(time),time)], [td1 td2 td3 td4 td5 td6]);

V_x=V-M_x*j_dot*theta_dot;
G_x=G;
