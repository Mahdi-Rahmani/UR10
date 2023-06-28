clc
clear
close all

syms a2 a3 d1 d4 d5 d6 t1 t2 t3 t4 t5 t6 
syms r11 r12 r13 r21 r22 r23 r31 r32 r33 x y z

Td=[r11 r12 r13 x;r21 r22 r23 y;r31 r32 r33 z;0 0 0 1];
DOF=6;
a = [0 0 612.9 571.6 0 0]*1e-3;      % Link lengths
alpha = [0 pi/2 0 0 -pi/2 pi/2];        % Link twists
d = [128 0 0 163.89 115.7 92.2]*1e-3;         % Link offsets
theta = [t1 t2 t3 t4 t5 t6];      % Joint angles
T = eye(4);

for i=1:DOF
    T_i = [cos(theta(i))  -sin(theta(i))  0  a(i);
           sin(theta(i))*round(cos(alpha(i)))  cos(theta(i))*round(cos(alpha(i)))  round(-sin(alpha(i)))  round(-sin(alpha(i)))*d(i);
           sin(theta(i))*round(sin(alpha(i)))  cos(theta(i))*round(sin(alpha(i)))  round(cos(alpha(i)))    round(cos(alpha(i)))*d(i);
           0 0 0 1];

    disp(['T', num2str(i)])
    T=T*T_i;
    T=simplify(T, 'IgnoreAnalyticConstraints', true);
    disp(T)
end

T=subs(T,[t1,t2,t3,t4,t5,t6],[0,0,0,0,0,0]);
T=vpa(T,4);
disp('T6')
disp(T)

