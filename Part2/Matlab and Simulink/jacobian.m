clc
clear
close all

syms a2 a3 d1 d4 d5 d6 t1 t2 t3 t4 t5 t6 

DOF=6;
a = [0 0 a2 a3 0 0];      % Link lengths
alpha = [0 pi/2 0 0 -pi/2 pi/2];        % Link twists
d = [d1 0 0 d4 d5 d6];          % Link offsets
theta = [t1 t2 t3 t4 t5 t6];      % Joint angles

T = eye(4);
z = [theta; theta; theta];
O_i= z;

for i=1:DOF
    T_i = [cos(theta(i))  -sin(theta(i))  0  a(i);
           sin(theta(i))*round(cos(alpha(i)))  cos(theta(i))*round(cos(alpha(i)))  round(-sin(alpha(i)))  round(-sin(alpha(i)))*d(i);
           sin(theta(i))*round(sin(alpha(i)))  cos(theta(i))*round(sin(alpha(i)))  round(cos(alpha(i)))    round(cos(alpha(i)))*d(i);
           0 0 0 1];
    disp(['T', num2str(i)])
    T=T*T_i;
    T=simplify(T, 'IgnoreAnalyticConstraints', true);
    disp(T)
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
        disp(['J', num2str(i)])
        disp(j(1:6,i))
end

dtr_j=det(j);
simplify(dtr_j, 'IgnoreAnalyticConstraints', true);



