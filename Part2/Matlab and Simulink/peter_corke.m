% clc
% clear
% close all

DOF=6;
a = [0 0 612.9 571.6 0 0]*1e-3;      % Link lengths
alpha = [0 pi/2 0 0 -pi/2 pi/2];        % Link twists
d = [128 0 0 163.89 115.7 92.2]*1e-3;         % Link offsets
theta = [input('enter theta1 (rad): ') input('enter theta2 (rad): ') input('enter theta3 (rad): ') input('enter theta4 (rad): ') input('enter theta5 (rad): ') input('enter theta6 (rad): ')];      % Joint angles
T = eye(4);

for i=1:DOF

    T_i = [cos(theta(i))  -sin(theta(i))  0  a(i);
           sin(theta(i))*round(cos(alpha(i)))  cos(theta(i))*round(cos(alpha(i)))  round(-sin(alpha(i)))  round(-sin(alpha(i)))*d(i);
           sin(theta(i))*round(sin(alpha(i)))  cos(theta(i))*round(sin(alpha(i)))  round(cos(alpha(i)))    round(cos(alpha(i)))*d(i);
           0 0 0 1];

    T=T*T_i;

end
disp('T_e')
disp(T)

startup_rvc 

t1 = theta(1);
t2 = theta(2);
t3 = theta(3);
t4 = theta(4);
t5 = theta(5);
t6 = theta(6);

L(1)= Link( [t1      d(1)     a(1)     alpha(1)]);
L(2)= Link( [t2      d(2)     a(2)     alpha(2)]);
L(3)= Link( [t3      d(3)     a(3)     alpha(3)]);
L(4)= Link( [t4      d(4)     a(4)     alpha(4)]);
L(5)= Link( [t5      d(5)     a(5)     alpha(5)]);
L(6)= Link( [t6      d(6)     a(6)     alpha(6)]);

Robot=SerialLink(L);
Robot.name='UR10';
Robot.fkine([t1 t2  t3 t4 t5 t6])
figure(1)
Robot.plot([t1 t2  t3 t4 t5 t6])

syms t1 t2 t3 t4 t5 t6
Robot.fkine([t1 t2  t3 t4 t5 t6])


%  for t = 0:0.1:pi/2
%      Rob.plot([t 0 0 0 0 0]);
%      pause(.25)
%  end
