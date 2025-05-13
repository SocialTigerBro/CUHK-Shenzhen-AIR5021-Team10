%2020.10.16，动力学建模，不再更改
clc
clear
syms   q1 q2 q3 q4 q5 q6 Q1 Q2 Q3 Q4 Q5 Q6 real %关节变量即广义变量，既可以指角度也可以指长度
syms   l1 l2 l3 l4 l5 l6 lc1 lc2 lc3 lc4 lc5 lc6 lcm real %连杆长度和连杆质心位置
syms   eq1 eq2 eq3 eq4 eq5 eq6 eeq1 eeq2 eeq3 eeq4 eeq5 eeq6 real %广义变量一二阶导数
syms   m m1 m2 m3 m4 m5 m6 d1 d4 d6 a2 a3 real
syms   Ic1x Ic1y Ic1z Ic2x Ic2y Ic2z Ic3x Ic3y Ic3z Ic4x Ic4y Ic4z Ic6x Ic6y Ic6z Icmx Icmy Icmz real %各连杆质量和转动惯量
syms   tq1 tq2 tq3 tq4 tq5 tq6 g s real %时间变量，重力加速度，占位符号
M = diag([s s s s s s]);%创建容器，准备赋值
D = diag([s s s s s s]);%创建容器，准备赋值
K = diag([s s s s s s]);%创建容器，准备赋值
G = [s s s s s s]';%创建容器，准备赋值

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%数值赋值%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%屏蔽这部分则产生的是不带字符的形式（方便计算），不屏蔽则产生完整的形式（方便写论文）
l1 = 0.05;  l2 = 0.36;  l3 = 0.23;  l4 = 0.08;  l5 = 0.05;  l6 = 0.05;  lm = 0.05;
lc1 = l1/2; lc2 = l2/2; lc3 = l3/2; lc4 = l4/2; lc5 = l5/2; lc6 = l6/2; lcm = lm/2;
m1 = 0.67;
m2 = 1.19;
m3 = 0.84;
m4 = 0.56;
m5 = 0.39;
m6 = 0.29;
m = 0.8;%执行器
g = 9.8;


%%根据自己建立的固连坐标系调整转动惯量对应的位置
Ic1x = 0.00128328; Ic1y = 0.00083936; Ic1z = 0.00071931;
Ic2x = 0.00102138; Ic2y = 0.02466114; Ic2z = 0.02429457;
Ic3x = 0.00108061; Ic3y = 0.00886621; Ic3z = 0.00954238;
Ic4x = 0.00031576; Ic4y = 0.00092996; Ic4z = 0.00097912;
Ic5x = 0.00055896; Ic5y = 0.00053860; Ic5z = 0.00017605;
Ic6x = 0.00014750; Ic6y = 0.00014680; Ic6z = 0.00018328;
Ic1 =[Ic1x    0       0;
      0       Ic1y    0;
      0       0       Ic1z];
Ic2 =[Ic2x    0       0;
      0       Ic2y    0;
      0       0       Ic2z];
Ic3 =[Ic3x    0       0;
      0       Ic3y    0;
      0       0       Ic3z];
Ic4 =[Ic4x    0       0;
      0       Ic4y    0;
      0       0       Ic4z];
Ic6 =[Ic6x    0       0;
      0       Ic6y    0;
      0       0       Ic6z];
Icm =[Icmx    0       0;
      0       Icmy    0;
      0       0       Icmz];

%%%%%%%%%%%%%%%%%%%%%%定义广义变量及其一、二阶导数的关系%%%%%%%%%%%%%%%%%%%
q1 = 0.5*eeq1*tq1^2;
q2 = 0.5*eeq2*tq2^2;
q3 = 0.5*eeq3*tq3^2;
q4 = 0.5*eeq4*tq4^2;
q5 = 0.5*eeq5*tq5^2;
q6 = 0.5*eeq6*tq6^2;
%给公式的时候一定要给出t的后缀，后面做手动矫正时有用

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%根据D-H参数表格列写齐次变换矩阵%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T01 = [cos(q1) -sin(q1) 0 0
       sin(q1) cos(q1) 0 0
       0 0 1 0.1
       0 0 0 1];

T12 = [cos(q2) -sin(q2) 0 0
       0 0 -1 0
       sin(q2) cos(q2) 0 0
       0 0 0 1];

T23 = [cos(q3+0.2191) -sin(q3+0.2191) 0 0.36
       -sin(q3+0.2191) -cos(q3+0.2191) 0 0
       0 0 -1 0
       0 0 0 1];

T34 = [cos(q4-0.2191) -sin(q4-0.2191) 0 -0.23
       sin(q4-0.2191) cos(q4-0.2191) 0 0
       0 0 1 0
       0 0 0 1];

T45 = [cos(q5+pi/2) -sin(q5+pi/2) 0 -0.08
       0 0 -1 0
       sin(q5+pi/2) cos(q5+pi/2) 0 0
       0 0 0 1];

T56 = [cos(q6) -sin(q6) 0 0
       0 0 1 0.05
       -sin(q6) -cos(q6) 0 0
       0 0 0 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算各部分质心坐标%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%这里所求的质心坐标全部是为了后面去求质心的线速度
C11 = [0 0 -lc1 1]';
C1E = simplify(T01*C11);

C22 = [lc2 0 0 1]';
C2E = simplify(T01*T12*C22);
%2号连杆可以同时绕着Z0和Z1轴旋转，所以把2号连杆放在2号坐标系中，就可以同时包含两个旋转运动

C33 = [-lc3 0 0 1]';%杆l3质心
C3E = simplify(T01*T12*T23*C33);

C44 = [-lc4 0 0 1]';%杆l4质心
C4E = simplify(T01*T12*T23*T34*C44);%这里4号连杆本身是坐标系的Z轴，所以表示方法发生了变化

C55 = [0 lc5 0 1]';
C5E = simplify(T01*T12*T23*T34*T45*C55);

C66 = [0 0 lc6 1]'; 
C6E = simplify(T01*T12*T23*T34*T45*T56*C66);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算各个连杆质心运动时的线速度%%%%%%%%%%%%%%%%%%%%%%%%%
c = [C1E(1:3),C2E(1:3),C3E(1:3),C4E(1:3),C5E(1:3),C6E(1:3)];%封装，便于循环
%这里没有对C1E求导，因为后面求动能时没有用到1号连杆质心的线速度
vc_t = sym(zeros(3,6));
tq = [tq1 tq2 tq3 tq4 tq5 tq6];
for i = 1:6
    for j = 1:6
        vc_t(:,i) = vc_t(:,i) + diff(c(:,i),tq(j));
    end
end

%对t求导，相当于对所有具有不同 t 后缀的 t分量 求导，再相加


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算各部分质心角速度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T65 = T56(1:3,1:3)';%先求一下变换矩阵的反变换，即转置
T54 = T45(1:3,1:3)';
T43 = T34(1:3,1:3)';
T32 = T23(1:3,1:3)';
T21 = T12(1:3,1:3)';
T10 = T01(1:3,1:3)';
%注意，角度和位移进行统一坐标系的时候，变换的方向不同

W11 = [0 0 eq1]';%连杆l1绕Z1轴转
W22 = T21*W11+[0 0 eq2]';
W33 = T32*W22+[0 0 eq3]';
W44 = T43*W33+[0 0 eq4]';
W55 = T54*W44+[0 0 eq5]';    
W66 = T65*W55+[0 0 eq6]';


%将角速度转化到同一个坐标系就可以进行正常的加减运算
%转化时需要乘以齐次变换矩阵左上角三维矩阵，它即为只涉及旋转不涉及位移的三维变换矩阵
%角速度的方向一定和自身坐标系的Z轴平行，而拉格朗日方程中的角速度是自身坐标系中的角速度（Z轴）
%所以需要将前面所有的角速度转化到目标连杆坐标系中再进行求和

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%中间变量替换%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vc = sym(zeros(3,6));%创建容器
wc_t = [W11,W22,W33,W44,W55,W66];%封装，便于循环替换,相比最初的版本增加了1号连杆的自转
wc = sym(zeros(3,6));
for i = 1:6
    vc(:,i) = subs(vc_t(:,i),{q1,q2,q3,q4,q5,q6,eeq1*tq1,eeq2*tq2,eeq3*tq3,eeq4*tq4,eeq5*tq5,eeq6*tq6},{Q1,Q2,Q3,Q4,Q5,Q6,eq1,eq2,eq3,eq4,eq5,eq6});
end
for i = 1:6
    wc(:,i) = subs(wc_t(:,i),{q1,q2,q3,q4,q5,q6,eeq1*tq1,eeq2*tq2,eeq3*tq3,eeq4*tq4,eeq5*tq5,eeq6*tq6},{Q1,Q2,Q3,Q4,Q5,Q6,eq1,eq2,eq3,eq4,eq5,eq6}); 
end

%对中间变量求导时要避免接触底层变量,所以要进行变量替换，把式子中所有的 t 消掉
%这种替换只能替换掉完全一样的部分，matlab不能智能的运用乘法交换律，无法完美的替换掉
%但是可以通过上式的结果进行手动矫正，因为t的后缀都给出了，所以矫正的过程很容易

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算各个连杆的动能%%%%%%%%%%%%%%%%%%%%%%%%%%%
E1 = simplify(0.5*m1*(vc(:,1)'*vc(:,1))+0.5*wc(:,1)'*Ic1*wc(:,1));
E2 = simplify(0.5*m2*(vc(:,2)'*vc(:,2))+0.5*wc(:,2)'*Ic2*wc(:,2));
E3 = simplify(0.5*m3*(vc(:,3)'*vc(:,3))+0.5*wc(:,3)'*Ic3*wc(:,3));
E4 = simplify(0.5*m4*(vc(:,4)'*vc(:,4))+0.5*wc(:,4)'*Ic4*wc(:,4));
E5 = simplify(0.5*m4*(vc(:,5)'*vc(:,5))+0.5*wc(:,5)'*Ic4*wc(:,5));
E6 = simplify(0.5*m5*(vc(:,6)'*vc(:,6))+0.5*wc(:,6)'*Ic6*wc(:,6));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算各个连杆的势能%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U1 = subs(m1*g*C1E(3),{q1,q2,q3,q4,q5,q6,eeq1*tq1,eeq2*tq2,eeq3*tq3,eeq4*tq4,eeq5*tq5,eeq6*tq6},{Q1,Q2,Q3,Q4,Q5,Q6,eq1,eq2,eq3,eq4,eq5,eq6});
U2 = subs(m2*g*C2E(3),{q1,q2,q3,q4,q5,q6,eeq1*tq1,eeq2*tq2,eeq3*tq3,eeq4*tq4,eeq5*tq5,eeq6*tq6},{Q1,Q2,Q3,Q4,Q5,Q6,eq1,eq2,eq3,eq4,eq5,eq6});
U3 = subs(m3*g*C3E(3),{q1,q2,q3,q4,q5,q6,eeq1*tq1,eeq2*tq2,eeq3*tq3,eeq4*tq4,eeq5*tq5,eeq6*tq6},{Q1,Q2,Q3,Q4,Q5,Q6,eq1,eq2,eq3,eq4,eq5,eq6});
U4 = subs(m4*g*C4E(3),{q1,q2,q3,q4,q5,q6,eeq1*tq1,eeq2*tq2,eeq3*tq3,eeq4*tq4,eeq5*tq5,eeq6*tq6},{Q1,Q2,Q3,Q4,Q5,Q6,eq1,eq2,eq3,eq4,eq5,eq6});
U5 = subs(m5*g*C5E(3),{q1,q2,q3,q4,q5,q6,eeq1*tq1,eeq2*tq2,eeq3*tq3,eeq4*tq4,eeq5*tq5,eeq6*tq6},{Q1,Q2,Q3,Q4,Q5,Q6,eq1,eq2,eq3,eq4,eq5,eq6});
U6 = subs(m6*g*C6E(3),{q1,q2,q3,q4,q5,q6,eeq1*tq1,eeq2*tq2,eeq3*tq3,eeq4*tq4,eeq5*tq5,eeq6*tq6},{Q1,Q2,Q3,Q4,Q5,Q6,eq1,eq2,eq3,eq4,eq5,eq6});



%质心到各个方向质量相同，重力势能又是一个线性变化，直接用质心的高度就可以了。
%因为我们先求的是L对eq的偏导数，不涉及 t，所以要直接使用 Q ，而不能用 q,因此需要矫正


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%拉格朗日函数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = E1 + E2 + E3 + E4 + E5 + E6  - U1 - U2 - U3 - U4 - U5 - U6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算拉格朗日方程中各项%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%L_eq1_t:
eq = [eq1 eq2 eq3 eq4 eq5 eq6];
L_eq_t = sym(zeros(1,6));
for i = 1:6
    temp = diff(L,eq(i));
    temp1 = simplify(subs(temp,{Q1,Q2,Q3,Q4,Q5,Q6,eq1,eq2,eq3,eq4,eq5,eq6},{q1,q2,q3,q4,q5,q6,eeq1*tq1,eeq2*tq2,eeq3*tq3,eeq4*tq4,eeq5*tq5,eeq6*tq6}));%对t求导之前，替换回所有变量
    temp2 = 0;
    for j = 1:6
        temp2 = temp2 + (diff(temp1,tq(j)));
    end
    L_eq_t(i) = subs(temp2,{q1,q2,q3,q4,q5,q6,eeq1*tq1,eeq2*tq2,eeq3*tq3,eeq4*tq4,eeq5*tq5,eeq6*tq6},{Q1,Q2,Q3,Q4,Q5,Q6,eq1,eq2,eq3,eq4,eq5,eq6});
end

L_eq_t = subs(L_eq_t,{eq1*tq1,eq2*tq2,eq3*tq3,eq4*tq4,eq5*tq5,eq6*tq6},{2*Q1,2*Q2,2*Q3,2*Q4,2*Q5,2*Q6});%不应该再出现tq的，需要手动矫正

L_Q1 = diff(L,Q1);
L_Q2 = diff(L,Q2);
L_Q3 = diff(L,Q3);
L_Q4 = diff(L,Q4);
L_Q5 = diff(L,Q5);
L_Q6 = diff(L,Q6);

%%%%%%%%%%%%%%%%%%%%%%代入拉格朗日方程中建立系统模型%%%%%%%%%%%%%%%%%%%
t1 = simplify(L_eq_t(1) - L_Q1);
t2 = simplify(L_eq_t(2) - L_Q2);
t3 = simplify(L_eq_t(3) - L_Q3);
t4 = simplify(L_eq_t(4) - L_Q4);
t5 = simplify(L_eq_t(5) - L_Q5);
t6 = simplify(L_eq_t(6) - L_Q6);

%这里的simplify函数非常非常非常非常非常重要，一定在这里加，并且不可遗漏。否则不容易验证结果！！！
%整个程序有三个地方需要手动矫正，一个是"中间变量替换"，一个是“势能计算”，一个是“勒格朗日各项”
%subs()函数表示将符号表达式中的某些符号变量替换为指定的新的变量
%为了对中间变量求导，就必须把底层变量替换掉，而在替换的过程中，一定要避免再使用第一部分的式子,所以要用Q
%对t求导时一定要把所有中间变量再替换为底层变量

%%%%%%%%%%%%%%%%%%%%%%%%按照机械臂动力学模型标准形式提取系数矩阵%%%%%%%%%%%%%%%%%%%%
t =[t1 t2 t3 t4 t5 t6];
q = {'Q1','Q2','Q3','Q4','Q5','Q6'};
eq = {'eq1','eq2','eq3','eq4','eq5','eq6'};
eeq = {'eeq1','eeq2','eeq3','eeq4','eeq5','eeq6'};
global rest
for i = 1:6
    if t(i) == 0
        rest =s;
    else
    rest = t(i);
    end
    for j = 1:6
        [temp] = poly_coeffs(rest,eeq{j});
        rest = temp(1,1);
        if length(temp) == 1
            M(i,j) = 0;
        else if length(temp) == 2
                M(i,j) = temp(1,2);
            else if length(temp) == 3
                M(i,j) = temp(1,2)+temp(1,3)*eeq{j};
                end
            end
        end
    end
    
    
    for k = 1:6
        [temp] = poly_coeffs(rest,eq{k});
        rest = temp(1,1);
        if length(temp) == 1
            D(i,k) = 0;
        else if length(temp) == 2
                D(i,k) = temp(1,2);
            else if length(temp) == 3
                D(i,k) = temp(1,2)+temp(1,3)*eq{k};
                end
            end
        end
        if rest == 0
            rest = s;
        end
    end
    
    
    
    if rest == 0
        K(i,1:6) = [0 0 0 0 0 0];
    else
        for m = 1:6
        [temp] = poly_coeffs(rest,q{m});
        rest = temp(1,1);
        if length(temp) == 1
            K(i,m) = 0;
        else if length(temp) == 2
                K(i,m) = temp(1,2);
            else if length(temp) == 3
                K(i,m) = temp(1,2)+temp(1,3)*q{m};
                end
            end
        end
        end
    end
          
    if rest == 's'
        G(i) = 0;
    else
        G(i) = rest;
    end
    
end
