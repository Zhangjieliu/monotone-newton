%% The admittance matrix Y, YSS, YSL, YLS, YLL
Y=zeros(123,123);
Y(1,23)=1; Y(23,24)=1;Y(21,22)=0.5; Y(22,24)=0.5; Y(24,25)=1; Y(25,26)=1;
Y(26,29)=1; Y(28,29)=1; Y(27,28)=1; Y(11,27)=0.2; Y(29,34)=0.5; Y(2,34)=0.25;
Y(34,35)=0.25; Y(35,44)=0.1; Y(43,44)=0.2; Y(16,43)=0.2; Y(44,57)=0.5; Y(35,36)=0.2;
Y(44,45)=0.2; Y(45,65)=0.5; Y(64,65)=0.2; Y(65,66)=1; Y(15,66)=0.25; Y(42,45)=0.5;
Y(40,42)=0.2; Y(40,41)=0.5; Y(37,40)=0.5; Y(37,38)=0.25; Y(38,39)=0.2; Y(31,37)=0.5;
Y(30,31)=1; Y(31,32)=0.2; Y(32,33)=0.5; Y(12,33)=0.2; Y(14,42)=0.2; Y(57,58)=0.5; Y(58,59)=0.5;
Y(4,59)=0.25; Y(59,63)=0.5; Y(55,57)=0.2; Y(55,56)=0.5; Y(51,55)=1; Y(48,51)=1; Y(47,48)=1;
Y(46,47)=1; Y(3,46)=0.5; Y(46,49)=1; Y(54,55)=1; Y(53,54)=1; Y(50,53)=1; Y(52,53)=1;
Y(53,60)=0.2; Y(5,60)=1; Y(60,61)=0.5; Y(61,62)=0.5; Y(57,74)=0.1; Y(74,75)=0.2; Y(75,76)=0.2;
Y(76,77)=0.5; Y(19,77)=0.5; Y(68,76)=0.1; Y(67,68)=0.25; Y(17,67)=0.25; Y(68,69)=0.2;
Y(69,89)=0.5; Y(89,90)=0.2; Y(69,70)=0.5; Y(69,92)=0.1; Y(92,93)=0.5; Y(93,94)=0.5;
Y(94,95)=0.5; Y(13,95)=0.25; Y(71,72)=0.5; Y(70,71)=0.5; Y(72,73)=0.5; Y(72,108)=0.25;
Y(73,112)=0.5; Y(112,113)=0.5; Y(8,113)=0.2; Y(108,109)=0.5; Y(110,111)=0.5; Y(105,108)=0.5;
Y(105,106)=0.5;Y(106,107)=0.5;Y(7,107)=0.5;Y(102,105)=0.5;Y(102,103)=0.5;Y(103,104)=0.5;
Y(101,102)=0.5;Y(100,101)=0.5;Y(96,100)=0.5;Y(97,100)=0.5;Y(96,97)=0.5;Y(97,98)=0.5;Y(98,99)=0.5;
Y(6,99)=0.2;Y(109,110)=0.25;Y(72,91)=0.25;Y(91,114)=0.5;Y(114,115)=0.5;Y(115,116)=0.5;
Y(9,116)=0.2;Y(88,91)=0.1;Y(88,119)=1;Y(87,88)=0.1;Y(86,87)=0.1;Y(85,86)=0.2;Y(84,86)=0.25;
Y(83,84)=0.2;Y(82,84)=0.25;Y(81,82)=0.2;Y(80,82)=0.25;Y(79,80)=0.2;Y(78,80)=0.25;
Y(20,78)=0.2;Y(118,119)=1;Y(18,118)=0.5;Y(118,120)=0.5;Y(120,122)=0.5;Y(121,122)=0.2;
Y(122,123)=0.2;Y(10,123)=0.2;Y(117,121)=0.2;
Y=Y+Y'; Y=diag(Y*ones(123,1))-Y;

YSS=Y(1:10, 1:10);YSL=Y(1:10, 11:123);YLS=Y(11:123, 1:10);YLL=Y(11:123,11:123);

%% voltage reverence, open-circuit voltage and load vecter 
v_S=ones(10,1); % voltage reverence
p=[-0.1*ones(10,1);0.006*ones(103,1)]; % load
v=-YLL^(-1)*YLS*v_S; %open-circuit voltage

%% The NR method based on the fractional formulation in (4)
 A=diag(v)^(-1)*YLL^(-1)*diag(v)^(-1)*diag(p); % A in equation (4)
 x=v*ones(1, 21); % initinalization 
 for k=1:20
     fx=eye(113)-A*diag(x(:,k))^(-2);   % Jacobian matrix f_x(x_k)
     f=x(:,k)+A*(1./x(:,k))-ones(113,1);  % residual  f(x_k)
     x(:,k+1)=x(:,k)-fx^(-1)*f;  % NR update
 end

 %% The NR method based on the quadratic formulation in (3)
 xx=v*ones(1, 10); % initinalization 
for k=1:10
     gx=diag(xx(:,k))*YLL+diag(YLL*xx(:,k)+YLS*v_S);  % Jacobian matrix g_x(x_k)
     g=diag(xx(:,k))*(YLL*xx(:,k)+YLS*v_S)+p;   % residual g(x_k)
     xx(:,k+1)=xx(:,k)-gx^(-1)*g;   % NR update
end

%% Monotone convergence whenever a solution exists
alpha=0.1:0.001:2;     % alpha starts at 0.1 and increases in steps of 0.001
vh=[];        % high-voltage equilibrium at each power point
for j=1:length(alpha)   % sweep over the power points
    pa=[-0.1*ones(10,1);0.006*alpha(j)*ones(103,1)];   % update the loads
    u=v;       % open-circuit initialization
    for k=1:20
        g=diag(u)*(YLL*u+YLS*v_S)+pa;  %  residual  g(x_k)
        gx=diag(u)*YLL+diag(YLL*u+YLS*v_S);   % Jacobian matrix g_x(x_k)
        d=gx\g;    % Newton correction
        if k>1 && min(d)<0, break; end;  % monotonicity check after the first update
        u=u-d; % NR update
    end
    if ~(norm(g,inf)<1e-10), break; end;   % convergence check
    vh(:,j)=u;   % store the high-voltage equilibrium
end 
alpha(j-1)     % display the estimated critical loading factor

%%  Plotting Fig. 1
figure(1); clf;
plot([0 NaN 1:20],[x(:,1) nan(113,1) x(:,2:end)]','--*');
xlabel('k'); ylabel('Voltage (p.u.)'); title('(a)'); ax(1)=gca;

figure(2); clf;
plot([0 NaN 1:size(xx,2)-1],[xx(:,1) nan(113,1) xx(:,2:end)]','--*');
xlabel('k'); ylabel('Voltage (p.u.)'); title('(b)'); ax(2)=gca;

figure(3); clf;
[~,n]=min(vh(:,end));
plot(alpha(1:size(vh,2)),vh(n,:));
xlabel('\alpha'); ylabel('Voltage (p.u.)'); title('(c)'); ax(3)=gca;


