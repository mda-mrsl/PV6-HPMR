function v = vfaConstant(n, a, e1)
%VFACONSTANT creates variable flip angle scheme with constant Mxy
%   
%   Usage: v = vfaConstant(n, a, e1)
%
%       where n is the number of excitations
%             a is the effective excitation angle of all pulses in degrees
%               a = acosd(prod(cosd(v))) without T1 losses
%               default is 90 (full Mz usage)
%             e1 is the factor describing T1 loss/recovery between pulses
%               e1 = exp(-TR/T1) for HP magnetization
%               default is 1 (no T1 loss) - not currently used. (TODO)
%             v is a vector of length n with excitation angles in degrees
%
%   Literature: 
%     Nagashima, K. JMR 2008
%
%   See also IDDECOMP
%
%   06/2019, Keith Michel

%% Parse inputs
if nargin<1,    help(mfilename), return; end
if nargin<2,    a = []; end
if isempty(a),  a = 90; end
if nargin<3,    e1 = []; end
if isempty(e1), e1 = 1; end
a = min(max(a, 1.0), 90.0);
e1 = 1;

%% Calculate Mz and excitation angles forward in time with constant Mx
v = zeros(1, n);
X = sind(a)/sqrt(n); % Mxy value
x = zeros(1, n+1);
z = ones(1, n+1);
for ii = 1:n
    v(ii)   = asind(X/z(ii)/e1);
    x(ii+1) = z(ii)*sind(v(ii))*e1;
    z(ii+1) = z(ii)*cosd(v(ii))*e1;
end

return

%% Plot results
figure
subplot 121
plot(0:n,x,'o-', 0:n,z,'o-')
legend('M_X', 'M_Z')
title('Magnetization')
grid on
subplot 122
plot(1:n,v,'-o')
ylim([0, 90])
title('Excitation Angles')
grid on

% keyboard

end
