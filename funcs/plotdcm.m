function plotdcm(DCM)

y  = DCM.xY.y;
p  = DCM.H;

% y = f(x) + e
y = spm_cat(y'); % y
p = spm_cat(p'); % f(x)
e = y - p;       % e

plot(y,':','LineWidth',2);hold on; plot(p,'LineWidth',2);

% error
offset = 1.3*min([y(:); p(:)]);
plot(e+offset,'g:','LineWidth',2);
grid on;

channel = DCM.xY.name;
trial   = DCM.xY.code;

legend(channel);
set(gca,'fontsize',18);