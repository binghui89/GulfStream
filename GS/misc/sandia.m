clear; clc;
cap = [4, 40, 200, 400];
c_11 = [10783335, 22271891, 42697218, 48954903]; % 1.1, development y = a + b*log(x)
c_12 = [34884000, 36125350, 74847360, 119694720];
c_13 = [2220884, 18353305,89913512, 179454125]; % 1.3, mooring, y = a + bx
c_14 = [8856020, 58581409, 249475395, 475762842]; % 1.4, device structural components, y = a + bx
c_15 = [21622845, 150396585, 600781220, 1098683153]; % 1.5, power take off, y = a + bx
c_16 = [3047886, 20897799, 85025662, 157444600]; % 1.6, Subsystem integration & profit margin, y = a + bx
c_17 = [12805055, 33025253, 95817735, 186450876]; % 1.7, Installation, y = a + bx

c_21 = [1668734 , 6347594 , 11958609 , 11087452];
c_22 = [1478750 , 1972900 , 2036150 , 2036150];
c_23 = [114995 , 1149954 , 5749770 , 11499540];
c_24 = [312129 , 373182 , 1003520 , 1773726];
c_25 = [731280 , 4582328 , 20577411 , 39994455];
c_26 = [17494 , 174936 , 874680 , 1749361];

p_11 = polyfit(log10(cap), c_11, 1);
p_12 = polyfit(cap, c_12, 1);
p_13 = polyfit(cap, c_13, 1);
p_14 = polyfit(cap, c_14, 1);
p_15 = polyfit(cap, c_15, 1);
p_16 = polyfit(cap, c_16, 1);
p_17 = polyfit(cap, c_17, 1);

c = polyval(p_11, log10(cap))+...
    polyval(p_12, cap)+...
    polyval(p_13, cap)+...
    polyval(p_14, cap)+...
    polyval(p_15, cap)+...
    polyval(p_16, cap)+...
    polyval(p_17, cap);
c = 1.1.*c; % 1.9, contingency, 1.8 decomissioning are 0

total = 1000.*[102500, 371400, 1358100, 2488300];
error = abs(c - total)./total;

