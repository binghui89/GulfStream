function TAC = ANNUALIZED_COSTf_sandia(X_grid, d_grid, data)
CCR = data.CCR;
% load depth_domain;
h = 800.*ones(1, 4);

cap_gen = 4E6.*X_grid;
cap_TL = 167000.*data.N+501558.*X_grid+329200.*d_grid;
cap_dply = 12496440+431853.*X_grid;
cap_dev = -13564308+10369850.*log(X_grid + 3.698963);
cap_moor = (475270 + 832.1.*h).*X_grid;
cap_sub = 4E5.*X_grid; % Subsystem and profit margin
% Contingency
cap_con = 0.1.*(cap_gen + cap_TL + cap_dply + cap_dev + cap_moor + cap_sub);

fix_OM = 2521577 + 136827.*X_grid;
% Insurance
fix_ins = 0.01.*(cap_gen + cap_TL + cap_dply + cap_dev + cap_moor + cap_sub + cap_con);

TAC = CCR.*(cap_gen + cap_TL + cap_moor + cap_dply + cap_dev + cap_sub...
    + cap_con)+ ...
    (fix_OM + fix_ins);
end