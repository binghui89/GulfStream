function MEP_combine()

load 2009MEP; MEP1 = ME_grid; clear ME_grid;
load 2010MEP; MEP2 = ME_grid; clear ME_grid;
load 2011MEP; MEP3 = ME_grid; clear ME_grid;
load 2012MEP; MEP4 = ME_grid; clear ME_grid;
load 2013MEP; MEP5 = ME_grid; clear ME_grid;
load 2014MEP; MEP6 = ME_grid; clear ME_grid;

MEP = cat(3, MEP1, MEP2, MEP3, MEP4, MEP5, MEP6);
clear MEP1 MEP2 MEP3 MEP4 MEP5 MEP6;
save('MEP.mat');
end