function [Re5e4,Re1e5,Re2e5,Re5e5,Re1e6] = ImportTables()
%Airfoil Data S1223-IL
Re5e4_tab = readtable('xf-s1223-il-50000-n5.csv');
Re1e5_tab = readtable('xf-s1223-il-100000-n5.csv');
Re2e5_tab = readtable('xf-s1223-il-200000-n5.csv');
Re5e5_tab = readtable('xf-s1223-il-500000-n5.csv');
Re1e6_tab = readtable('xf-s1223-il-1000000-n5.csv');

Re5e4 = table2array(Re5e4_tab);
Re1e5 = table2array(Re1e5_tab);
Re2e5 = table2array(Re2e5_tab);
Re5e5 = table2array(Re5e5_tab);
Re1e6 = table2array(Re1e6_tab);

end

