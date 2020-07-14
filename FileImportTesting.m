%% Parameter values from .txt file
opts = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines = [1, Inf];
opts.Delimiter = ";";
opts.VariableNames = ["Var1", "cc"];
opts.SelectedVariableNames = "cc";
opts.VariableTypes = ["string", "string"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts = setvaropts(opts, ["Var1", "cc"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "cc"], "EmptyFieldRule", "auto");

UserPars = readmatrix("MonroeProjectParameters.txt", opts);

clear opts

apred = str2num(UserPars(1)); %spanning tree, 7 is pretty much the limit for length due to concerns about testing the exponentially growing number of combinations
ligand_intro = str2num(UserPars(2)); %each reaction/edge's order of x(the ligand)
mons = str2double(UserPars(3)); %number or receptors in the oligomer
mon_ks = str2num(UserPars(4)); %equilibrium constants for the monomer
resp_meas = char(UserPars(5,1)); %which state is being measured
Simulate = str2num(UserPars(6)); %simulate data based on input parameters or use data input by user?
allopars = str2num(UserPars(7)); %if simulate=true, these are parameters used to simulate
%% Data from excel spreadsheet, which was generated from allopars [10,1,1] and noise = .01
if Simulate(1) == false
    dimerdata = readmatrix('Exceldimerdata.xlsx');
    dimerxdata = dimerdata(:,1).';
    dimerydata = dimerdata(:,2).';
end

