%% 

%I must have done something wrong with the function definitions: 

%Formatting user inputs 
if Simulate(1) == false
    pts = length(dimerxdata);%number of data points
end

if Simulate(1) == true
    pts = Simulate(2);
    noise = Simulate(3);
end
%% 
% name vertices and edges of T_G 
m = length(apred);
if m <= 7 % use first few letters of alphabet
    syms b c d e f g
    edges = [ 0 b c d e f g ];
else % use e_i 
    edges = sym('e%d', [1 m]); edges(1)=NaN; 
end

v = length(apred); % number of vertices in T_G
T = sym(zeros(m)); % T_G becomes symbolic adjacency matrix of tree

for i = 2:v
    estr{i-1} = ['e' num2str(i) ]; 
    T(i,apred(i))=edges(i); % T_G is lower triangular 
end
% T % show T_G
% mons % show mons 

% P(i) is the sum of edges from a_i to root; P(1)=0
P = sym(zeros(1,m));
for i=2:m
    Timo = T^(i-1); % paths of length i-1 from row to col
    for j=1:m
        if Timo(j,1)~=0 % AA(j,1) gives paths from aj to root a1
            P(j) = sum(children(Timo(j,1))); % sum of edges from aj to root
        end
    end
end

% create matrix C where each row is state of homo-oligomer
CC = nchoosek(repmat(1:v,[1,mons]),mons); 
C = [];
for i=1:size(CC,1)   %column length of CC
    if all(diff(CC(i,:))>=0)  %difference between elements
        C = [ C ; CC(i,:) ];
    end
end
C = unique(C,'rows'); % done

% [N,M]=meshgrid(1:m,1:m);
% temp1=cat(2,N',M');
% CC=reshape(temp1,[],2);
% for i = 1:length(CC)
%     CC(i,:)=sort(CC(i,:));
% end
% C=unique(CC,'rows');

V = size(C,1);

Q = sym(zeros(V,1));
for i = 1:V % each state of oligomer 
    for n=2:mons % all the n-way interactions   
        % calculate all ways to choose n monomer states 
        % from the i-th oligomer configuration
        ways = nchoosek(C(i,:),n); 
        for r=1:size(ways,1) % for each way 
            % calculate the product of paths 
            q = 1;
            for s=1:n 
                q = q*P(ways(r,s));
            end
            Q(i) = Q(i)+q;
        end
    end
end

for i=1:V
    Q(i) = expand(Q(i));
end

for i=1:V
    disp([ sprintf('%d',C(i,:)) ' : ' sprintf('%s',Q(i))])
end
%% create symbolic matrix of monomer equilibrium constants
mon_eq_consts = sym(size(ligand_intro));
for i = 2:length(ligand_intro)+1
    varb = str2sym(['k' char(edges(i))]);
    mon_eq_consts(i-1)=varb;
end
%%
%obtain list of parameter names to pass into fitting section 
parnames = cell.empty();
for i = 1:length(Q)
    edit = arrayfun(@char, Q(i), 'uniform', 0); %make into string from sym
    edit = regexprep(edit,'*',''); %remove * signs
    edit = regexprep(edit,'+',''); %remove + signs
    edit = split(edit).';
    edit = regexprep(edit,'(^)\s*\d*',''); %remove any digits that start a term (coefficients but not exponents)
    for j = 1:length(edit) %if the parameter does not already exist in the list, add it
        if any(strcmp(parnames,edit{j})) == false
            parnames{end+1}=edit{j};          
        end
    end
end
parnames = parnames(2:end); %take off the first element, from the definition
k = length(parnames); %number of allosteric parameters

for i = 1:length(parnames)
   foo = extractBefore(parnames{i},'^2');
   parnames{i} = replace(parnames{i},'^2',foo);
end

for i = 1:length(parnames)
    parnames{i} = str2sym(parnames{i});
end
%% 
if Simulate(1) == true && isempty(allopars)
    disp('Order:')
    celldisp(parnames)
    error('Add allosteric parameter matrix in the order specified above and run again')
end
%%
%specify the response measure, either in terms of monomer or dimer states
%for example, in the two state model, [1 0 0] ([a b c]) would indicate that
%you want to calculate the percentage of monomer in state a
%the program will convert this to [1 .5 .5 0 0 0] ([aa ab ac bb bc cc])
monnames = ['a' 'b' 'c']; %find way to make this from spanning tree
dim_resp_mat = zeros(1,length(C));
G = [];
for i=1:length(C)
    G = [G {sprintf('%d',C(i,:))}];
end
G = char(G);

if length(resp_meas)==2
    [tf1, idx1] = ismember(resp_meas(1),monnames);
    [tf2, idx2] = ismember(resp_meas(2),monnames);
    converted = [num2str(idx1) num2str(idx2)];
    for i = 1:length(G)
        if G(i,:) == converted
            idx = i;
        end
    end
    dim_resp_mat(1,idx) = 1;
end

if length(resp_meas) == 1
    for i = 1:length(C)
        [tf, idx] = ismember(resp_meas,monnames);
        repeats = sum(C(i,:)==idx);
        dim_resp_mat(1,i) = repeats/mons;
    end
end

if ismember('custom',resp_meas)==true
    temp = split(resp_meas,":");
    if length(temp)==1
        disp('Order:')
        disp(char(C+'a'-1))
        error('Add custom response variable matrix in the order above and run again')
    end
    dim_resp_mat = cell2mat(temp(2));
    dim_resp_mat = str2num(dim_resp_mat);
end
%% write the function that is characterized by the given response measure
%first, use ligand_intro to determine what pi_monomerstate's are
syms x
lig = x.^ligand_intro;
lig_coeff = mon_eq_consts.*lig;
monstate_fracs = sym(size(P));
for i = 1:length(P)
    edit = char(P(i));
    for j = 2:length(ligand_intro)+1
        edit = strrep(edit,char(edges(j)),char(lig_coeff(j-1)));
    end
    edit = strrep(edit, '+', '*');
    monstate_fracs(i) = str2sym(edit);
end
monstate_fracs(1) = 1;
monstate_fracs = monstate_fracs./sum(monstate_fracs);
%next, reformat the allosteric parameters from Q to include in the dimer state fraction matrix
dimstate_fracs = sym(zeros(size(Q.')));
for i = 1:length(Q)
    edit = arrayfun(@char, Q(i), 'uniform', 0); %make into string from sym
    edit = regexprep(edit,'*',''); %remove * signs
    edit = regexprep(edit,'+','*'); %replace + with *
    edit = split(edit).';
    for j = 1:length(edit)
        movepower=edit{j};
        foo = extractBefore(movepower,'^2');
        if movepower ~= '0'
            movepower=replace(movepower,'^2',foo);
            if isstrprop(movepower(1), 'digit') == true
                power = extractBefore(movepower,2);
                movepower = movepower(2:end);
                movepower = insertAfter(movepower,length(movepower),['^' power]);
            end
        end
        if movepower == '0'
            movepower = '1';
        end
        edit{j} = movepower;
    end
        
    edit = join(edit);
    dimstate_fracs(1,i) = str2sym(edit);
end
%now multiply in appropriate monstate_fracs and coefficients to account for
%number of ways to attain each state

CC = num2str(CC);
%CC = 
frac_coeffs = zeros(size(dimstate_fracs));
for i = 1:length(C)
    count = zeros(size(C(1,:)));
    for j = 1:m
        count(j)=sum(C(i,:)==j);
    end
    ways = factorial(m)/prod(factorial(count));
    frac_coeffs(1,i)=ways;
end
frac_coeffs=frac_coeffs/gcd(sym(ways));
dimstate_fracs = dimstate_fracs .* frac_coeffs; %multiply each fractional part by the number of ways that state can be obtained
for i = 1:length(dimstate_fracs)
    dimstate_fracs(i) = dimstate_fracs(i)*monstate_fracs(C(i,1))*monstate_fracs(C(i,2)); %multiply by each of the appropriate pi_monomer's
end
dimstate_fracs;  
%% Define the function to be fit
temp = cell2sym(parnames);
pi_sym = sum(dimstate_fracs .* dim_resp_mat) %to what extent does each state contribute to the response measure
pi_temp_fun=matlabFunction(pi_sym,'vars',{mon_eq_consts temp 'x'});
pi_fun=@(allostericpars,x)pi_temp_fun(mon_ks,allostericpars,x);

%% Could simulate data based on provided parameters...
%pi_c=@(pars,x)(kb*kc*x.^2+pars(2)*pars(1)*kb^2*kc*x.^3+pars(2)*pars(1)^2*pars(3)*kb^2*kc^2*x.^4)./(1+2*kb*x+2*kb*kc*x.^2+pars(2)*kb^2*x.^2+2*pars(2)*pars(1)*kb^2*kc*x.^3+pars(2)*pars(1)^2*pars(3)*kb^2*kc^2*x.^4);

if Simulate(1) == true 
    dimerxdata = logspace(-3,3,pts);
    pics = pi_fun(allopars,dimerxdata);
    dimerydata = pics + noise*normrnd(0,1,1,pts);
end

%% plot the data to be fit

figure 
hold on
set(gca, 'XScale', 'log')
scatter(dimerxdata,dimerydata,'*','b')
if Simulate(1) == true
    x = logspace(-3,3,1000);
    plot(x,pi_fun(allopars,x),'b')
end
%% 
base10combos = linspace(0,(2^k)-1,2^k);
Models = dec2bin(base10combos);
combos=table(Models);
%% 
options = optimset('Display','off');
dimerfits_SSQs = zeros(2^k,k+1);
AICs = zeros(2^k,1);
for i = 1:2^k
    lb = ones(1,k);
    ub = ones(1,k);
    guess = ones(1,k);
    K = 0; %number of parameters used in the fit
    for j = 1:k
        if combos.Models(i,j) == '1' %if the parameter should be fit,
            lb(1,j) = 0;             %allow it to vary within wider bounds
            ub(1,j) = 1e10;
            K = K+1;
        end
    syms dimerfit SSQ    
    [dimerfit,SSQ] = lsqcurvefit(pi_fun, guess, dimerxdata, dimerydata,lb,ub,options);
    dimerfits_SSQs(i,:) = [dimerfit,SSQ];
    AICs(i,1) = pts*log(dimerfits_SSQs(i,k+1)/pts)+2*K+((2*K*(K+1))/(pts-K-1));
    end
end
%% 
x = logspace(-3,3,1000);
for i = 1:2^k
    plot(x,pi_fun(dimerfits_SSQs(i,1:length(allopars)),x))
end

labels = ["sample data", "original function"];
for i = 1:2^k
    labels = [labels ,Models(i,:)];
end
legend(labels,'Location','eastoutside')
hold off

%%
for i = 1:length(parnames)
    parnames{i} = char(parnames{i});
end
modeltable = table(base10combos.',Models);
modeltable.Properties.VariableNames = {'Base 10','Binary'};
parvals = array2table(dimerfits_SSQs(:,1:k));
parvals.Properties.VariableNames = parnames;
output = table(modeltable, AICs,parvals);
output.Properties.VariableNames = {'Model','AIC','Parameter Values'};

disp(output)
[minAIC,min_index] = min(output{:,2});
disp(['The best model is:  ' Models(min_index,:) ',  ' num2str(base10combos(min_index))])
disp(['With an AIC value of:  ' num2str(AICs(min_index))])
disp('And parameter values of:')
disp(output{min_index,3})
