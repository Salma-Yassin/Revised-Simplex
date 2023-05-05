
function [f,soln ]= revised_simplex(cj,A,b)
f = '';
Soln = '';

% Takes a LPP minimization problem in the standard form (phase I)
% Get a set of equality equation in A
% Takes the linear programming problem as a minimization problem
% Detects unbounded problem and problems with no feasible solution 

%-----------------------------------------------------------------------------------------------------

% get number of constarints and number of original constraints (before adding the artificial variables)
[m,n1] = size(A);

% Check that all RHS are non-negative
for i = 1:1:m
    if b(i)<0
        b(i) = -b(i);
        A(i,:) = -A(i,:);
    end
end

% check that a canonical form exist and thus go to phase II directly 
A_check = A(:,end-m+1:end);
if A_check == eye(m)
    disp('starting Phase II')
    [f,soln]= revised_phase2(cj,A,b); %%% Continue with the precedure of phase II 
end

% Add artificial variables to each of the m equations
% to make sure it is in canonical form and to start phase I
dj = [-sum(A) zeros(1,m)];
A = [A eye(m)];
cj = [cj zeros(1,m)];
w0 = sum(b);


% get number of constarints and total number of variables
[m,n] = size(A);

init_tablue = [A b';cj 0;dj -w0];

% initialize
basic_tablue= init_tablue(:,end-m:end); % identify the initial basic variables
soln = zeros(1,n);
basic_var = (n-m+1):1:n ;
non_basic_var = 1:1:(n-m);

while(1)
    % simplex multipliers
    pi_= - basic_tablue(end-1, 1:m);
    sigma = - basic_tablue(end, 1:m);
    
    cj_updated = cj - pi_*A;
    dj_updated = dj - sigma*A; % work only with 4 decimal places, more cause a problem 
    
    if dj_updated >= 0 & round(basic_tablue(end,end),4) > 0
        disp('There is no feasible solution to the LPP')
        return 
    elseif  dj_updated >= 0 & round(basic_tablue(end,end),4) == 0
        disp('Starting Phase II')
        
        % dropping all variables with dj_updated > 0
         A_new =[];
         C_new =[];
         index = [];
         new_basic_var = zeros(1,length(basic_var));
         new_non_basic_var = zeros(1,length(non_basic_var));
         
        for i = 1:1:length(dj_updated)
            if dj_updated(i) <= 0
                A_new =[A_new A(:,i)];
                C_new =[C_new cj(i)];
                index=[index i];
                if any(basic_var == i)
                    ind = find(basic_var == i);
                    [R,L] = size(A_new);
                    new_basic_var(ind)= L;
                else
                    ind = find(non_basic_var == i);
                    [R,L] = size(A_new);
                    new_non_basic_var(ind)= L;
                end
            end
        end
        
        new_non_basic_var=new_non_basic_var(new_non_basic_var>0);
        new_basic_var=new_basic_var(new_basic_var>0);
        
        % dropping -w row from the tablue  
        basic_tablue = basic_tablue(1:end-1,:);
        
        [f,soln_phase2] = revised_phase2(C_new,A_new,b,basic_tablue,new_basic_var,new_non_basic_var);
        soln(index) = soln_phase2;
        soln = soln(1:n1) % return only the basic variables 
        return 
    else
        
        B_inv = basic_tablue(1:end - 2 ,1:end - 1);
        
        [ms,is]= min(dj_updated); % ith variable is the entering variable
        cs = cj_updated(is);
        As = A(:,is);
        As_updated = B_inv * As;
        
        % ratio test to determine leaving variable
        if As_updated < 0
            disp('Unbounded Problem')
            return
        else
            [mr,ir]= ratio_test(As_updated,basic_tablue(1:end-2,end));
            % exchange entering and leaving variables
            l = basic_var(ir);
            basic_var(ir) = is;
            ind = find(non_basic_var == is);
            non_basic_var(ind)= l;
        end
        
        % pivot basic_tablue on the pivot element
        new_col = [As_updated ; cs ;ms];
        basic_tablue = [basic_tablue new_col];
        % m+2 is always the last column
        basic_tablue = pivot_table(ir , m+2 , basic_tablue);
        basic_tablue = basic_tablue(:,1:end-1); % Remove last column to start a new cycle
        
    end
end
end



