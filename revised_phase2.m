
function [f,soln ]= revised_phase2(cj,A,b,basic_tablue,basic_var,non_basic_var) %% takes a LPP minimization problem in the standard form (phase II)
f = '';
Soln = '';
% get number of constarints and total number of variables
[m,n] = size(A);

init_tablue = [A b';cj 0];

% initialize
cost_tablue = [];

if nargin <4
    basic_tablue= init_tablue(:,end-m:end); % identify the initial basic variables and includes the constants column
    basic_var = (n-m+1):1:n;
    non_basic_var = 1:1:(n-m);
end
soln = zeros(1,n);


while(1)
    % simplex multipliers
    pi_= - basic_tablue(end, 1:m);
    cj_updated = cj - pi_*A;
    cost_tablue=[cost_tablue;cj_updated];
    B_inv = basic_tablue(1:end-1,1:end-1);
    
    
    if cj_updated >= 0
        disp('Optimal solution has been reached')
        soln(basic_var)= basic_tablue(1:end-1,end);
        f = - basic_tablue(end,end);
        return 
    else %......
        [ms,is]= min(cj_updated); % ith variable is the entering variable 
        As = A(:,is);
        As_updated = B_inv * As;
                
        % ratio test to determine leaving variable 
        if As_updated < 0
            disp('Unbounded Solution')
            return 
        else
            [mr,ir]= ratio_test(As_updated,basic_tablue(1:end-1,end));
            % exchange entering and leaving variables 
            l = basic_var(ir);
            basic_var(ir) = is;
            ind = find(non_basic_var == is);
            non_basic_var(ind)= l;
        end
        
        % pivot basic_tablue on the pivot element
        new_col = [As_updated ; ms];
        basic_tablue = [basic_tablue new_col];
        % m+2 is always the last column      
        basic_tablue = pivot_table(ir , m+2 , basic_tablue);
        basic_tablue = basic_tablue(:,1:end-1); % Remove last column to start a new cycle 
        
    end
end
end



