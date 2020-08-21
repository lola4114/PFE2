function [res] = load_condition_determine( N, L, U)

% N is the number of PUEs
%L is the number of licenced RBs
%U is the number of unlicenced RBs

    if L + U >= N,
        
        res = input ( 'light load condition'),
        
    else
        
        res = input ( 'heavy load condition'),
        
    end
    
end