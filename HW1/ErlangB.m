%% Function of Erlang-B

%Blocking Prob. B(p,m) = (p^m) / ((m!)*sum(p^k/k!))  where k = 0...m

function B = ErlangB(p,m)

B = 1 %initial value of blocking rate
for k = 1:m
    
    B = (B.*p)./(B.*p+k);
end