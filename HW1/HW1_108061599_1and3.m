%% 
clear;
close all;
clc;

%% variable
% m : the total number of channels in the trunk
% p : the total offered traffic

Blocking_rate = [0.01 0.03 0.05 0.1];
%1-1 m 1~20
p1 = zeros(20,4);
m1 = [1:1:20]
for i = 1:20
    for j = 1:4
        while Blocking_rate(j) > ErlangB(p1(i,j) ,m1(i))
            p1(i,j) = p1(i,j) + 0.001;
        end
        p1(i,j) = p1(i,j) - 0.001;
    end
end

%1-2 m 200~220
p2 = zeros(21,4);
m2 = [200:1:220];
for x = 1:21
    for y = 1:4
         while Blocking_rate(y) > ErlangB(p2(x,y), m2(x))
             p2(x,y) = p2(x,y) + 0.001;
         end
         p2(x,y) = p2(x,y) - 0.001;
    end
end

%3
%total channels = 600
%reuse factore N = 5
% 3 operators
p3 = zeros(3,4)
Gc = zeros(3,4) %Trunking efficiency
divided_channel = [120 60 40];
for a = 1:3
    for b = 1:4
         while Blocking_rate(b) > ErlangB(p3(a,b), divided_channel(a))
             p3(a,b) = p3(a,b) + 0.001;
         end
         p3(a,b) = p3(a,b) - 0.001;
    end
    
end
%Trunking efficiency
for c = 1:3
    Gc(c,:) = p3(c,:)./divided_channel(c)
end
%Maximum offered traffic load per cell = the total offered traffic * share
%number
max_p(1,:) = p3(1,:);
max_p(2,:) = p3(2,:).*2
max_p(3,:) = p3(3,:).*3