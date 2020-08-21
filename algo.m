%algorithme 

Pv =23;
s1 = 10;
s2 = 8;
%a = zeros(1,2);
A = zeros(2,6);
B = zeros(2,6);
C = zeros(2,6);
D = zeros(2,6);
%sup1 = zeros(2,6);
%sup2 = zeros(2,6);
%K = zeros(3,2);
%L = zeros(3,2);
%Z = zeros(2,3);
%X = zeros(2,3);
Pc = 23;
%a = -174;
priority1 = randi (8);
priority2 = randi (8);
%G = rand(1,9);
G = ones(1,9);
for j = 1:1:2
for i=1:1:6
[x, y] = optimization(Pv,s1,s2);
F1 = ((x * G(1)) / (y * G(3) + Pc * G(2) ));%+ a));
F2 = ((y * G(4)) / (x * G(6) + Pc * G(5) ));%+ a));
F3 = ((Pc * G(7)) / (x * G(8) + y * G(9) ));%+ a)); 
 %F4 = (1+F1) * (1+F2);
 F5 = priority1 * log2(1+F1);
 F6 = priority2 * log2(1+F2);
 F = log2(1+F3);
 % matrix pour la somme des valeurs d'informations des différents couples
 % formés par toutes les combinaisons 
 
 A(j,i) = abs(F5) + abs(F6);
 %B(j,i) = abs(F6);
 
 %pour remplir le tableau des pedestriçans avec les rates
 C(j,i) = F;
 %disp(abs(F));
end
end
%for j = 1:1:2
%for i=1:1:6
  %  disp('A = ',num2str(A(i)));
   % disp('C =',num2str(C(i)));
%end 
%
for i=1:1:3 , j=6:-1:4;
    %for the first pedestrian
    if(and((A(1,i)>A(2,i)), (A(1,j) > A(2,j))))
        if(C(1,i) < C(1,j))
           B(1,i) = A(1,j);
            B(2,i) = j;
        
        elseif(C(1,i) > C(1,j))
            B(1,i) = A(1,i);
            B(2,i) = i;
        end
       %for the second pedestrian
    elseif(and((A(1,i) <A(2,i)), (A(1,j) < A(2,j))))
        if (C(2,i) < C(2,j))
            D(1,i) = A(2,j);
            D(2,i) = j;
        elseif (C(2,i) > C(2,j))
               D(1,i) = A(2,i);
            D(2,i) = i;
        end
    end
end
% plotting the total thr
lambda = 0 : 1 : 5;
delta = 0:1:5;
    plot(lambda, B(1,:), 'ks', delta,D(1,:), 'ks')

