 %optimization function    
   %optimization function

function [x, y] = optimization (Pvmax, sv_min,sp_min) 
Pc = 23;
%gaussian AWGN
a = -174;
%SINR min for vehicles and pedestrian
%sv_min = 10;
%sp_min = 8;
%Pvmax = 11.5;
% channel gains array
%G = [g1 g21 gc1 g2 g12 gc2 gcB g1B g2B]
%G = rand(1,9);
G = ones(1,9);
% awgn pour ajouter le AWGN au lieu du nombre a
S1 = ((Pvmax * G(1)) / (Pvmax * G(3) + Pc * G(2) +a));
S2 =  ((Pvmax * G(4)) / (Pvmax * G(6) + Pc * G(5) +a));
S3 = ((Pc * G(7)) / (Pvmax * G(8) + Pvmax * G(9) +a ));

%F1 = ((Pvmax * G(1)) / (Pvmax * G(3) + Pc * G(2) + a));
%F2 = ((Pvmax * G(4)) / (Pvmax * G(6) + Pc * G(5) + a));
% F3 = (1+F1) * (1+F2);
 
 %F = log2(((x * G(1)) / (y * G(3) + Pc * G(2) +a)) + ((y * G(4)) / (x * G(6) + Pc * G(5) +a)));
 
 P11 = (Pc * G(7) - sp_min * (Pvmax * G(9) +a)) / (sp_min * G(8));
 P13 = (Pvmax * G(4) - sv_min * (Pc * G(6) +a)) / (sv_min * G(5));
 P12 =  (sv_min * (Pc * G(3) + Pvmax * G(2) +a)) / G(1);
P21 = (Pc * G(7) - sp_min * (Pvmax * G(8) +a)) / (sp_min * G(9));
 P22 = (Pvmax * G(1) - sv_min * (Pc * G(3) +a)) / (sv_min * G(2));
 P23 =  (sv_min * (Pc * G(6) + Pvmax * G(5) +a)) / G(4);
 
 switch (S3 >= sp_min && S1 >= sv_min && S2 >= sv_min)
         case 1
          % if (log2(((Pvmax * G(1)) / (P23 * G(3) + Pc * G(2) +a)) + ((P23 * G(4)) / (Pvmax * G(6) + Pc * G(5) +a))) > log2(((P12 * G(1)) / (Pvmax * G(3) + Pc * G(2) +a)) + ((¨Pvmax * G(4)) / (P12 * G(6) + Pc * G(5) +a))) && log2(((Pvmax * G(1)) / (P23 * G(3) + Pc * G(2) +a)) + ((P23 * G(4)) / (Pvmax * G(6) + Pc * G(5) +a))) > log2(((Pvmax * G(1)) / (Pvmax * G(3) + Pc * G(2) +a)) + ((¨Pvmax * G(4)) / (Pvmax * G(6) + Pc * G(5) +a))))
             
         x = Pvmax; y = P23;
         % end
     case 2
   %  if (log2(((P12 * G(1)) / (Pvmax * G(3) + Pc * G(2) +a)) + ((Pvmax * G(4)) / (P12 * G(6) + Pc * G(5) +a))) > log2(((Pvmax * G(1)) / (P23 * G(3) + Pc * G(2) +a)) + ((¨P23 * G(4)) / (Pvmax * G(6) + Pc * G(5) +a))) && log2(((P12 * G(1)) / (Pvmax * G(3) + Pc * G(2) +a)) + ((Pvmax * G(4)) / (P12 * G(6) + Pc * G(5) +a))) > log2(((Pvmax * G(1)) / (Pvmax * G(3) + Pc * G(2) +a)) + ((¨Pvmax * G(4)) / (Pvmax * G(6) + Pc * G(5) +a))))
         x = P12; y = Pvmax;
    % end
     case 3
   % if (log2(((Pvmax * G(1)) / (Pvmax * G(3) + Pc * G(2) +a)) + ((Pvmax * G(4)) / (Pvmax * G(6) + Pc * G(5) +a))) > log2(((Pvmax * G(1)) / (P23 * G(3) + Pc * G(2) +a)) + ((¨P23 * G(4)) / (Pvmax * G(6) + Pc * G(5) +a))) && log2(((Pvmax * G(1)) / (Pvmax * G(3) + Pc * G(2) +a)) + ((Pvmax * G(4)) / (Pvmax * G(6) + Pc * G(5) +a))) > log2(((P12 * G(1)) / (Pvmax * G(3) + Pc * G(2) +a)) + ((¨Pvmax * G(4)) / (Pvmax * G(6) + Pc * G(5) +a))))
         x = Pvmax; y = Pvmax;
          % end
 end
 switch (S3 < sp_min && S1 >= sv_min && S2 >= sv_min)
     case 0
           x = P11; y = Pvmax;
           case 1
         x = Pvmax; y = P21;
            case 2
         x = Pvmax; y = P22;
          case 3
          x = Pvmax; y = P23;
     case 4
               x = P12; y = Pvmax;
     case 5
            x = P13; y = Pvmax;
 end      
   
     switch (S3 >= sp_min && S1 < sv_min && S2 >= sv_min)
         case 1
         x = Pvmax; y = P21;
     case 2
         x = Pvmax; y = P22;
     case 3
         x = Pvmax; y = P23;
     end
      switch (S3 >= sp_min && S1 < sv_min && S2 < sv_min)
                 case 1
         x = P13; y = Pvmax;
     case 2
         x = P12; y = Pvmax;
     case 3
         x = P11; y = Pvmax;
     end 
 
 
 
    %end 
 %elseif (S3 < sp_min && S1 >= sv_min && S2 >= sv_min)
   %  x == P14 && y == Pvmax || x ==Pvmax && y == P21 || x == Pvmax && y == P22 || x == Pvmax && y == P23 || x == P13 && y ==  Pvmax ||  x == P12 && y ==  Pvmax;     
   %switch(S3 >= sp_min && S1 < sv_min && S2 >= sv_min)
   %  case 0
      %   x = Pvmax; y = P21;
     %case 1
     %    x = Pvmax; y = P22;
    % case 2
     %    x = Pvmax; y = P23;
   %end 
% elseif (S3 >= sp_min && S1 < sv_min && S2 >= sv_min)
 %   x == Pvmax && y ==  P21 || x == Pvmax && y == P22 || x == Pvmax && y== P23;
%    switch(S3 >= sp_min && S1 < sv_min && S2 < sv_min)
    % case 0
       %  x = P13; y = Pvmax;
    % case 1
      %   x = P12; y = Pvmax;
  %   case 2
      %   x = P11; y = Pvmax;
  % end 
% elseif (S3 >= sp_min && S1 < sv_min && S2 < sv_min)
 %  x == P13 && y == Pvmax ||  x == P12 && y == Pvmax ||  x == P11 && y == Pvmax; 
    %  if ( [x,y] == argmax(F))
       %   disp(['P1 = ', num2str(x) , .....
              %    'P2 = ', num2str(y)]);
    %  end
% end
 %x= 0; y=0;
end 
   
 