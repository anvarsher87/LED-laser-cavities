function abs_length=absorption_length_cm(abs_coef)

%   global max_abs_length


 
 if abs_coef>0 %then goto 115;
     
     %{*************** ABSORPTION *****************}
     r1=rand(1);
     if r1>0; r1=log(r1); else; r1=-1.0e+05; end
     
     if abs_coef<1.0e-10; abs_coef=1.0e-10; end
     r1=-r1/(abs_coef);
     
     abs_length=r1;
     
     %{*****************************************}
 else
     abs_length=2000;  % 2000, 20 000
 end
  
    
end