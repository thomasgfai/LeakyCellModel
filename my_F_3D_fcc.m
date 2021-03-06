% calculate free volume in FCC lattice
function F = my_F_3D_fcc(v,RR)
a_vec = (sqrt(2)*v).^(1/3);

F = zeros(1,numel(a_vec));
for i=1:numel(a_vec)
    a = a_vec(i);
    
    if (a < 2*RR)
        F(i) = 0;
        
    elseif (a <= 2*sqrt(2)*RR)
        Nt = 8; No = 6;
        Vt = a.^3/(6*sqrt(2)); Vo = sqrt(2)*a.^3/3;
        
        solid_angle_t = acos(23/27)/(4*pi);
        solid_angle_o = 4*asin(1/3)/(4*pi);

        dihedral_angle_t = acos(1/3)/(2*pi);
        dihedral_angle_o = acos(-1/3)/(2*pi);

        double_int = pi/12*(4*(2*RR)+a).*(2*(2*RR)-a).^2;
        long_double_int = pi/12*(4*(2*RR)+sqrt(2)*a).*(2*(2*RR)-sqrt(2)*a).^2;

        w = sqrt((2*RR)^2*3*a^4-a^6);
        q1 = a^3;
        
        triple_int = w/6-3*a/2*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            +4*(2*RR)^3*myatan(a*w/(2*RR*q1));
        
        w = sqrt((2*RR)^2*4*a^4-2*a^6);
        q1 = 2*a^3;
        square_triple_int = w/6-a*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            -a*sqrt(2)/2*(2*(2*RR)^2-a^2/3)*pi/2+...
            +4*(2*RR)^3/3*(myatan(a*w/(2*RR*q1))+pi/2)+...
            +2*(2*RR)^3/3*(2*myatan(a*w/(2*RR*q1)));
        
        quad_int = 2*square_triple_int-long_double_int;
        
           F(i) = Nt*(Vt-3*4/3*pi*(2*RR).^3*solid_angle_t...
            +3*double_int*dihedral_angle_t)+No*(Vo/2 ...
            -2*4/3*pi*(2*RR).^3*solid_angle_o...
            +2*double_int*dihedral_angle_o...
            +long_double_int-2*square_triple_int+0.5*quad_int)-4*triple_int;
        
    elseif (a <= 2*sqrt(3)*RR)
        Nt = 8; No = 6;
        Vt = a.^3/(6*sqrt(2)); Vo = sqrt(2)*a.^3/3;
        
        solid_angle_t = acos(23/27)/(4*pi);
        solid_angle_o = 4*asin(1/3)/(4*pi);

        dihedral_angle_t = acos(1/3)/(2*pi);
        dihedral_angle_o = acos(-1/3)/(2*pi);

        double_int = pi/12*(4*(2*RR)+a).*(2*(2*RR)-a).^2;
        
        w = sqrt((2*RR)^2*3*a^4-a^6);
        q1 = a^3;
        
        triple_int = w/6-3*a/2*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            +4*(2*RR)^3*myatan(a*w/(2*RR*q1));
        
        F(i) = Nt*(Vt-3*4/3*pi*(2*RR).^3*solid_angle_t...
            +3*double_int*dihedral_angle_t)+No*(Vo...
            -5*4/3*pi*(2*RR).^3*solid_angle_o...
            +8*double_int*dihedral_angle_o)-16*triple_int;

    elseif (a <= 4*RR)
        Nt = 8; No = 6;
        Vt = a.^3/(6*sqrt(2)); Vo = sqrt(2)*a.^3/3;
        
        solid_angle_t = acos(23/27)/(4*pi);
        solid_angle_o = 4*asin(1/3)/(4*pi);

        dihedral_angle_t = acos(1/3)/(2*pi);
        dihedral_angle_o = acos(-1/3)/(2*pi);
        
        double_int = pi/12*(4*(2*RR)+a).*(2*(2*RR)-a).^2;
        
        F(i) = Nt*(Vt-3*4/3*pi*(2*RR).^3*solid_angle_t...
            +3*double_int*dihedral_angle_t)+No*(Vo...
            -5*4/3*pi*(2*RR).^3*solid_angle_o...
            +8*double_int*dihedral_angle_o);
    else
        Nt = 8; No = 6;
        Vt = a.^3/(6*sqrt(2)); Vo = sqrt(2)*a.^3/3;
        
        solid_angle_t = acos(23/27)/(4*pi);
        solid_angle_o = 4*asin(1/3)/(4*pi);
                
        F(i) = Nt*(Vt-3*4/3*pi*(2*RR).^3*solid_angle_t) ...
            +No*(Vo-5*4/3*pi*(2*RR).^3*solid_angle_o);

    end
end
end

function y = myatan(x)
  
   y = atan(x);
   if (y < 0)
       y = y+pi;
   end
end