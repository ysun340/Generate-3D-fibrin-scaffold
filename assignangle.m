function [output111,output222,output33]= assignangle(my1, my2, my3,output22222,Leq,gelsize, H,angles_counter,new_bonds_counter_2)

output22 = my2;
output11 = my1;
output33 = my3;

Lx = gelsize;
Ly = H;
Lz = gelsize;

vecArr = [];
number_divisions = 0;
added_points_counter_2 = 0;
output2=[];

angle_type = 2;
atom_type = 3;
mol_type = 3;
bond_type = 3;

output222 = output22222;
output111 = output11;

for i=1:length(output22(:,1))
    number_divisions = round( output22(i,5) / Leq );
    vecArr = [vecArr; number_divisions];
    if ( number_divisions <= 1)
        new_bonds_counter_2 = new_bonds_counter_2 + 1;
        output222(new_bonds_counter_2,1) = new_bonds_counter_2;
        output222(new_bonds_counter_2,2) = bond_type;
        output222(new_bonds_counter_2,3) = output22(i,3);
        output222(new_bonds_counter_2,4) = output22(i,4);
        output222(new_bonds_counter_2,5) = output22(i,5);    
        continue;
    end
            
    for j=1:number_divisions
        
        %New ponits

        if (j < number_divisions)
           added_points_counter_2 = added_points_counter_2 +1;
        
            output111(length(output11(:,1)) + added_points_counter_2,1) = length(output11(:,1)) + added_points_counter_2;
            output111(length(output11(:,1)) + added_points_counter_2,2) = atom_type;   
                      
            if ( abs( output11(output22(i,3),3) - output11(output22(i,4),3) ) > (Lx/2) )            
            
                xdis = output11(output22(i,3),3) ...
                       + j * sign(output11(output22(i,3),3)) * ( Lx - abs( output11(output22(i,3),3) - output11(output22(i,4),3) ) ) / number_divisions;
            
                if ( xdis < Lx/2 && xdis > -Lx/2)
                    output111(length(output11(:,1)) + added_points_counter_2,3) = xdis;
                else
                    output111(length(output11(:,1)) + added_points_counter_2,3) = xdis - sign(xdis) * Lx;
                end
            else
                output111(length(output11(:,1)) + added_points_counter_2,3) = output11(output22(i,3),3) ...
                                                                      + j * ( output11(output22(i,4),3) - output11(output22(i,3),3) ) / number_divisions;
            end
         
        
         if ( abs( output11(output22(i,3),4) - output11(output22(i,4),4) ) > (Ly/2) )            
            
                ydis = output11(output22(i,3),4)...
                       + j * sign(output11(output22(i,3),4)) *( Ly - abs ( output11(output22(i,3),4) - output11(output22(i,4),4) ) ) / number_divisions;
               
              if ( ydis < Ly/2 && ydis > -Ly/2)
                   output111(length(output11(:,1)) + added_points_counter_2,4) = ydis;
              else
                   output111(length(output11(:,1)) + added_points_counter_2,4) = ydis - sign(ydis) * Ly;
              end
         else
                output111(length(output11(:,1)) + added_points_counter_2,4) = output11(output22(i,3),4) ...
                                                                      + j * ( output11(output22(i,4),4) - output11(output22(i,3),4) ) / number_divisions;
         end
        
        
          if ( abs( output11(output22(i,3),5) - output11(output22(i,4),5) ) > (Lz/2) )            
            
               zdis = output11(output22(i,3),5)...
                     + j * sign(output11(output22(i,3),5)) * ( Lz - abs ( output11(output22(i,3),5) - output11(output22(i,4),5) ) ) / number_divisions;
            
              if ( zdis < Lz/2 && zdis > -Lz/2)
                 output111(length(output11(:,1)) + added_points_counter_2,5) = zdis;
              else
                 output111(length(output11(:,1)) + added_points_counter_2,5) = zdis - sign(zdis) * Lz;
              end
          else
              output111(length(output11(:,1)) + added_points_counter_2,5) = output11(output22(i,3),5) ...
                                                                      + j * ( output11(output22(i,4),5) - output11(output22(i,3),5) ) / number_divisions;
         end 
        end
              
        output111(length(output11(:,1)) + added_points_counter_2,6) = mol_type;
        
        
        %New bonds
        if (j == 1)
            new_bonds_counter_2 = new_bonds_counter_2 + 1;
            output222(new_bonds_counter_2,1) = new_bonds_counter_2;
            output222(new_bonds_counter_2,2) = bond_type;
            output222(new_bonds_counter_2,3) = output22(i,3);
            output222(new_bonds_counter_2,4) = length( output11(:,1) ) + added_points_counter_2;
            output222(new_bonds_counter_2,5) = sqrt ( ( min( abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)), Lx - abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)) ) )^2 ...
                                                 + ( min( abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)), Ly - abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)) ) )^2 ...
                                                 + ( min( abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)), Lz - abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)) ) )^2 );
        elseif (j == number_divisions)
            new_bonds_counter_2 = new_bonds_counter_2 + 1;
            output222(new_bonds_counter_2,1) = new_bonds_counter_2;
            output222(new_bonds_counter_2,2) = bond_type;
            output222(new_bonds_counter_2,3) = length( output11(:,1) ) + added_points_counter_2;
            output222(new_bonds_counter_2,4) = output22(i,4);
            output222(new_bonds_counter_2,5) = sqrt ( ( min( abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)), Lx - abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)) ) )^2 ...
                                                 + ( min( abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)), Ly - abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)) ) )^2 ...
                                                 + ( min( abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)), Lz - abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)) ) )^2 );
        else
            new_bonds_counter_2 = new_bonds_counter_2 + 1;
            output222(new_bonds_counter_2,1) = new_bonds_counter_2;
            output222(new_bonds_counter_2,2) = bond_type;
            output222(new_bonds_counter_2,3) = length( output11(:,1) ) + added_points_counter_2 - 1;
            output222(new_bonds_counter_2,4) = length( output11(:,1) ) + added_points_counter_2;
            output222(new_bonds_counter_2,5) = sqrt ( ( min( abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)), Lx - abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)) ) )^2 ...
                                                 + ( min( abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)), Ly - abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)) ) )^2 ...
                                                 + ( min( abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)), Lz - abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)) ) )^2 );
         end
        
        %New angles
        if ( (j == 1) && (number_divisions == 2) )
            angles_counter = angles_counter + 1;
            output33(angles_counter,1) = angles_counter;
            output33(angles_counter,2) = angle_type;
            output33(angles_counter,3) = output22(i,3);
            output33(angles_counter,4) = length( output11(:,1) ) + added_points_counter_2;
            output33(angles_counter,5) = output22(i,4);
        elseif ( (j == 1) && (number_divisions ~= 2) )
            angles_counter = angles_counter + 1;
            output33(angles_counter,1) = angles_counter;
            output33(angles_counter,2) = angle_type;
            output33(angles_counter,3) = output22(i,3);
            output33(angles_counter,4) = length( output11(:,1) ) + added_points_counter_2;
            output33(angles_counter,5) = length( output11(:,1) ) + added_points_counter_2 + 1;
        elseif ( (j == number_divisions - 1) && (number_divisions ~= 2) )
            angles_counter = angles_counter + 1;
            output33(angles_counter,1) = angles_counter;
            output33(angles_counter,2) = angle_type;
            output33(angles_counter,3) = length( output11(:,1) ) + added_points_counter_2 - 1;
            output33(angles_counter,4) = length( output11(:,1) ) + added_points_counter_2;
            output33(angles_counter,5) = output22(i,4);
        elseif ( (1 < j) && ( j < number_divisions - 1) )
            angles_counter = angles_counter + 1;
            output33(angles_counter,1) = angles_counter;
            output33(angles_counter,2) = angle_type;
            output33(angles_counter,3) = length( output11(:,1) ) + added_points_counter_2 - 1;
            output33(angles_counter,4) = length( output11(:,1) ) + added_points_counter_2;
            output33(angles_counter,5) = length( output11(:,1) ) + added_points_counter_2 + 1;
                                
        end         
        
    end  
    
 end

    