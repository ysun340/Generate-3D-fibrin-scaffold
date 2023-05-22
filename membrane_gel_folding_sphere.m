function [ Natom, Matom, Tmol, Tatom, Cx, Cy, Cz, Tbond, N1bond, N2bond, Tang, N1ang, N2ang, N3ang, output111, output11, output1, output1111, output2222, output222] = ...
    membrane_gel_folding_sphere( Smap, Box, initial_number_points, Leq, ave_connectivity, r_cutoff, H, ~, MatomOne, periodic_x, periodic_y, periodic_z, atom_type, mol_type, bond_type, angle_type, gelsize,...
    crosslinkData,r_cutoff_inner)

Cx = [];
Cy = [];
Cz = [];
Natom = [];
Matom = [];
Tmol = [];
Tatom = [];
Tbond = [];
N1bond = [];
N2bond = [];
Tang = [];
N1ang = [];
N2ang = [];
N3ang = [];

% cutoff distance - periodic - more organized - new point, bonds and angles
L.x = gelsize;L.y = H;L.z = gelsize;
% L.x = gelsize;L.y = H;L.z = gelsize;

%initial_number_points = round(initial_number_points * L.x * L.y * L.z);
initial_number_points = length(crosslinkData);

output1 = [];           %initial points

% initial_number_points = round ( (L.x * L.y * L.z) * Initial_Network_Rho );

% output1 = textread(sprintf('nodes_%d_1.dat', initial_number_points));
% output1 = textread(sprintf('nodes_%d_1.dat', 1000));
% output1(:,2) = atom_type;
% output1(:,6) = mol_type;
output1(:,1) = crosslinkData(:,1);
output1(:,2) = crosslinkData(:,3);
output1(:,3) = crosslinkData(:,4);
output1(:,4) = crosslinkData(:,5);
output1(:,5) = crosslinkData(:,6);
output1(:,6) = crosslinkData(:,2);

output11 = [];
output111 = [];          %final points
output2 = [];           %initial bonds
output22 = [];
output222 = [];          %final bonds
output33 = [];          %angles

initial_connections = [];
global_bonds_counter = 0;
local_bonds_counter = zeros(length(output1(:,1)),1);
distance = 0;
number_divisions = 0;
added_points_counter_1 = 0;
added_points_counter_2 = 0;
new_bonds_counter_1 = 0;
new_bonds_counter_2 = 0;
angles_counter = 0;

rg_molecule = mol_type;

for i=1:length(output1(:,1))
    if ( local_bonds_counter(i,1) >= ave_connectivity )
        continue;
    end

    local_connections = zeros ( (ave_connectivity - local_bonds_counter(i,1)) , 2 );
    local_connections(:,1) = local_connections(:,1) + 2 * sqrt(L.x^2+L.y^2+L.z^2) ;

    for j=1:length(output1(:,1))
        if (j <= i)
            continue;
        end

        [C,I] = max ( local_connections(:,1) );

        distance = sqrt ( ( min( abs(output1(i,3) - output1(j,3)), L.x - abs(output1(i,3) - output1(j,3)) ) )^2 ...
            + ( min( abs(output1(i,4) - output1(j,4)), L.y - abs(output1(i,4) - output1(j,4)) ) )^2 ...
            + ( min( abs(output1(i,5) - output1(j,5)), L.z - abs(output1(i,5) - output1(j,5)) ) )^2 );

        if ( ( distance < C ) && ( local_bonds_counter(j,1) < ave_connectivity ) && ( distance < r_cutoff )&& ( distance > r_cutoff_inner))
            local_connections(I,1) = distance;
            local_connections(I,2) = j;
        end
    end
    initial_connections ( i , ( local_bonds_counter(i,1) + 3 ) : ave_connectivity + 2) = local_connections(:,2)';

    for k=1:length( local_connections(:,2) )
        if ( local_connections(k,2) == 0 )
            continue;
        end
        global_bonds_counter = global_bonds_counter + 1;
        local_bonds_counter(i,1) = local_bonds_counter(i,1) +1;
        local_bonds_counter( local_connections(k,2),1) = local_bonds_counter( local_connections(k,2),1 ) + 1;
        initial_connections ( local_connections(k,2) , local_bonds_counter( local_connections(k,2) ) + 2 ) = i;
        output2( global_bonds_counter, 1) = global_bonds_counter;
        output2( global_bonds_counter, 2) = bond_type;
        output2( global_bonds_counter, 3) = i;
        output2( global_bonds_counter, 4) = local_connections(k,2);
        output2( global_bonds_counter, 5) = local_connections(k,1);
    end

end

initial_connections(:,1) = (1:length(output1(:,1)))';
initial_connections(:,2) = local_bonds_counter(:,1);
avg_fib_length = mean(output2(:,end));


for i=1:length(output1(:,1))
    output11(i,:) = output1(i,:);
end

for ii=1:length(output2(:,1))
    flagx = 0;
    flagy = 0;
    flagz = 0;
    if ( abs( output1(output2(ii,3),3) - output1( output2(ii,4),3) ) > L.x/2 )
        flagx = 1;
    end
    if ( abs( output1(output2(ii,3),4) - output1( output2(ii,4),4) ) > L.y/2 )
        flagy = 1;
    end
    if ( abs( output1(output2(ii,3),5) - output1( output2(ii,4),5) ) > L.z/2 )
        flagz = 1;
    end


    if ( ( (flagx == 1) && (periodic_x==0) ) || ( (flagy == 1) && (periodic_y==0) ) || ( (flagz == 1) && (periodic_z==0) ) )

        added_points_counter_1 = added_points_counter_1 + 2;
        added_points = zeros(8,3);

        added_points(1,1) = output1( output2(ii,4),3 ) - sign(output1( output2(ii,4),3) ) * L.x * flagx;
        added_points(1,2) = output1( output2(ii,4),4 ) - sign(output1( output2(ii,4),4) ) * L.y * flagy;
        added_points(1,3) = output1( output2(ii,4),5 ) - sign(output1( output2(ii,4),5) ) * L.z * flagz;

        added_points(5,1) = output1( output2(ii,3),3 ) - sign(output1( output2(ii,3),3) ) * L.x * flagx;
        added_points(5,2) = output1( output2(ii,3),4 ) - sign(output1( output2(ii,3),4) ) * L.y * flagy;
        added_points(5,3) = output1( output2(ii,3),5 ) - sign(output1( output2(ii,3),5) ) * L.z * flagz;



        if ( (periodic_x == 0) && ( flagx == 1 ) )

            if ( abs(added_points(1,1)) > (L.x/2) )
                added_points(2,1) = sign( added_points(1,1) ) * L.x / 2;
            else
                added_points(2,1) = added_points(1,1);
            end

            added_points(2,2) = output1( output2(ii,3),4 ) + ...
                ( ( added_points(1,2) - output1( output2(ii,3),4 ) ) / ( added_points(1,1) - output1( output2(ii,3),3 ) ) ) ...
                * ( added_points(2,1) - output1( output2(ii,3),3 ) );
            added_points(2,3) = output1( output2(ii,3),5 ) + ...
                ( ( added_points(1,3) - output1( output2(ii,3),5 ) ) / ( added_points(1,1) - output1( output2(ii,3),3 ) ) ) ...
                * ( added_points(2,1) - output1( output2(ii,3),3 ) );

            if ( abs(added_points(5,1)) > (L.x/2) )
                added_points(6,1) = sign( added_points(5,1) ) * L.x / 2;
            else
                added_points(6,1) = added_points(5,1);
            end

            added_points(6,2) = output1( output2(ii,4),4 ) + ...
                ( ( added_points(5,2) - output1( output2(ii,4),4 ) ) / ( added_points(5,1) - output1( output2(ii,4),3 ) ) ) ...
                * ( added_points(6,1) - output1( output2(ii,4),3 ) );
            added_points(6,3) = output1( output2(ii,4),5 ) + ...
                ( ( added_points(5,3) - output1( output2(ii,4),5 ) ) / ( added_points(5,1) - output1( output2(ii,4),3 ) ) ) ...
                * ( added_points(6,1) - output1( output2(ii,4),3 ) );
        else
            added_points(2,1) = added_points(1,1);
            added_points(2,2) = added_points(1,2);
            added_points(2,3) = added_points(1,3);

            added_points(6,1) = added_points(5,1);
            added_points(6,2) = added_points(5,2);
            added_points(6,3) = added_points(5,3);
        end


        if ( (periodic_y == 0) && ( flagy == 1 ) )

            if ( abs(added_points(2,2)) > (L.y/2) )
                added_points(3,2) = sign( added_points(2,2) ) * L.y / 2;
            else
                added_points(3,2) = added_points(2,2);
            end

            added_points(3,1) = output1( output2(ii,3),3 ) + ...
                ( ( added_points(2,1) - output1( output2(ii,3),3 ) ) / ( added_points(2,2) - output1( output2(ii,3),4 ) ) ) ...
                * ( added_points(3,2) - output1( output2(ii,3),4 ) );
            added_points(3,3) = output1( output2(ii,3),5 ) + ...
                ( ( added_points(2,3) - output1( output2(ii,3),5 ) ) / ( added_points(2,2) - output1( output2(ii,3),4 ) ) ) ...
                * ( added_points(3,2) - output1( output2(ii,3),4 ) );

            if ( abs(added_points(6,2)) > (L.y/2) )
                added_points(7,2) = sign( added_points(6,2) ) * L.y / 2;
            else
                added_points(7,2) = added_points(6,2);
            end

            added_points(7,1) = output1( output2(ii,4),3 ) + ...
                ( ( added_points(6,1) - output1( output2(ii,4),3 ) ) / ( added_points(6,2) - output1( output2(ii,4),4 ) ) ) ...
                * ( added_points(7,2) - output1( output2(ii,4),4 ) );
            added_points(7,3) = output1( output2(ii,4),5 ) + ...
                ( ( added_points(6,3) - output1( output2(ii,4),5 ) ) / ( added_points(6,2) - output1( output2(ii,4),4 ) ) ) ...
                * ( added_points(7,2) - output1( output2(ii,4),4 ) );
        else
            added_points(3,1) = added_points(2,1);
            added_points(3,2) = added_points(2,2);
            added_points(3,3) = added_points(2,3);

            added_points(7,1) = added_points(6,1);
            added_points(7,2) = added_points(6,2);
            added_points(7,3) = added_points(6,3);
        end




        if ( (periodic_z == 0) && ( flagz == 1 ) )

            if ( abs(added_points(3,3)) > (L.z/2) )
                added_points(4,3) = sign( added_points(3,3) ) * L.z / 2;
            else
                added_points(4,3) = added_points(3,3);
            end

            added_points(4,1) = output1( output2(ii,3),3 ) + ...
                ( ( added_points(3,1) - output1( output2(ii,3),3 ) ) / ( added_points(3,3) - output1( output2(ii,3),5 ) ) ) ...
                * ( added_points(4,3) - output1( output2(ii,3),5 ) );
            added_points(4,2) = output1( output2(ii,3),4 ) + ...
                ( ( added_points(3,2) - output1( output2(ii,3),4 ) ) / ( added_points(3,3) - output1( output2(ii,3),5) ) ) ...
                * ( added_points(4,3) - output1( output2(ii,3),5 ) );

            if ( abs(added_points(7,3)) > (L.z/2) )
                added_points(8,3) = sign( added_points(7,3) ) * L.z / 2;
            else
                added_points(8,3) = added_points(7,3);
            end

            added_points(8,1) = output1( output2(ii,4),3 ) + ...
                ( ( added_points(7,1) - output1( output2(ii,4),3 ) ) / ( added_points(7,3) - output1( output2(ii,4),5 ) ) ) ...
                * ( added_points(8,3) - output1( output2(ii,4),5 ) );
            added_points(8,2) = output1( output2(ii,4),4 ) + ...
                ( ( added_points(7,2) - output1( output2(ii,4),4 ) ) / ( added_points(7,3) - output1( output2(ii,4),5 ) ) ) ...
                * ( added_points(8,3) - output1( output2(ii,4),5 ) );
        else
            added_points(4,1) = added_points(3,1);
            added_points(4,2) = added_points(3,2);
            added_points(4,3) = added_points(3,3);

            added_points(8,1) = added_points(7,1);
            added_points(8,2) = added_points(7,2);
            added_points(8,3) = added_points(7,3);
        end

        output11(length(output1(:,1)) + added_points_counter_1 - 1,1) = length(output1(:,1)) + added_points_counter_1 - 1;
        output11(length(output1(:,1)) + added_points_counter_1 - 1,2) = atom_type;
        output11(length(output1(:,1)) + added_points_counter_1 - 1,6) = mol_type;

        if ( abs(added_points(4,1)) > (L.x/2) )
            output11(length(output1(:,1)) + added_points_counter_1 - 1,3) = added_points(4,1) - sign(added_points(4,1)) * L.x;
        else
            output11(length(output1(:,1)) + added_points_counter_1 - 1,3) = added_points(4,1);
        end

        if ( abs(added_points(4,2)) > (L.y/2) )
            output11(length(output1(:,1)) + added_points_counter_1 - 1,4) = added_points(4,2) - sign(added_points(4,2)) * L.y;
        else
            output11(length(output1(:,1)) + added_points_counter_1 - 1,4) = added_points(4,2);
        end

        if ( abs(added_points(4,3)) > (L.z/2) )
            output11(length(output1(:,1)) + added_points_counter_1 - 1,5) = added_points(4,3) - sign(added_points(4,3)) * L.z;
        else
            output11(length(output1(:,1)) + added_points_counter_1 - 1,5) = added_points(4,3);
        end


        output11(length(output1(:,1)) + added_points_counter_1,1) = length(output1(:,1)) + added_points_counter_1;
        output11(length(output1(:,1)) + added_points_counter_1,2) = atom_type;
        output11(length(output1(:,1)) + added_points_counter_1,6) = mol_type;

        if ( abs(added_points(8,1)) > (L.x/2) )
            output11(length(output1(:,1)) + added_points_counter_1,3) = added_points(8,1) - sign(added_points(8,1)) * L.x;
        else
            output11(length(output1(:,1)) + added_points_counter_1,3) = added_points(8,1);
        end

        if ( abs(added_points(8,2)) > (L.y/2) )
            output11(length(output1(:,1)) + added_points_counter_1,4) = added_points(8,2) - sign(added_points(8,2)) * L.y;
        else
            output11(length(output1(:,1)) + added_points_counter_1,4) = added_points(8,2);
        end

        if ( abs(added_points(8,3)) > (L.z/2) )
            output11(length(output1(:,1)) + added_points_counter_1,5) = added_points(8,3) - sign(added_points(8,3)) * L.z;
        else
            output11(length(output1(:,1)) + added_points_counter_1,5) = added_points(8,3);
        end

        new_bonds_counter_1 = new_bonds_counter_1 + 2;
        output22(new_bonds_counter_1 - 1,1) = new_bonds_counter_1 - 1;
        output22(new_bonds_counter_1 - 1,2) = bond_type;
        output22(new_bonds_counter_1 - 1,3) = output2(ii,3);
        output22(new_bonds_counter_1 - 1,4) = length(output1(:,1)) + added_points_counter_1 - 1;
        output22(new_bonds_counter_1 - 1,5) = norm ( output1(output2(ii,3),3:5) - added_points(4,:) );

        output22(new_bonds_counter_1,1) = new_bonds_counter_1;
        output22(new_bonds_counter_1,2) = bond_type;
        output22(new_bonds_counter_1,3) = output2(ii,4);
        output22(new_bonds_counter_1,4) = length(output1(:,1)) + added_points_counter_1;
        output22(new_bonds_counter_1,5) = norm ( output1(output2(ii,4),3:5) - added_points(8,:) );

    else

        new_bonds_counter_1 = new_bonds_counter_1 + 1;
        output22(new_bonds_counter_1,1) = new_bonds_counter_1;
        output22(new_bonds_counter_1,2) = bond_type;
        output22(new_bonds_counter_1,3) = output2(ii,3);
        output22(new_bonds_counter_1,4) = output2(ii,4);
        output22(new_bonds_counter_1,5) = output2(ii,5);
    end

end


for i=1:length(output11(:,1))
    output111(i,:) = output11(i,:);
end

vecArr = [];

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

            if ( abs( output11(output22(i,3),3) - output11(output22(i,4),3) ) > (L.x/2) )

                xdis = output11(output22(i,3),3) ...
                    + j * sign(output11(output22(i,3),3)) * ( L.x - abs( output11(output22(i,3),3) - output11(output22(i,4),3) ) ) / number_divisions;

                if ( xdis < L.x/2 && xdis > -L.x/2)
                    output111(length(output11(:,1)) + added_points_counter_2,3) = xdis;
                else
                    output111(length(output11(:,1)) + added_points_counter_2,3) = xdis - sign(xdis) * L.x;
                end
            else
                output111(length(output11(:,1)) + added_points_counter_2,3) = output11(output22(i,3),3) ...
                    + j * ( output11(output22(i,4),3) - output11(output22(i,3),3) ) / number_divisions;
            end


            if ( abs( output11(output22(i,3),4) - output11(output22(i,4),4) ) > (L.y/2) )

                ydis = output11(output22(i,3),4)...
                    + j * sign(output11(output22(i,3),4)) *( L.y - abs ( output11(output22(i,3),4) - output11(output22(i,4),4) ) ) / number_divisions;

                if ( ydis < L.y/2 && ydis > -L.y/2)
                    output111(length(output11(:,1)) + added_points_counter_2,4) = ydis;
                else
                    output111(length(output11(:,1)) + added_points_counter_2,4) = ydis - sign(ydis) * L.y;
                end
            else
                output111(length(output11(:,1)) + added_points_counter_2,4) = output11(output22(i,3),4) ...
                    + j * ( output11(output22(i,4),4) - output11(output22(i,3),4) ) / number_divisions;
            end


            if ( abs( output11(output22(i,3),5) - output11(output22(i,4),5) ) > (L.z/2) )

                zdis = output11(output22(i,3),5)...
                    + j * sign(output11(output22(i,3),5)) * ( L.z - abs ( output11(output22(i,3),5) - output11(output22(i,4),5) ) ) / number_divisions;

                if ( zdis < L.z/2 && zdis > -L.z/2)
                    output111(length(output11(:,1)) + added_points_counter_2,5) = zdis;
                else
                    output111(length(output11(:,1)) + added_points_counter_2,5) = zdis - sign(zdis) * L.z;
                end
            else
                output111(length(output11(:,1)) + added_points_counter_2,5) = output11(output22(i,3),5) ...
                    + j * ( output11(output22(i,4),5) - output11(output22(i,3),5) ) / number_divisions;
            end

            if j==1
                rg_molecule = rg_molecule+1;
            end
        end

        if number_divisions ==1;
            output111(length(output11(:,1)) + added_points_counter_2,6) = mol_type;
        else
            output111(length(output11(:,1)) + added_points_counter_2,6) = rg_molecule;
        end


        %New bonds
        if (j == 1)
            new_bonds_counter_2 = new_bonds_counter_2 + 1;
            output222(new_bonds_counter_2,1) = new_bonds_counter_2;
            output222(new_bonds_counter_2,2) = bond_type;
            output222(new_bonds_counter_2,3) = output22(i,3);
            output222(new_bonds_counter_2,4) = length( output11(:,1) ) + added_points_counter_2;
            output222(new_bonds_counter_2,5) = sqrt ( ( min( abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)), L.x - abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)) ) )^2 ...
                + ( min( abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)), L.y - abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)) ) )^2 ...
                + ( min( abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)), L.z - abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)) ) )^2 );
        elseif (j == number_divisions)
            new_bonds_counter_2 = new_bonds_counter_2 + 1;
            output222(new_bonds_counter_2,1) = new_bonds_counter_2;
            output222(new_bonds_counter_2,2) = bond_type;
            output222(new_bonds_counter_2,3) = length( output11(:,1) ) + added_points_counter_2;
            output222(new_bonds_counter_2,4) = output22(i,4);
            output222(new_bonds_counter_2,5) = sqrt ( ( min( abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)), L.x - abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)) ) )^2 ...
                + ( min( abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)), L.y - abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)) ) )^2 ...
                + ( min( abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)), L.z - abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)) ) )^2 );
        else
            new_bonds_counter_2 = new_bonds_counter_2 + 1;
            output222(new_bonds_counter_2,1) = new_bonds_counter_2;
            output222(new_bonds_counter_2,2) = bond_type;
            output222(new_bonds_counter_2,3) = length( output11(:,1) ) + added_points_counter_2 - 1;
            output222(new_bonds_counter_2,4) = length( output11(:,1) ) + added_points_counter_2;
            output222(new_bonds_counter_2,5) = sqrt ( ( min( abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)), L.x - abs(output111(output222(new_bonds_counter_2,3),3) - output111(output222(new_bonds_counter_2,4),3)) ) )^2 ...
                + ( min( abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)), L.y - abs(output111(output222(new_bonds_counter_2,3),4) - output111(output222(new_bonds_counter_2,4),4)) ) )^2 ...
                + ( min( abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)), L.z - abs(output111(output222(new_bonds_counter_2,3),5) - output111(output222(new_bonds_counter_2,4),5)) ) )^2 );
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

atom_flag = zeros ( length(output111(:,1)),2 );
bond_flag = zeros ( length(output222(:,1)),1 );
angle_flag = zeros ( length(output33(:,1)),1 );
atom_flag_counter = 0;
bond_flag_counter = 0;
angle_flag_counter = 0;

output1111 = [];
output11111 = []; % I added this for DFS
output2222 = [];
output22222 = [];
output33333 = [];
output333 = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BELOW FOR DFS TEST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(output111(:,1))
    %     if ( sqrt( output111(i,3)^2 + output111(i,5)^2 + output111(i,4)^2) - (L.x./2) >= 0)
    %         continue;
    %     end
    atom_flag_counter = atom_flag_counter + 1;
    atom_flag(i,1) = 1;
    atom_flag(i,2) = atom_flag_counter;
    output11111(atom_flag_counter,1) = atom_flag_counter;
    output11111(atom_flag_counter,3) = output111(i,3);
    output11111(atom_flag_counter,4) = output111(i,4);
    output11111(atom_flag_counter,5) = output111(i,5);
    output11111(atom_flag_counter,6) = mol_type;
    output11111(atom_flag_counter,2) = atom_type;
end

for i=1:length(output222(:,1))
    if ( (atom_flag(output222(i,3),1) == 1) && ( atom_flag(output222(i,4),1) == 1 ) )
        bond_flag(i) = 1;
        bond_flag_counter = bond_flag_counter + 1;
        output22222(bond_flag_counter,1) = bond_flag_counter;
        output22222(bond_flag_counter,2) = output222(i,2);
        output22222(bond_flag_counter,3) = atom_flag(output222(i,3),2);
        output22222(bond_flag_counter,4) = atom_flag(output222(i,4),2);
        output22222(bond_flag_counter,5) = output222(i,5);
        %         if (sqrt(output111(output22222(bond_flag_counter,3),3))^2 + output111(output22222(bond_flag_counter,3),4)^2 + output111(output22222(bond_flag_counter,3),5)^2) <= (gelsize./3);
        %             output22222(bond_flag_counter,2) = output22222(bond_flag_counter,2) + 1;
        %         end
    end
end

for i=1:length(output33(:,1))
    if ( (atom_flag(output33(i,3),1) == 1) && ( atom_flag(output33(i,4),1) == 1 ) && ( atom_flag(output33(i,5),1) == 1 ) )
        angle_flag(i) = 1;
        angle_flag_counter = angle_flag_counter + 1;
        output33333(angle_flag_counter,1) = angle_flag_counter;
        output33333(angle_flag_counter,2) = output33(i,2);
        output33333(angle_flag_counter,3) = atom_flag(output33(i,3),2);
        output33333(angle_flag_counter,4) = atom_flag(output33(i,4),2);
        output33333(angle_flag_counter,5) = atom_flag(output33(i,5),2);
    end
end

output1111 = output111;
avg_fibrinL = mean(output22(:,5));sprintf('average fibrin length=%.2fdpd',avg_fibrinL)
num_fibrin=length(output22);







%%

Cx = output1111(:,3);
Cy = output1111(:,4);
Cz = output1111(:,5);
Natom = output1111(:,1);
Matom (1:length(output1111(:,1)))= MatomOne;
Tmol = output1111(:,6);
Tatom = output1111(:,2);


Tbond = output22222(:,2);
N1bond = output22222(:,3);
N2bond = output22222(:,4);


Tang = output33333(:,2);
N1ang = output33333(:,3);
N2ang = output33333(:,4);
N3ang = output33333(:,5);




end




