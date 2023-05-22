clear all
close all
clc

Lx =2500;
Ly = Lx;
Lz = Lx;

Box.x.lo = -Lx / 2.;
Box.x.hi = Lx / 2.;
Box.y.lo = -Ly / 2.;
Box.y.hi = Ly / 2.;
Box.z.lo = -Lz / 2.;
Box.z.hi = Lz / 2.;

L.x = Box.x.hi - Box.x.lo;L.y = Box.y.hi - Box.y.lo;L.z = Box.z.hi - Box.z.lo;

Rho = 3;
periodic_x = 0;periodic_y = 0;periodic_z = 0;
Nmap = 1;Smap = 2;
% Map = zeros(ceil(Smap*L.x+1), ceil(Smap*L.y+1), ceil(Smap*L.z+1));
MatomOne = 1;Mone = 1;
num_cluster = 200;R_particle = 0.209;
post_length = 40;post_inter_layer = 10;post_outer_layer = 11;post_angle = -0*pi/180;num_post = 1;
[Natom Matom Tmol Tatom Cx Cy Cz Tbond N1bond N2bond Tang N1ang N2ang N3ang] = create_arrays();

%%%%%%%%%%%%%% add gel_shell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gelsize = Lx/2;
gelsize = Lx./2;Initial_Network_Rho = 0.1;
Leq = 1;

H = gelsize;
r1 = H - 1;
psize = 0.75;
atom_type = 2;
mol_type = 2;
bond_type = 2;
angle_type = 2;

% platelet data  (mol_type(3) atom_type(3) x y z)
% crosslink data (id mol_type(2) atom_type(2) x y z)

ave_connectivity = 5;
r_cutoff =44;
r_cutoff_inner = 15;

%% new crosslink and pdaata
clotsize_ref = 200;
num_cl_ref   = 6000;  

% set up ----  sphere clot (fib )
clot_diameter = 120 ;  %(dpd)     
clot_sidelength = clot_diameter ; %(dpd) 
box_sidelength = clot_sidelength+8;

sprintf('clot diameter=%ddpd, box size=%ddpd,c_cl=%d, c_plt=%d, c_rbc=%d',clot_diameter,box_sidelength,num_cl_ref)

cubicV = clot_sidelength^3 ;
boxV = box_sidelength^3 ;

% calculating #fib crosslink, and generating location
num_CL= ceil(num_cl_ref *cubicV/(clotsize_ref^3));
xup = clot_sidelength /2; xlow =- clot_sidelength/2;   yup = clot_sidelength /2; ylow = -clot_sidelength/2;    zup = clot_sidelength/2;  zlow = -clot_sidelength/2;
new_crosslinkdata=zeros(num_CL,6);
new_crosslinkdata(:,2)=2; new_crosslinkdata(:,3)=2;
new_crosslinkdata(:,4)= xlow +(xup-xlow).*rand(num_CL,1);  new_crosslinkdata(:,5)= ylow +(yup-ylow).*rand(num_CL,1);   new_crosslinkdata(:,6)= zlow +(zup-zlow).*rand(num_CL,1);

cal1 = sqrt(new_crosslinkdata(:,4).^2+new_crosslinkdata(:,5).^2 + new_crosslinkdata(:,6).^2  );
todelete = find(cal1>(clot_diameter/2));          new_crosslinkdata(todelete,:)=[];
new_crosslinkdata(:,1)=1:length(new_crosslinkdata(:,1));
crosslinkData=new_crosslinkdata;
num_CL=length(crosslinkData(:,1));


% figure
% scatter3(crosslinkData(:,4),crosslinkData(:,5),crosslinkData(:,6),'*b')
% hold on
% axis equal


%%
[ NatomN, MatomN, TmolN, TatomN, CxN, CyN, CzN, TbondN, N1bondN, N2bondN, TangN, N1angN, N2angN, N3angN, Spacing_nodes, FilamentDistribution] ...
    = membrane_gel_folding_sphere( Smap, Box, Initial_Network_Rho, Leq, ave_connectivity, r_cutoff, H, r1, MatomOne, periodic_x, periodic_y, periodic_z, ...
    atom_type, mol_type, bond_type, angle_type, gelsize, crosslinkData,r_cutoff_inner);

% [Map, NatomN, MatomN, TmolN, TatomN, CxN, CyN, CzN, TbondN, N1bondN, N2bondN, TangN, N1angN, N2angN, N3angN, FilamentDistribution, Spacing_nodes] ...
%     = membrane_gel_folding_sphereNEW_ANC(Map, Smap, Box, Initial_Network_Rho, Leq, ave_connectivity, r_cutoff, H, r1, MatomOne, periodic_x, periodic_y, periodic_z, atom_type, mol_type, bond_type, angle_type, psize, gelsize);

[Natom Matom Tmol Tatom Cx Cy, Cz Tbond N1bond N2bond Tang N1ang N2ang N3ang]...
    = add_arrays( Natom, Matom, Tmol, Tatom, Cx, Cy, Cz, Tbond, N1bond, N2bond, Tang, N1ang, N2ang, N3ang,...
    NatomN, MatomN, TmolN, TatomN, CxN, CyN, CzN, TbondN, N1bondN, N2bondN, TangN, N1angN, N2angN, N3angN);

%%

FileNameNew = sprintf('clot_%d_%d_%g_%g.dat',Leq, ave_connectivity, Initial_Network_Rho,r_cutoff);

create_output_file_new(FileNameNew, Box, Natom, Matom, Tmol, Tatom, Cx, Cy, Cz, Tbond, N1bond, N2bond, Tang, N1ang, N2ang, N3ang,box_sidelength);

%plot_atoms(Box, Natom, Matom, Tmol, Tatom, Cx, Cy, Cz, Tbond, N1bond, N2bond, Tang, N1ang, N2ang, N3ang);
