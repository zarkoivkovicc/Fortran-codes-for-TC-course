clear all;

%=== Input section =======================================================================
system = 'cyclopentyl';
%=== XY
XY = dlmread(strcat(system,'.xyz'));
XY = XY(3:7,2:3);
%=== Hamiltonian Matrix
H= dlmread(strcat(system,'.bo'));
H = H(2:end,:);

%=== GNUPLOT Control

LWIDTH = 20.0;         %=== Line width of C-C bonds
SCALE  = 0.6;          %=== Scaling factor for the MOs
CMINUS = 'blue';        %=== Color of -ve coeff. of MO 
CPLUS  = 'red';       %=== Color of +ve coeff. of MO
X_MAX  = 5.0;          %=== Maximum X coordinate
X_MIN  =-5.0;          %=== Minimum X coordinate
Y_MAX  = 5.0;          %=== Maximum Y coordinate
Y_MIN  =-5.0;          %=== Minimum Y coordinate

%=========================================================================================







%=== V - eigenvector, E - eigenvalue =====================================================
%=== No. of C atoms
V = dlmread(strcat(system,'.evec'));
V = V(2:end,:);
E = dlmread(strcat(system,'.eval'));
E = E(2:end,:);
N_atm = size(V,1);

file='huckel.out';
file_id=fopen(file,'w');

for i = 1: N_atm
  fprintf(file_id, ' Root: %5d Eigen value: %15.8f\n\n', i, E(i));
  fprintf(file_id, ' Eigen vector\n');
  for j = 1: N_atm
    fprintf(file_id, ' %15.8f\n', V(j,i));
  endfor
  fprintf(file_id,'\n');
endfor

fclose(file_id);


%=== Generate plot data
N_m(1:N_atm) = 0;
N_p(1:N_atm) = 0;
for i = 1: N_atm
  filename_p=strcat('plot/data/p_plot_',num2str(i),'.dat'); 
  filename_m=strcat('plot/data/m_plot_',num2str(i),'.dat'); 
  file_id_p = fopen(filename_p, 'w');
  file_id_m = fopen(filename_m, 'w');
  for j = 1: N_atm
    if ( sign(V(j,i)) == -1 )
      N_m(i) = N_m(i) + 1;
      fprintf(file_id_m, '%15.8f %15.8f %15.8f\n', XY(j,1:2), abs(V(j,i)));
    else
      N_p(i) = N_p(i) + 1;
      fprintf(file_id_p, '%15.8f %15.8f %15.8f\n', XY(j,1:2), V(j,i));
    endif
  endfor
  fclose(file_id_p);
  fclose(file_id_m);
endfor


%=== Generate GNUPLOT scripts

for ilevel = 1: N_atm
  
  filename = strcat('plot/data/plot_',num2str(ilevel),'.gnu');
  filename_p=strcat('plot/data/p_plot_',num2str(ilevel),'.dat'); 
  filename_m=strcat('plot/data/m_plot_',num2str(ilevel),'.dat'); 
  epsfilename= strcat('plot/images/plot_',num2str(ilevel),'.eps');
  file_id = fopen(filename, 'w');
  fprintf(file_id,'set term postscript enhanced color "Arial" 24\n');
  fprintf(file_id,'set output "%s"\n',epsfilename);
  fprintf(file_id,'set xrange [%15.8f:%15.8f]\n',X_MIN, X_MAX);
  fprintf(file_id,'set yrange [%15.8f:%15.8f]\n',Y_MIN, Y_MAX);
  fprintf(file_id,'set lmargin 0\n');
  fprintf(file_id,'set rmargin 0\n');
  fprintf(file_id,'set tmargin 0\n');
  fprintf(file_id,'set bmargin 0\n');
  fprintf(file_id,'unset key\n');
  fprintf(file_id,'unset tics\n');
  fprintf(file_id,'unset border\n');
  fprintf(file_id,'set size square\n');
  for i = 1: N_atm
    for j = i+1: N_atm
      if ( abs( H(i,j) != 0) )
        fprintf(file_id, "set arrow from %15.8f, %15.8f to %15.8f, %15.8f nohead lt 1 lc rgb 'black' lw %15.8f \n", XY(i,1:2), XY(j,1:2), LWIDTH);
      endif
    endfor
  endfor
  if ( N_m(ilevel) == 0 )
    fprintf(file_id, "plot '%s' u 1:2:($3*%15.8f) notitle with circles lc rgb '%s' lw %15.8f fill solid noborder\n",filename_p, SCALE, CPLUS, LWIDTH);
  elseif ( N_p(ilevel) == 0 )
    fprintf(file_id, "plot '%s' u 1:2:($3*%15.8f) notitle with circles lc rgb '%s' lw %15.8f fill solid noborder\n",filename_m, SCALE, CMINUS, LWIDTH);
  else
    fprintf(file_id, "plot '%s' u 1:2:($3*%15.8f) notitle with circles lc rgb '%s' lw %15.8f fill solid noborder, \\\n",filename_m, SCALE, CMINUS, LWIDTH);
    fprintf(file_id, "     '%s' u 1:2:($3*%15.8f) notitle with circles lc rgb '%s' lw %15.8f fill solid noborder\n",    filename_p, SCALE, CPLUS, LWIDTH);
  endif
  fclose(filename);
endfor
%=========================================================================================
