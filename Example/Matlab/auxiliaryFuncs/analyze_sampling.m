clear all;
clc;
close all;

% sampling
% pathname='/home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/src_eclipse/output'
pathname='/data/scc/thobi/200213_XYModel/sqrtN_64/T_1.34/run_1/output';
% /home/newton/thobi/Code/1911_XYModelBasedOnMHoefler/sqrtN_200_testrun/'
% filename='sampling.out'
filename='sampling_output_samp.m';
file=fullfile(pathname,filename);
run(file)

indexmat = zeros(length(TCF_times),length(qbin));
% rowmat = zeros(length(TCF_times),length(qbin));
% columnmat = zeros(length(TCF_times),length(qbin));
for i = 1:length(TCF_times)
%     rowmat(i,:) = matrix_row(length(TCF_times),length(qbin),i);
    for j = 1:length(qbin)
        indexmat(i,j) = find_matrix_index(length(qbin),i,j);
    end
end
% for j = 1:length(qbin)
%     columnmat(:,j) = matrix_column(length(TCF_times),length(qbin),j);
% end
% for i = 1:length(TCF_times)
%     if ( max(abs(matrix_row(length(TCF_times),length(qbin),i) - testmat(i,:))) ~= 0 )
%            fprintf('i=%i ',i); 
%     end
% end
% for j = 1:length(qbin)
%     if ( max(abs(matrix_column(length(TCF_times),length(qbin),j) - testmat(:,j))) ~= 0 )
%            fprintf('j=%i ',j); 
%     end
% end

% param_string = num2str(T,'T = %10.5e');
% test =  num2str('%3$s %2$s %1$s %2$s','A','B','C');

figure(1000)
subplot(3,2,1); semilogx(averaging_times,H); xlabel('time'); ylabel('H'); 
subplot(3,2,2); semilogx(averaging_times,temperature); xlabel('time'); ylabel('kT'); 
subplot(3,2,3); semilogx(averaging_times,M(1,:)); xlabel('time'); ylabel('M_x'); ylim([-1 1]);
subplot(3,2,4); semilogx(averaging_times,M(2,:)); xlabel('time'); ylabel('M_y'); ylim([-1 1]);
subplot(3,2,5); semilogx(averaging_times,absM); xlabel('time'); ylabel('|M|'); ylim([0 1]);
subplot(3,2,6); semilogx(averaging_times,M(1,:).^2+M(2,:).^2); xlabel('time'); ylabel('M^2'); ylim([0 1]);
sgtitle('Macroscopic quantities during runtime, static XY model');

figure(2000)
qindex = 1;
% indices = matrix_column(length(TCF_times),length(qbin),tindex);
indices = indexmat(:,qindex);
subplot(4,2,1); loglog(TCF_times,abs(real(gxx(indices)))); xlabel('time'); ylabel('gxx'); 
    hold on; loglog(TCF_times,abs(imag(gxx(indices))));
subplot(4,2,2); loglog(TCF_times,abs(gxx(indices))); xlabel('time'); ylabel('|gxx|'); 
subplot(4,2,3); loglog(TCF_times,abs(real(gyy(indices)))); xlabel('time'); ylabel('gyy'); 
    hold on; loglog(TCF_times,abs(imag(gyy(indices))));
subplot(4,2,4); loglog(TCF_times,abs(gyy(indices))); xlabel('time'); ylabel('|gyy|'); 
subplot(4,2,5); loglog(TCF_times,abs(real(gww(indices)))); xlabel('time'); ylabel('gww'); 
    hold on; loglog(TCF_times,abs(imag(gww(indices))));
subplot(4,2,6); loglog(TCF_times,abs(gww(indices))); xlabel('time'); ylabel('|gww|'); 
subplot(4,2,7); loglog(TCF_times,abs(real(gee(indices)))); xlabel('time'); ylabel('gee'); 
    hold on; loglog(TCF_times,abs(imag(gee(indices))));
subplot(4,2,8); loglog(TCF_times,abs(gee(indices))); xlabel('time'); ylabel('|gee|'); 
% sgtitle({'For q close to 0, time correlation functions of static XY model';param_string});
sgtitle('For q close to 0, time correlation functions of static XY model');

    
figure(3000)
tindex = 1;
indices = indexmat(tindex,:);
subplot(3,2,1); loglog(qbin,abs(real(gxx(indices)))); xlabel('q'); ylabel('gxx'); 
    hold on; loglog(qbin,abs(imag(gxx(indices))));
subplot(3,2,2); loglog(qbin,abs(real(gyy(indices)))); xlabel('q'); ylabel('gyy'); 
    hold on; loglog(qbin,abs(imag(gyy(indices))));
subplot(3,2,3); loglog(qbin,abs(real(gww(indices)))); xlabel('q'); ylabel('gww'); 
    hold on; loglog(qbin,abs(imag(gww(indices))));
subplot(3,2,4); loglog(qbin,abs(real(gee(indices)))); xlabel('q'); ylabel('gee'); 
    hold on; loglog(qbin,abs(imag(gee(indices))));
subplot(3,2,5); loglog(qbin,abs(real(gxy(indices)))); xlabel('q'); ylabel('gxy'); 
    hold on; loglog(qbin,abs(imag(gxy(indices))));
sgtitle('Susceptibilities (TCF at t = 0) for static XY model');
    
    
    
    
    
    
    