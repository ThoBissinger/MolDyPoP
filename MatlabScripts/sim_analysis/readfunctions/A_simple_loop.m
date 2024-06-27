basedir='/data/scc/thobi/210715_LinearTimeSampling/mxy_3.00';
infile_name='samp_Dynamics';
% outfile_name='samp_Dynamics_mag_errors';
outfile_name='samp_Dynamics_SCF_errors';
outfile_name='samp_Dynamics_jackknife_gmperp';
sqrtN_vals=[16 32 64 128 256];
L_vals=[9.25 18.5 37 74 148];
T_str = [".01" ".03" ".05" ".07" ".09" ".11" ".13" ".14" ".15" ".155" ".16" ".165" ".17" ".175" ".18" ".185" ".19" ".195" ".20" ".205" ".21" ".22" ".23" ".24" ".25" ".27" ".29" ".31" ".33" ".35" ".37" ".40" ".43" ".46" ".49" ".52"];
N_N=numel(sqrtN_vals);
N_T=numel(T_str);
xx=str2double(T_str);
yy=zeros(N_N,N_T);
yy2=yy;
for i_N = 1:N_N
    sqrtN=sqrtN_vals(i_N);
    L=L_vals(i_N);
    fprintf("%d   ",sqrtN);
    for i_T = 1:N_T
        T=T_str(i_T);
        fprintf("%s ",T);
        infile_base=sprintf('%s/sqrtN_%d/T_%s/%s',basedir,sqrtN,T,infile_name);
        outfile=sprintf('%s/sqrtN_%d/T_%s/%s',basedir,sqrtN,T,outfile_name);
        
        calculate_jackknife_cmperp(infile_base, outfile);
%         calculate_errors_etafit(infile_base, outfile,L);
%         calculate_errors_mag(infile_base, outfile);
%         
%         S=load(outfile);
%         yy(i_N,i_T)=S.SCF_mean;
%         yy2(i_N,i_T)=S.SCF_err;
%         
    end
    fprintf("\n");
end