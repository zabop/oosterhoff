clear
close all
format long
cd /home/ps738/5vardata/data/selector_g_lc/
files=dir;
n_files=length(files)
 
cd '/home/ps738/4vardata/data/selector_overall'
QSOlist = importdata('new_selectedQSOs');
cd '/home/ps738/4vardata/matlab'

number_of_qsos = length(QSOlist.textdata);
 
% G2spectrafrom=56400
% G2spectrauntil=56800
% 
% G3spectrafrom=56800
% G3spectrauntil=57100
% 
% G4spectrafrom=57100
% G4spectrauntil=57500
% 
% %for k=3:n_files
% 
for k=42:42
    cd /home/ps738/5vardata/data/selector_g_lc/
    h=figure; 
    set(gca,'fontsize',14)
    hold on
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    set(h,'visible','on');
    s=dir('*');
    file_list = {s.name}';
%     objectname_cell=file_list(k,1);
%     objectname=objectname_cell{1,1};
%     LC=importdata(objectname)
%     mjd_g=LC(:,1)
%     mag_g=LC(:,2)
%     plot(mjd_g,mag_g,'+b')
%     hold on
%     cd /home/ps738/5vardata/data/selector_r_lc/
%     s=dir('*');
%     file_list = {s.name}';
%     objectname_cell=file_list(k,1);
%     objectname=objectname_cell{1,1}
%     LC=importdata(objectname)
%     mjd_r=LC(:,1)
%     mag_r=LC(:,2)
%     plot(mjd_r,mag_r,'+g')
%     hold on
%     cd /home/ps738/5vardata/data/selector_i_lc/
%     s=dir('*');
%     file_list = {s.name}';
%     objectname_cell=file_list(k,1);
%     objectname=objectname_cell{1,1}
%     LC=importdata(objectname)
%     mjd_i=LC(:,1)
%     mag_i=LC(:,2)
%     plot(mjd_i,mag_i,'+r')
%     hold on           
%     objectname_nospaces=strtrim(objectname)
%     cd /home/ps738/5vardata/data/plottogether_wherearespeclines/
%     spectralines_mjd=importdata(objectname_nospaces)
%     spectralines_coord=spectralines_mjd.data
%     yL = get(gca,'YLim');
%     numberofspectralines=size(spectralines_coord,1)
%     for n=1:numberofspectralines
%         if ((G2spectrafrom <= spectralines_coord(n,1)) & (spectralines_coord(n,1) <= G2spectrauntil))
%             line([spectralines_coord(n,1) spectralines_coord(n,1)],yL,'Color','c','LineStyle','--');
%         end
%     end
%     
%     for n=1:numberofspectralines
%         if ((G3spectrafrom <= spectralines_coord(n,1)) & (spectralines_coord(n,1) <= G3spectrauntil))
%             line([spectralines_coord(n,1) spectralines_coord(n,1)],yL,'Color','m','LineStyle','--');
%         end
%     end
%     
%     for n=1:numberofspectralines
%         if ((G4spectrafrom <= spectralines_coord(n,1)) & (spectralines_coord(n,1) <= G4spectrauntil))
%             line([spectralines_coord(n,1) spectralines_coord(n,1)],yL,'Color','k','LineStyle','--');
%         end
%     end
%     
%     cd /home/ps738/5vardata/matlab/
%     
%     xlabel('mjd')
%     ylabel('mag')
%     mod_objectname=strrep(objectname,'_','\_')
%     title(mod_objectname)
%     legend('g band','r band','i band','Location','best')
%     
%     mod_objectname1=strtrim(objectname)
%  %   cd /home/ps738/5vardata/data/flux_and_spectra/
%  %   kk=int2str(k-2)
%  %   saveas(h,kk,'epsc')
%  %   cd /home/ps738/5vardata/matlab/
%     
%     
%   cd '/home/ps738/6vardata/data/errorcalc2' 
%   EWdata = importdata('showEWchange');
%   cd '/home/ps738/4vardata/matlab'
% 
%   for i=1:length(EWdata.data(:,1))
%           times1=EWdata.data(i,6) 
%           times=ones(3).*times1
%           timesVsEWs=[times(1), EWdata.data(i,3);times(2),EWdata.data(i,4);times(3), EWdata.data(i,5)] 
%           yyaxis right
%           plot(timesVsEWs(:,1),timesVsEWs(:,2),'.r','markersize',40)
%   end
end
