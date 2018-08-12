close all
clc
format long
clear
echo on

cd /home/ps738/7vardata/data/
fileID = fopen('resultsofeachbin_ext_inVarD7.txt','w');
cd /home/ps738/7vardata/matlab/

cd /home/ps738/7vardata/TR/
united=importdata('united_headerontop.txt');
cd /home/ps738/7vardata/matlab/

qsolist1 = united.textdata(2:end,1)
qsolist = string(qsolist1)
numberofqsos1=size(qsolist);
numberofqsos =numberofqsos1(1);

for k=30:numberofqsos
    
    lower_window_default=[2200 2700];
    higher_window_default=[2900 3090];
    
    [skipthisone,redshift]=datainpeakrange(lower_window_default,higher_window_default,united,qsolist,k)
    if skipthisone==0
        continue
    end
    
    if sum(united.data(k,18:25)) ~= 0
        
        fluxdatabins1=united.data(k,18:25);
        [row,col,fluxdatabins] = find(fluxdatabins1);
        fluxwindows1=zeros(length(fluxdatabins)/2,2);
        fluxwindow=fluxwindows1-99;
        lowerwindowsarray=fluxdatabins(1:2:length(fluxdatabins)-1);
        higherwindowsarray=fluxdatabins(2:2:length(fluxdatabins));
        fluxwindow=[lowerwindowsarray' higherwindowsarray']
        clear fluxwindows1
        
        numberoffluxwindows=length(fluxwindow(1:end,1))    
        for w=1:numberoffluxwindows
            spectrafrom=fluxwindow(w,1)
            spectrauntil=fluxwindow(w,2)
            
            ozdes_path='/data/desardata4/OzDES/Spectra/170130/';
            specname=qsolist(k)
            info=fitsinfo(char(strcat(ozdes_path,specname,'.fits')));
            next=length(info.Image);
            nepochs=(next-3)/2;
            
            spectot=fitsread(char(strcat(ozdes_path,specname,'.fits')),'Primary');
            errtot=fitsread(char(strcat(ozdes_path,specname,'.fits')),'Image',1);
            npix=length(spectot);
            pix=[1:npix];
            
            iwav1=find(ismember(info.PrimaryData.Keywords(:,1),'CRVAL1'));
            wav1=cell2mat(info.PrimaryData.Keywords(iwav1,2));
            ipix1=find(ismember(info.PrimaryData.Keywords(:,1),'CRPIX1'));
            pix1=cell2mat(info.PrimaryData.Keywords(ipix1,2));
            idelta=find(ismember(info.PrimaryData.Keywords(:,1),'CDELT1'));
            deltaw=cell2mat(info.PrimaryData.Keywords(idelta,2));
            
            wav=wav1-(pix1-pix)*deltaw;
            wav_rest=wav./(1+redshift);
            
            WindowSpectot=NaN(nepochs,5000);
            
            spectracounter=0
            for j=1:nepochs
                spec_here=fitsread(char(strcat(ozdes_path,specname,'.fits')),'Image',2*j);
                %err_here=fitsread(char(strcat(ozdes_path,specname,'.fits')),'Image',2*j+1);
                keywords_here=info.Image(2*j).Keywords;
                imjd=find(ismember(keywords_here(:,1),'UTMJD'));
                mjd=keywords_here(imjd,2);
                mjd_value=cell2mat(mjd);
                NUMBEROFNANSINSPEC = sum(isnan(spec_here));
                
                s = ((spectrafrom <= mjd_value) & (mjd_value <= spectrauntil));
                if s==1
                    spectracounter=spectracounter+1
                    WindowSpectot(j,:)=spec_here;
                    mjds(spectracounter,1)=mjd_value
                end
                
            end
            meanmjd=mean(mjds)
            clear mjds
            
            
            
            %% HERE IS THE COAD SPECTRA AND ITS ERROR IN THE MJD REGION SPECIFIED IN THE WTH ROW OF FLUXWINDOW MATRIX
            meanspec   = nanmean(WindowSpectot,1);
            stdspec = nanstd(WindowSpectot,1); %THIS ERROR DOESNT SEEM GOOD, IT DOESNT USE ERR_HERE
            %% EW WITH DEFAULT WINDOWS
            [continuum, transformedspace_flambda_original, zeros1]=totransformedspace(lower_window_default,higher_window_default,wav_rest,meanspec)
            [a1, b1, c1, x,gof]=equivalentwidth(lower_window_default, higher_window_default, wav_rest,transformedspace_flambda_original,1);
            func1 = @(x) a1*exp(-((x-b1)/c1).^2);
            ew_gauss1 = integral(func1,2700,2900)
            %% EW WITH NEW WINDOWS
            [new_lower_window,new_higher_window]=decideonwindows(united,lower_window_default,higher_window_default,k)
            [new_continuum, new_transformedspace_flambda_original, new_zeros1]=totransformedspace(new_lower_window,new_higher_window,wav_rest,meanspec)
            [aa1, bb1, cc1, xx,gof]=equivalentwidth(new_lower_window, new_higher_window, wav_rest,new_transformedspace_flambda_original,2);
            ffuunncc1 = @(xx) aa1*exp(-((xx-bb1)/cc1).^2);
            new_ew_gauss1 = integral(ffuunncc1,2700,2900)
            %% EW WITH LORENTZIAN, NEW WINDOWS
            [nanfiltered_wav_rest_MgII,lorentzianIntegralNT,params2]=equivalentwidth(new_lower_window, new_higher_window, wav_rest,new_transformedspace_flambda_original,3);
            %% PLOT DATA OF EACH BIN
            localdat = figure('units','normalized','outerposition',[0 0 1 1]);
            set(gcf,'Visible','Off')
            locd1 = plot(wav_rest,meanspec)
            hold on
            locd2 = plot(wav_rest,continuum,'--k')
            hold on
            new_locd2 = plot(wav_rest,new_continuum,'r')
            hold on
            new_locd3 = plot(wav_rest,new_transformedspace_flambda_original,'r')
            hold on
            locd3 = plot(wav_rest,transformedspace_flambda_original,':b')
            hold on
            locd4 = plot(wav_rest,zeros1,'y')
            hold on
            locd5 = plot(x,a1*exp(-((x-b1)/c1).^2),'LineWidth',3,'color','b')
            hold on
            locd6 = plot(xx,aa1*exp(-((xx-bb1)/cc1).^2),'LineWidth',3,'color','r');
            hold on
            locd7 = plot(nanfiltered_wav_rest_MgII,(2*params2(1)./pi)*(params2(3)./(4*(nanfiltered_wav_rest_MgII-params2(2)).^2 + params2(3).^2)),'g-','LineWidth',3);
            thisisfig = get(gca)
            
            title([specname,spectrafrom,spectrauntil])
            legend([locd1, locd2, new_locd2 new_locd3 locd3 locd4 locd5 locd6 locd7],'meanspec','continuum','new\_continuum','transforedspace\_neww','transformedspace\_oldw','zeros1','G1\_oldw','G2\_neww','lorentzian')
            localdatgca = gca
            
            specnamefortitle=strrep(specname,'-','_')
            %            thisistitle1=strcat(specnamefortitle,'_',num2str(spectrafrom),'_',num2str(spectrauntil))
            cd /home/ps738/7vardata/data/
            saveas(localdat,char(strcat(specnamefortitle,'_',num2str(spectrafrom),'_',num2str(spectrauntil))),'epsc')
            cd /home/ps738/7vardata/matlab/
            %            location1=strcat('/home/ps738/6vardata/data/errorcalc2/x12',thisistitle1)
            
            zoomin=copyobj(localdat,0)
            set(gcf,'Visible','Off')
            xlim([2650 2950])
            ylim([min(real(new_transformedspace_flambda_original(2650<wav_rest & wav_rest<2950))) max(real(new_transformedspace_flambda_original(2650<wav_rest & wav_rest<2950)))])
            
            %            thisistitle2=['zoomedin',specname,spectrafrom,spectrauntil]
            cd /home/ps738/7vardata/data/
            saveas(zoomin,char(strcat('zoomedin_',specnamefortitle,'_',num2str(spectrafrom),'_',num2str(spectrauntil))),'epsc')
            cd /home/ps738/7vardata/matlab/
            %            location2=strcat('/home/ps738/6vardata/data/errorcalc2/x12',thisistitle2)
            close all
            %% WRITE OBTAINED DATA TO OUTPUT FILE
            
            cd /home/ps738/7vardata/data/
            fprintf(fileID,'%21s %f %f %f %f %f %f\n',specname,spectrafrom,spectrauntil,ew_gauss1,new_ew_gauss1,lorentzianIntegralNT,meanmjd);
            cd /home/ps738/7vardata/matlab/
            
        end
    else
        continue
    end    
end


cd /home/ps738/7vardata/data/
fclose(fileID);
cd /home/ps738/7vardata/matlab/


function [skipthisone,redshift]=datainpeakrange(lower_window,higher_window,united,qsolist,k)
[specname,wav_rest,flambda,errortot,redshift]=readin(united,qsolist,k);
no_data_in_MgII_range=isempty(find(lower_window(2) <= wav_rest & wav_rest <= higher_window(1)))
if no_data_in_MgII_range==1
    skipthisone=0
end
if no_data_in_MgII_range==0
    skipthisone=1
end
end
function [specname,wav_rest,flambda,errortot,redshift]=readin(qsolistALLdat,qsolist,k)
ozdes_path='/data/desardata4/OzDES/Spectra/170130/';
objid_cell=qsolist(k,1);
objid1=objid_cell{1};
objid = strrep(objid1,'#','')

objid_cell
%%

specname=strcat(objid,'.fits');

disp('lehet hogy engem is megtalalsz')



redshift1=qsolistALLdat.data(:,1);
redshift=redshift1(k);
info=fitsinfo(strcat(ozdes_path,specname));

spectot=fitsread(strcat(ozdes_path,specname),'Primary');
errortot=fitsread(strcat(ozdes_path,specname),'Image',1);
npix=length(spectot);
pix=[1:npix];

iwav1=find(ismember(info.PrimaryData.Keywords(:,1),'CRVAL1'));
wav1=cell2mat(info.PrimaryData.Keywords(iwav1,2));
ipix1=find(ismember(info.PrimaryData.Keywords(:,1),'CRPIX1'));
pix1=cell2mat(info.PrimaryData.Keywords(ipix1,2));
idelta=find(ismember(info.PrimaryData.Keywords(:,1),'CDELT1'));
deltaw=cell2mat(info.PrimaryData.Keywords(idelta,2));

wav=wav1-(pix1-pix)*deltaw;
wav_rest=wav./(1+redshift);
flambda=spectot;
end
function [continuum, transformedspace_flambda, zeros1]=totransformedspace(lower_window,higher_window,wav_rest,flambda)
lambda_rest=wav_rest;
lwspeclam_index1=find(lower_window(1) <= lambda_rest & lambda_rest <= lower_window(2));
lwspeclam=lambda_rest(lwspeclam_index1);
lwflambda =flambda(lwspeclam_index1);

higherspeclam_index1=find(higher_window(1) <= lambda_rest & lambda_rest <= higher_window(2));
higherspeclam=lambda_rest(higherspeclam_index1);
higherflambda =flambda(higherspeclam_index1);

log_lambda_rest=log(lambda_rest);
log_flambda=log(flambda');

log_lw_window=log(lower_window);
log_higher_window=log(higher_window);

log_lw_window(1);
log_lw_window(2);
log_lwspeclam_index1=find(log_lw_window(1) <= log_lambda_rest & log_lambda_rest <= log_lw_window(2));
log_lwspeclam=log_lambda_rest(log_lwspeclam_index1);
log_lwflambda =log_flambda(log_lwspeclam_index1);

log_higherspeclam_index1=find(log_higher_window(1) <= log_lambda_rest & log_lambda_rest <= log_higher_window(2));
log_higherspeclam=log_lambda_rest(log_higherspeclam_index1);
log_higherflambda =log_flambda(log_higherspeclam_index1);

log_lw_avg_lambda_rest=nanmean(log_lwspeclam);
log_lw_avg_flambda=nanmean(log_lwflambda);

log_higher_avg_lambda_rest=nanmean(log_higherspeclam);
log_higher_avg_flambda=nanmean(log_higherflambda);

log_x = [log_lw_avg_lambda_rest log_higher_avg_lambda_rest];
log_y = [log_lw_avg_flambda log_higher_avg_flambda];
log_p = polyfit(log_x,log_y,1);

log_yy= log_p(1).*log_x+log_p(2);

yyy = exp(log_p(2)).*lambda_rest.^(log_p(1));

divided_spec=(flambda./yyy);
divided_lwflambda=(lwflambda./yyy(lwspeclam_index1));
divided_higherflambda=(higherflambda./yyy(higherspeclam_index1));

ones1=ones(size(lambda_rest));

int_divided_spec = ones1-divided_spec;
int_divided_lwflambda=1-divided_lwflambda;
int_divided_higherflambda=1-divided_higherflambda;

zeros1=zeros(size(lambda_rest));

MgIIspeclam_index1=find(lower_window(2) <= lambda_rest & lambda_rest <= higher_window(1));
MgIIspeclam=lambda_rest(MgIIspeclam_index1);
MgIIflambda =flambda(MgIIspeclam_index1);
divided_MgIIflambda=MgIIflambda./yyy(MgIIspeclam_index1);
int_divided_MgIIflambda=1-divided_MgIIflambda;

continuum = yyy;
transformedspace_flambda=int_divided_spec;

end

function [varargout]=equivalentwidth(lower_window, higher_window, wav_rest,transformedspace_flambda,wewantplot)
MgIIspeclam_index1=find(lower_window(2) < wav_rest & wav_rest < higher_window(1));
wav_rest_MgII=wav_rest(MgIIspeclam_index1)';
transformedspace_flambda_MgII=transformedspace_flambda(MgIIspeclam_index1);

[nanfiltered_wav_rest_MgII,nanfiltered_transformedspace_flambda_MgII]=nanfilterLocal(wav_rest_MgII,transformedspace_flambda_MgII);

if wewantplot==1
    [f1,gof] = fit(nanfiltered_wav_rest_MgII,nanfiltered_transformedspace_flambda_MgII','gauss1');
    x=[lower_window(2):higher_window(1)];
    coeffvals1 = coeffvalues(f1);
    a1=coeffvals1(1);
    b1=coeffvals1(2);
    c1=coeffvals1(3);
    varargout={a1,b1,c1,x,gof};
end
if wewantplot==2
    [f1,gof] = fit(nanfiltered_wav_rest_MgII,nanfiltered_transformedspace_flambda_MgII','gauss1');
    xx=[lower_window(2):higher_window(1)];
    ggooff=gof
    coeffvals1 = coeffvalues(f1);
    aa1=coeffvals1(1);
    bb1=coeffvals1(2);
    cc1=coeffvals1(3);
    varargout={aa1,bb1,cc1,xx,ggooff};
end
if wewantplot==3
    e=[nanfiltered_wav_rest_MgII,real(nanfiltered_transformedspace_flambda_MgII')]
    [lorentzianIntegralNT,params2] = INTEGloretz(e)
    varargout = {nanfiltered_wav_rest_MgII,lorentzianIntegralNT,params2}
end
    function [lorentzianIntegral,params2]=INTEGloretz(e_in)
        
        [yprime2 params2 resnorm2 residual2] = lorentzfit(e_in(1:end,1),e_in(1:end,2),[1 1 1],[-Inf 2700 -Inf; Inf 2900 Inf],'3');
        
        func_lorentz = @(x) (2*params2(1)./pi)*(params2(3)./(4*(x-params2(2)).^2 + params2(3).^2))
        lorentzianIntegral = integral(func_lorentz,2700,2900)
        function varargout = lorentzfit(x,y,varargin)
            % LORENTZFIT fits a single- or multi-parameter Lorentzian function to data
            %
            %   LORENTZFIT(X,Y) returns YPRIME(X), a Lorentzian fit to the data
            %   found using LSQCURVEFIT. The function Y(X) is fit by the model:
            %       YPRIME(X) = P1./((X - P2).^2 + P3) + C.
            %
            %   [YPRIME PARAMS RESNORM RESIDUAL] = LORENTZFIT(X,Y) returns YPRIME(X)
            %   values in addition to fit-parameters PARAMS = [P1 P2 P3 C]. The RESNORM
            %   and RESIDUAL outputs from LSQCURVEFIT are also returned.
            %
            %   [...] = LORENTZFIT(X,Y,P0) can be used to provide starting
            %   values (P0 = [P01 P02 P03 C0]) for the parameters in PARAMS.
            %
            %   [...] = LORENTZFIT(X,Y,P0,BOUNDS) may be used to define lower
            %   and upper bounds for the possbile values for each parameter in PARAMS.
            %       BOUNDS = [LB1 LB2 LB3 LB4;
            %                 UB1 UB2 UB3 UB4].
            %   If the user does not wish to manually define values for P0, it may be
            %   enetered as an empty matrix P0 = []. In this case, default values will
            %   be used. The default bounds for all parameters are (-Inf,Inf).
            %
            %   [...] = LORENTZFIT(X,Y,P0,BOUNDS,NPARAMS) may be used to specify the
            %   number of parameters used in the Lorentzian fitting function. The
            %   number of parameters defined in P0 and BOUNDS must match the function
            %   specified by NPARAMS. If the user does not wish to manually define
            %   values for P0 or BOUNDS, both may be enetered as empty matricies:
            %   P0 = []; BOUNDS = [].
            %
            %   -NPARAMS options
            %
            %           '1'     - Single parameter Lorentzian (no constant term)
            %                     L1(X) = 1./(P1(X.^2 + 1))
            %
            %           '1c'    - Single parameter Lorentzian (with constant term)
            %                     L1C(X) = 1./(P1(X.^2 + 1)) + C
            %
            %           '2'     - Two parameter Lorentzian (no constant term)
            %                     L2(X) = P1./(X.^2 + P2)
            %
            %           '2c'    - Two parameter Lorentzian (with constant term)
            %                     L2C(X) = P1./(X.^2 + P2) + C
            %
            %           '3'     - Three parameter Lorentzian (no constant term)
            %                     L3(X) = P1./((X - P2).^2 + P3)
            %
            % [DEFAULT] '3c'    - Three parameter Lorentzian (with constant term)
            %                     L3C(X) = P1./((X - P2).^2 + P3) + C
            %
            %   [...] = LORENTZFIT(X,Y,P0,BOUNDS,NPARAMS,OPTIONS) defines the OPTIONS
            %   array for the MATLAB function LSQCURVEFIT. OPTIONS can be set using the
            %   following command:
            %
            %       OPTIONS = optimset('PARAM1',VALUE1,'PARAM2',VALUE2,...);
            %
            %   See the help documentation for OPTIMSET for more details.
            %
            %
            %   X and Y must be the same size, numeric, and non-complex. P0 and BOUNDS
            %   must also be numeric and non-complex. NPARAMS is a character array.
            %
            %   Examples:
            %       x = -16:0.1:35;
            %       y = 19.4./((x - 7).^2 + 15.8) + randn(size(x))./10;
            %       [yprime1 params1 resnorm1 residual1] = lorentzfit(x,y,[20 10 15 0]);
            %       figure; plot(x,y,'b.','LineWidth',2)
            %       hold on; plot(x,yprime1,'r-','LineWidth',2)
            %
            %       [yprime2 params2 resnorm2 residual2] = lorentzfit(x,y,[],[],'3');
            %       figure; plot(x,y,'b.','LineWidth',2)
            %       hold on; plot(x,yprime2,'r-','LineWidth',2)
            %
            %   See also: lsqcurvefit.
            
            % Jered R Wells
            % 11/15/11
            % jered [dot] wells [at] duke [dot] edu
            %
            % v1.6 (2015/07/31)
            %
            % REF: http://www.home.uos.de/kbetzler/notes/fitp.pdf
            %
            %   UPDATES
            %   v1.5 - 2015/07/22 - jrw
            %       Added INPUTCHECK and OPTIONS. Rearranged SWITCH loop.
            %   v1.6 - 2015/07/31 - jrw
            %       Added stopping thresholds for fitting based on magnitude of input data.
            
            % CHECK INPUTS
            narginchk(2,6);
            nargoutchk(0,4);
            fname = 'lorentzfit';
            
            % Checked required inputs
            inputcheck(x,{'numeric'},{'real','nonnan','nonempty','finite'},fname,'X',1);
            inputcheck(y,{'numeric'},{'real','nonnan','nonempty','finite','size',size(x)},fname,'Y',2);
            
            % Set defaults for optional inputs
            p3 = ((max(x(:))-min(x(:)))./10).^2;
            p2 = (max(x(:))+min(x(:)))./2;
            p1 = max(y(:)).*p3;
            c = min(y(:));
            
            optargs = {[],[],'3c',optimset('TolFun',max(mean(y(:))*1e-6,1e-15),'TolX',max(mean(x(:))*1e-6,1e-15))};
            
            numvarargs = length(varargin);
            for ii = 1:numvarargs; if isempty(varargin{ii}); varargin{ii} = optargs{ii}; end; end
            
            % Now put these defaults into the valuesToUse cell array,
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            % or ...
            % [optargs{1:numvarargs}] = varargin{:};
            
            % Place optional args in memorable variable names
            [p0,bounds,nparams,options] = optargs{:};
            
            % Check optional inputs
            if ~isempty(p0)
                inputcheck(p0,{'numeric'},{'real','nonnan','vector'},fname,'P0',3);
            end
            if ~isempty(bounds)
                inputcheck(bounds,{'numeric'},{'real','nonnan','nrows',2},fname,'BOUNDS',4);
                lb = bounds(1,:);
                ub = bounds(2,:);
            else
                lb = [];
                ub = [];
            end
            inputcheck(nparams,{'char'},{},fname,'NPARAMS',5);
            inputcheck(options,{'struct'},{},fname,'OPTIONS',6);
            
            % PROCESS
            switch lower(nparams)
                case '1'
                    % Define P0, LB, UB
                    if isempty(p0);
                        p1 = max(y(:)); p0 = p1;
                    elseif numel(p0)~=1
                        error 'P0 must be empty or have one element for NPARAMS = ''1''';
                    end
                    
                    if isempty(bounds)
                        lb = -Inf; ub = Inf;
                    elseif ~all(size(bounds)==[2 1])
                        error 'BOUNDS must be empty or it must be a 2x1 matrix for NPARAMS = ''1''';
                    else
                        lb = bounds(1,:); ub = bounds(2,:);
                    end
                    
                    if any(lb>=ub)
                        error 'Lower bounds must be less than upper bounds';
                    end
                    
                    [params,resnorm,residual] = lsqcurvefit(@lfun1,p0,x,y,lb,ub,options);
                    yprime = lfun1(params,x);
                    
                case '1c'
                    % Define P0, LB, UB
                    if isempty(p0);
                        p1 = max(y(:)); p0 = [p1 c];
                    elseif numel(p0)~=2
                        error 'P0 must be empty or have two elements for NPARAMS = ''1c''';
                    end
                    
                    if isempty(bounds)
                        lb = [-Inf,-Inf]; ub = [Inf,Inf];
                    elseif ~all(size(bounds)==[2 2])
                        error 'BOUNDS must be empty or it must be a 2x2 matrix for NPARAMS = ''1c''';
                    else
                        lb = bounds(1,:); ub = bounds(2,:);
                    end
                    
                    if any(lb>=ub)
                        error 'Lower bounds must be less than upper bounds';
                    end
                    
                    [params,resnorm,residual] = lsqcurvefit(@lfun1c,p0,x,y,lb,ub,options);
                    yprime = lfun1c(params,x);
                case '2'
                    % Define P0, LB, UB
                    if isempty(p0);
                        p2 = ((max(x(:))-min(x(:)))./10).^2;
                        p1 = max(y(:)).*p2;
                        p0 = [p1 p2];
                    elseif numel(p0)~=2
                        error 'P0 must be empty or have two elements for NPARAMS = ''2''';
                    end
                    
                    if isempty(bounds)
                        lb = [-Inf,-Inf]; ub = [Inf,Inf];
                    elseif ~all(size(bounds)==[2 2])
                        error 'BOUNDS must be empty or it must be a 2x2 matrix for NPARAMS = ''2''';
                    else
                        lb = bounds(1,:); ub = bounds(2,:);
                    end
                    
                    if any(lb>=ub)
                        error 'Lower bounds must be less than upper bounds';
                    end
                    
                    [params,resnorm,residual] = lsqcurvefit(@lfun2,p0,x,y,lb,ub,options);
                    yprime = lfun2(params,x);
                case '2c'
                    % Define P0, LB, UB
                    if isempty(p0);
                        p2 = ((max(x(:))-min(x(:)))./10).^2;
                        p1 = max(y(:)).*p2;
                        p0 = [p1 p2 c];
                    elseif numel(p0)~=3
                        error 'P0 must be empty or have three elements for NPARAMS = ''2c''';
                    end
                    
                    if isempty(bounds)
                        lb = [-Inf,-Inf,-Inf]; ub = [Inf,Inf,Inf];
                    elseif ~all(size(bounds)==[2 3])
                        error 'BOUNDS must be empty or it must be a 2x3 matrix for NPARAMS = ''2c''';
                    else
                        lb = bounds(1,:); ub = bounds(2,:);
                    end
                    
                    if any(lb>=ub)
                        error 'Lower bounds must be less than upper bounds';
                    end
                    
                    [params,resnorm,residual] = lsqcurvefit(@lfun2c,p0,x,y,lb,ub,options);
                    yprime = lfun2c(params,x);
                case '3'
                    % Define P0, LB, UB
                    if isempty(p0);
                        p0 = [p1 p2 p3];
                    elseif numel(p0)~=3
                        error 'P0 must be empty or have three elements for NPARAMS = ''3''';
                    end
                    
                    if isempty(bounds)
                        lb = [-Inf,-Inf,-Inf]; ub = [Inf,Inf,Inf];
                    elseif ~all(size(bounds)==[2 3])
                        error 'BOUNDS must be empty or it must be a 2x3 matrix for NPARAMS = ''3''';
                    else
                        lb = bounds(1,:); ub = bounds(2,:);
                    end
                    
                    if any(lb>=ub)
                        error 'Lower bounds must be less than upper bounds';
                    end
                    
                    [params,resnorm,residual] = lsqcurvefit(@lfun3,p0,x,y,lb,ub,options);
                    yprime = lfun3(params,x);
                case '3c'
                    % Define P0, LB, UB
                    if isempty(p0);
                        p0 = [p1 p2 p3 c];
                    elseif numel(p0)~=4
                        error 'P0 must be empty or have four elements for NPARAMS = ''3c''';
                    end
                    
                    if isempty(bounds)
                    elseif ~all(size(bounds)==[2 4])
                        error 'BOUNDS must be empty or it must be a 2x4 matrix for NPARAMS = ''3c''';
                    else
                        lb = bounds(1,:); ub = bounds(2,:);
                    end
                    
                    if any(lb>=ub)
                        error 'Lower bounds must be less than upper bounds';
                    end
                    
                    [params,resnorm,residual] = lsqcurvefit(@lfun3c,p0,x,y,lb,ub,options);
                    yprime = lfun3c(params,x);
                otherwise
                    warning('MATLAB:lorentzfit:default','Fitting default (''3c'') Lorentzian model.')
                    [params,resnorm,residual] = lsqcurvefit(@lfun3c,p0,x,y,lb,ub,options);
                    yprime = lfun3c(params,x);
            end
            
            varargout = {yprime,params,resnorm,residual};
            
        end % MAIN
        function F = lfun3(p,x)
            %F = p(1)./((x-p(2)).^2+p(3));
            F = (2*p(1)./pi)*(p(3)./(4*(x-p(2)).^2 + p(3).^2))
        end % LFUN3
        function varargout = inputcheck(A,varargin)
            % INPUTCHECK checks the validity of input array with VALIDATEATTRIBUTES
            %
            % INPUTCHECK(A,CLASSES,ATTRIBUTES) validates that array A belongs
            % to at least one of the specified CLASSES and has all of the specified
            % ATTRIBUTES. If A does not meet the criteria, MATLAB issues a formatted
            % error message. INPUTCHECK(A,CLASSES,ATTRIBUTES) is equivalent to
            % VALIDATEATTRIBUTES(A,CLASSES,ATTRIBUTES) with the exception that custom
            % CLASSES and ATTRIBUTES can be easily implemented by editing the function.
            %
            % Some additional input error-checking is also provided. Also, empty arrays
            % may be defined if any of the input arguments is irrelevant.
            %
            % INPUTCHECK(A,CLASSES,ATTRIBUTES,ARGINDEX) includes the
            % position of the input in your function argument list as part of any
            % generated error messages.
            %
            % INPUTCHECK(A,CLASSES,ATTRIBUTES,FUNCNAME) includes the
            % specified function name in generated error identifiers.
            %
            % INPUTCHECK(A,CLASSES,ATTRIBUTES,FUNCNAME,VARNAME) includes the
            % specified variable name in generated error messages.
            %
            % INPUTCHECK(A,CLASSES,ATTRIBUTES,FUNCNAME,VARNAME,ARGINDEX)
            % includes the specified information in the generated error messages or
            % identifiers.
            %
            % [V,ME] = INPUTCHECK(A,CLASSES,ATTRIBUTES,FUNCNAME,VARNAME,ARGINDEX,VERBOSE)
            % toggles the error message generation state. VERBOSE = TRUE will produce
            % error messages which halt function execution and return control to the
            % command line if array A is not valid. Otherwise, the validity state V
            % will be output as logical TRUE and the error catch variable ME is empty.
            % VERBOSE = FALSE will output validity V as either TRUE or FALSE in
            % addition to reporting errors in ME without halting function execution
            % when A is not valid. By default, VERBOSE = TRUE.
            %
            %   Input Arguments:
            %
            %   A          Any class of array.
            %
            %   CLASSES    Cell array of strings that specify valid classes for array A.
            %              For example, if CLASSES = {'logical','cell'}, A must be a
            %              logical array or a cell array. The string 'numeric' is an
            %              abbreviation for the classes uint8, uint16, uint32, uint64,
            %              int8, int16, int32, int64, single, double. CLASSES may
            %              include MATLAB built-in or custom classes:
            %
            %             'numeric'         Any value for which the isnumeric function
            %                               returns true, including int8, int16, int32,
            %                               int64, uint8, uint16, uint32, uint64,
            %                               single, or double.
            %             'single'          Single-precision number.
            %             'double'          Double-precision number.
            %             'int8'            Signed 8-bit integer.
            %             'int16'           Signed 16-bit integer.
            %             'int32'           Signed 32-bit integer.
            %             'int64'           Signed 64-bit integer.
            %             'uint8'           Unsigned 8-bit integer.
            %             'uint16'          Unsigned 16-bit integer.
            %             'uint32'          Unsigned 32-bit integer.
            %             'uint64'          Unsigned 64-bit integer.
            %             'logical'         Logical true or false.
            %             'char'            Character or string.
            %             'struct'          Structure array.
            %             'cell'            Cell array.
            %             'function_handle'	Scalar function handle.
            %
            %             Examples:
            %                   % Define CLASSES using the following syntax
            %                   classes = {};
            %                   classes = {'double'};
            %                   classes = {'int8','int16'};
            %
            %   ATTRIBUTES Cell array that contains descriptions of valid attributes
            %              for array A. For example, if ATTRIBUTES = {'real','finite'}
            %              A must contain only real and finite values.
            %
            %              Supported attributes include:
            %
            %             Attributes that describe the size and shape of array A:
            %             '2d'              Two-dimensional array, including scalars,
            %                               vectors, matrices, and empty arrays.
            %             'column'          Column vector, N-by-1.
            %             'row'             Row vector, 1-by-N.
            %             'scalar'          Scalar value, 1-by-1.
            %             'vector'          Row or column vector, or a scalar value.
            %             'size', [d1,...,dN]	Array with dimensions d1-by-...-by-dN.
            %                                   If you do not want to check a
            %                                   particular dimension, specify NaN for
            %                                   that dimension, such as [3,4,NaN,2].
            %               'ndims', N      Array with N dimensions.
            %             'numel', N        Array with N elements.
            %               'numel>', N     Array with more than N elements.
            %               'numel<', N     Array with fewer than N elements.
            %             'ncols', N        Array with N columns.
            %               'ncols>', N     Array with more than N columns.
            %               'ncols<', N     Array with fewer than N columns.
            %             'nrows', N        Array with N rows.
            %               'nrows>', N     Array with more than N rows.
            %               'nrows<', N     Array with fewer than N rows.
            %             'nonempty'        No dimensions equal zero.
            %             'nonsparse'       Array that is not sparse.
            %
            %             Attributes that specify valid ranges for values in A:
            %             '>', N            All values are greater than N.
            %             '>=', N           All values are greater than or equal to N.
            %             '<', N            All values are less than N.
            %             '<=', N           All values are less than or equal to N.
            %
            %             Attributes that check types of values in A, where A is a
            %             numeric or logical array:
            %             'binary'          Array of ones and zeros.
            %             'even'            All elements are even integers
            %                               (includes zero).
            %             'odd'             All elements are odd integers.
            %             'integer'         All elements are integer-valued.
            %             'real'            All elements are real.
            %             'finite'          All elements are finite.
            %             'nonnan'          No elements equal to NaN (Not a Number).
            %             'nonnegative'     All elements are nonnegative.
            %             'nonzero'         All elements are nonzero.
            %             'positive'        All elements are positive.
            %
            %              Some attributes also require numeric values. For those
            %              attributes, the numeric value or vector must immediately
            %              follow the attribute name string. For example,
            %              {'>=', 5, '<=', 10, size, [3 4 2]} checks that all
            %              values of A are between 5 and 10, and that A is 3-by-4-by-2.
            %
            %             Examples:
            %                   % Define ATTRIBUTES using the following syntax
            %                   attributes = {};
            %                   attributes = {'2d'};
            %                   attributes = {'2d','>',5,'numel',1000};
            %
            %   ARGINDEX   Positive integer that specifies the position of the input
            %              argument.
            %
            %   FUNCNAME   String that specifies the function name. If you specify an
            %              empty string, '', FUNCNAME is ignored.
            %
            %   VARNAME    String that specifies input argument name. If you specify an
            %              empty string, '', VARNAME is ignored.
            %
            %   VERBOSE    Logical that toggles the error-generation state of
            %              INPUTCHECK.
            %
            %   Output Arguments:
            %
            %   V          (Logical) validity state of A.
            %
            %   ME         Class MException object for error capture.
            %
            %   Example: Create a three dimensional array and then check for the
            %            attribute '2d'.
            %
            %       A = [ 1 2 3; 4 5 6 ];
            %       B = [ 7 8 9; 10 11 12];
            %       C = cat(3,A,B);
            %       inputcheck(C,{'numeric'},{'2d'},'my_func','my_var',2)
            %
            %   This code throws an error and displays a formatted message:
            %
            %       Expected input number 2, my_var, to be two-dimensional.
            %
            % [REF: http://www.mathworks.com/help/techdoc/ref/validateattributes.html]
            %
            % See also: validateattributes.
            
            %
            % Jered R Wells
            % 08/15/2012
            %
            % v1.0 (08/17/2012)
            
            % Initialize variables
            V = true;
            ME = MException('','');
            %#ok<*SPERR>
            
            classes = {};
            attributes = {};
            funcname = '';
            varname = '';
            argindex = [];
            verbose = true;
            
            % Input check
            if nargin==1
                % This case is not very useful as no error checking actually occurs!
                warning('MATLAB:inputcheck:noCrit','No criteria for array validation')
                classes = {};
                attributes = {};
                validateattributes(A,classes,attributes);
            elseif nargin>1
                for ii = 1:length(varargin)
                    switch ii
                        case 1
                            classes = varargin{ii};
                            if isempty(classes); classes = {}; end
                            if ischar(classes); classes = {classes}; end
                            validateattributes(classes,{'cell'},{},'INPUTCHECK','CLASSES',2);
                        case 2
                            attributes = varargin{ii};
                            if isempty(attributes); attributes = {}; end
                            validateattributes(attributes,{'cell'},{},'INPUTCHECK','ATTRIBUTES',3);
                        case 3
                            funcname = varargin{ii};
                            if isempty(funcname); funcname = ''; end
                            if isnumeric(funcname)&&length(varargin)==3;
                                funcname = '';
                                argindex = varargin{ii};
                                if isempty(argindex)
                                    argindex = [];
                                else
                                    validateattributes(argindex,{'numeric'},{'positive','integer','numel',1},'INPUTCHECK','ARGINDEX',3);
                                end
                            end
                            validateattributes(funcname,{'char'},{},'INPUTCHECK','FUNCNAME',4);
                        case 4
                            varname = varargin{ii};
                            if isempty(varname); varname = ''; end
                            validateattributes(varname,{'char'},{},'INPUTCHECK','VARNAME',5);
                        case 5
                            argindex = varargin{ii};
                            if isempty(argindex)
                                argindex = [];
                            else
                                validateattributes(argindex,{'numeric'},{'positive','integer','numel',1},'INPUTCHECK','ARGINDEX',6);
                            end
                        case 6
                            verbose = varargin{ii};
                            validateattributes(verbose,{'logical'},{'numel',1},'INPUTCHECK','VERBOSE',7);
                        otherwise
                            error 'Too many input arguments'
                    end % SWITCH
                end % ii
            else
                error 'Array A must be defined'
            end % IF
            
            % Format strings
            if isempty(funcname); funcnameIC = '';
            else funcnameIC = ['Error using ',funcname];
            end
            if isempty(argindex); argindexIC = '';
            else argindexIC = ['number ',num2str(argindex),',',' '];
            end
            if isempty(varname); varnameIC = '';
            else varnameIC = [varname,',',' '];
            end
            
            
            % Process
            % Check for custom CLASSES
            % N/A
            
            % Check for custom ATTRIBUTES
            ii = 1;
            customAttributes = true(length(attributes));
            while ii<=length(attributes)
                if ischar(attributes{ii})
                    tmpName = lower(attributes{ii});
                    switch tmpName
                        % Create a new case for each custom ATTRIBUTE
                        % Be sure to check for the requirement of N
                        case 'ndims'
                            % Check for N
                            ii = ii + 1;                % Advance ii if N is required
                            tmpN = attributes{ii};
                            validateattributes(tmpN,{'numeric'},{'integer','nonnegative','numel',1},'INPUTCHECK','ndims',3);
                            % Validate array A
                            if ndims(A)~=tmpN
                                % Formatted string
                                MSGID = 'MATLAB:inputcheck:incorrectNdims';
                                ERRMSG = sprintf('%s\nExpected input %s%sto be an array with number of dimensions equal to %i.',funcnameIC,argindexIC,varnameIC,tmpN);
                                if verbose
                                    error(ERRMSG)
                                else
                                    ME = MException(MSGID, ERRMSG);
                                    V = false;
                                    break
                                end
                            end
                        case 'numel>'
                            % Check for N
                            ii = ii + 1;                % Advance ii if N is required
                            tmpN = attributes{ii};
                            validateattributes(tmpN,{'numeric'},{'integer','nonnegative','numel',1},'INPUTCHECK','numel>',3);
                            % Validate array A
                            if ~(numel(A)>tmpN)
                                % Formatted string
                                MSGID = 'MATLAB:inputcheck:incorrectNumel';
                                ERRMSG = sprintf('%s\nExpected input %s%sto be an array with number of elements greater than %i.',funcnameIC,argindexIC,varnameIC,tmpN);
                                if verbose
                                    error(ERRMSG)
                                else
                                    ME = MException(MSGID, ERRMSG);
                                    V = false;
                                    break
                                end
                            end
                        case 'numel<'
                            % Check for N
                            ii = ii + 1;                % Advance ii if N is required
                            tmpN = attributes{ii};
                            validateattributes(tmpN,{'numeric'},{'integer','nonnegative','numel',1},'INPUTCHECK','numel<',3);
                            % Validate array A
                            if ~(numel(A)<tmpN)
                                % Formatted string
                                MSGID = 'MATLAB:inputcheck:incorrectNumel';
                                ERRMSG = sprintf('%s\nExpected input %s%sto be an array with number of elements less than %i.',funcnameIC,argindexIC,varnameIC,tmpN);
                                if verbose
                                    error(ERRMSG)
                                else
                                    ME = MException(MSGID, ERRMSG);
                                    V = false;
                                    break
                                end
                            end
                        case 'ncols>'
                            % Check for N
                            ii = ii + 1;                % Advance ii if N is required
                            tmpN = attributes{ii};
                            validateattributes(tmpN,{'numeric'},{'integer','nonnegative','numel',1},'INPUTCHECK','ncols>',3);
                            % Validate array A
                            if ~(size(A,2)>tmpN)
                                % Formatted string
                                MSGID = 'MATLAB:inputcheck:incorrectNcols';
                                ERRMSG = sprintf('%s\nExpected input %s%sto be an array with number of columns greater than %i.',funcnameIC,argindexIC,varnameIC,tmpN);
                                if verbose
                                    error(ERRMSG)
                                else
                                    ME = MException(MSGID, ERRMSG);
                                    V = false;
                                    break
                                end
                            end
                        case 'ncols<'
                            % Check for N
                            ii = ii + 1;                % Advance ii if N is required
                            tmpN = attributes{ii};
                            validateattributes(tmpN,{'numeric'},{'integer','nonnegative','numel',1},'INPUTCHECK','ncols<',3);
                            % Validate array A
                            if ~(size(A,2)<tmpN)
                                % Formatted string
                                MSGID = 'MATLAB:inputcheck:incorrectNcols';
                                ERRMSG = sprintf('%s\nExpected input %s%sto be an array with number of columns less than %i.',funcnameIC,argindexIC,varnameIC,tmpN);
                                if verbose
                                    error(ERRMSG)
                                else
                                    ME = MException(MSGID, ERRMSG);
                                    V = false;
                                    break
                                end
                            end
                        case 'nrows>'
                            % Check for N
                            ii = ii + 1;                % Advance ii if N is required
                            tmpN = attributes{ii};
                            validateattributes(tmpN,{'numeric'},{'integer','nonnegative','numel',1},'INPUTCHECK','nrows>',3);
                            % Validate array A
                            if ~(size(A,1)>tmpN)
                                % Formatted string
                                MSGID = 'MATLAB:inputcheck:incorrectNrows';
                                ERRMSG = sprintf('%s\nExpected input %s%sto be an array with number of rows greater than %i.',funcnameIC,argindexIC,varnameIC,tmpN);
                                if verbose
                                    error(ERRMSG)
                                else
                                    ME = MException(MSGID, ERRMSG);
                                    V = false;
                                    break
                                end
                            end
                        case 'nrows<'
                            % Check for N
                            ii = ii + 1;                % Advance ii if N is required
                            tmpN = attributes{ii};
                            validateattributes(tmpN,{'numeric'},{'integer','nonnegative','numel',1},'INPUTCHECK','nrows<',3);
                            % Validate array A
                            if ~(size(A,1)<tmpN)
                                % Formatted string
                                MSGID = 'MATLAB:inputcheck:incorrectNrows';
                                ERRMSG = sprintf('%s\nExpected input %s%sto be an array with number of rows less than %i.',funcnameIC,argindexIC,varnameIC,tmpN);
                                if verbose
                                    error(ERRMSG)
                                else
                                    ME = MException(MSGID, ERRMSG);
                                    V = false;
                                    break
                                end
                            end
                        otherwise
                            % If TMPNAME is not a custom attribute, assign FALSE to the
                            % appropriate element in CUSTOMATTRIBUTES
                            customAttributes(ii) = false;
                    end % SWITCH
                else
                    % If ATTRIBUTES{ii} is not a custom attribute, assign FALSE to the
                    % appropriate element in CUSTOMATTRIBUTES
                    customAttributes(ii) = false;
                end % IF
                ii = ii + 1;
            end % ii
            attributes = attributes(~customAttributes);
            
            if verbose
                validateattributes(A,classes,attributes,funcname,varname,argindex)
            else
                try
                    validateattributes(A,classes,attributes,funcname,varname,argindex)
                catch ME
                    V = false;
                end
            end % IF
            
            varargout = {V,ME};
            
        end % MAIN
    end
end
function [nanfiltered_wav_rest_MgII,nanfiltered_transformedspace_flambda_MgII]=nanfilterLocal(wav_rest_MgII,transformedspace_flambda_MgII)

Bnanfiltered_wav_rest_MgII=wav_rest_MgII(~isnan(wav_rest_MgII));
Bnanfiltered_transformedspace_flambda_MgII=transformedspace_flambda_MgII(~isnan(wav_rest_MgII));

nanfiltered_wav_rest_MgII=Bnanfiltered_wav_rest_MgII(~isnan(Bnanfiltered_transformedspace_flambda_MgII));
nanfiltered_transformedspace_flambda_MgII=Bnanfiltered_transformedspace_flambda_MgII(~isnan(Bnanfiltered_transformedspace_flambda_MgII));
end
function [new_lower_window,new_higher_window]=decideonwindows(united,lower_window_default,higher_window_default,k)
if united.data(k,29) == 0
    new_lower_window=[united.data(k,31) united.data(k,32)];
else
    new_lower_window=lower_window_default;
end
if united.data(k,30) == 0
    new_higher_window=[flagdocs.data(k,33) flagdocs.data(k,34)];
else
    new_higher_window=higher_window_default;
end
end
