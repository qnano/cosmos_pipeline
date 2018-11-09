classdef guiDarkCorrector < handle
   properties 
    Figure
    Axes
    CLimit
    
    status
    handels
    
    gadObj
    gdrObj  
    fitcfd = false
   end
   
   properties (Access = private)
   end
   events
       DataChange
   end
   methods       
       function whoAmI(obj,src,~)
            basevars = evalin('base','whos');
            testClassvars = basevars(strcmp({basevars.class},class(obj)));

            found = false;
            for i = 1:length(testClassvars)
                if(eq(evalin('base',testClassvars(i).name),obj))
                    found = true;
                    obj.name =testClassvars(i).name;
                end
            end
       end
       
     
       
    function [k_s_est,k_s_est_conf,...
            k_ns_1_est,k_ns_1_est_conf,...
            k_ns_2_est,k_ns_2_est_conf,...
            h_1_est,h_1_est_conf,...
            h_2_est,h_2_est_conf,h_est,h_est_conf] = getCorrectedOffrateFun(obj,initialValueMLE)
      
        dwellTimes = obj.gdrObj.getDwellTimes;

        if nargin < 2 || isempty(initialValueMLE)
            k_ns_intial = min(1,max(1./mean(dwellTimes),1e-6));
            h_initial = 0.5;
        else
            k_ns_intial = initialValueMLE(1);
            h_initial = initialValueMLE(2);
        end
        x = sym('x');
        k_ns_1 = sym('k_ns_1');
        k_ns_2 = sym('k_ns_2');
        h = sym('h');

        cdffun = 1-(1-h)*exp(-k_ns_1*x)-h*exp(-k_ns_2*x);
        if obj.fitcfd
            cdffun =  matlabFunction(cdffun);
            [f,d] = ecdf(dwellTimes);
            mu = sym('mu');

            obj_func = @(mu) norm(cdffun(mu(1),mu(2),mu(3),d) - f);
            LB = [0; eps*10; eps*10];
            UB = [1; 1; 1];
            mu = fmincon(obj_func, [0.5,1./mean(dwellTimes),1./mean(dwellTimes)],[],[],[],[],LB,UB);  
            theta_est = mu;
            conf = 0.*mu;
        else
            a = matlabFunction(diff(cdffun,x));
            b = @(x,h,k_ns_1,k_ns_2) max(a(max(min(h,1),0),k_ns_1,k_ns_2,x),1e-2);

            [theta_est,conf] = mle(dwellTimes,'pdf',b,'start',[h_initial 2*k_ns_intial k_ns_intial/2]);
        end
        h_est  = max(min(theta_est(1),1),0);
        h_est_conf = conf(:,1);
        k_ns_1_est  = theta_est(2);
        k_ns_1_est_conf  = conf(:,2);
        k_ns_2_est  = theta_est(3);
        k_ns_2_est_conf  = conf(:,3);     
        
        % It's a single non specific distribution 
        if h_est == 0 || any(any(isnan(conf)))
            k_ns_1_est = min(1,max(1./mean(dwellTimes),1e-6));
            dwellTimes = obj.gadObj.getDwellTimes;
            k_ns_2_est=0;
            k_ns_2_est_conf=0;
            h_est_conf=0;
            h_est=0;
            k_s = sym('k_s');
            h_1 = sym('h_1');
            cdffun = h_1*(1-exp(-k_ns_1*x))+(1-h_1)*(1-exp(-k_s*x));
            if obj.fitcfd
                cdffun =  matlabFunction(cdffun);
                [f,d] = ecdf(dwellTimes);
                mu = sym('mu');

                obj_func = @(mu) norm(cdffun(mu(1),mu(2),mu(3),d) - f);
                LB = [0; eps*10; eps*10];
                UB = [1; 1; 1];
                mu = fmincon(obj_func, [0.5 k_ns_1_est min(1,max(1./mean(dwellTimes),1e-6))],[],[],[],[],LB,UB);  
                theta_est = mu;
                theta_conf = 0.*mu;
            else
                a = matlabFunction(diff(cdffun,x));
                b = @(x,k_ns_1,k_s,h_1) max(a(min(max(h_1,0),1),k_ns_1,k_s,x),1e-2);
                params = statset('mlecustom');
                params.TolX = 1e-12;
                params.TolFun = 1e-12;
                params.Display='iter';

                [theta_est,conf] = mle(dwellTimes,'pdf',b,'start',[ k_ns_1_est min(1,max(1./mean(dwellTimes),1e-6)) 0.5],'options',params);
            end
            h_1_est = min(max(theta_est(3),0),1);
            h_1_est_conf =  conf(:,3);

            h_2_est = 0;
            h_2_est_conf = 0;

            k_ns_1_est = theta_est(1);
            k_ns_1_est_conf = conf(:,1);

            k_s_est = theta_est(2);
            k_s_est_conf = conf(:,2);
            if h_1_est == 0 || any(any(isnan(conf)))
                h_1_est=0;
                k_ns_1_est=0;
                k_ns_1_est_conf=0;
            end
        else     
            dwellTimes = obj.gadObj.getDwellTimes;
            x = sym('x');
            k_ns_1 = sym('k_ns_1');
            k_ns_2 = sym('k_ns_2');
            k_s = sym('k_s');
            h_1 = sym('h_1');
            h_2 = sym('h_2');

            cdffun = h_1*(1-exp(-k_ns_1*x))+h_2*(1-exp(-k_ns_2*x))+(1-h_1-h_2)*(1-exp(-k_s*x));
            if obj.fitcfd
                cdffun =  matlabFunction(cdffun);
                [f,d] = ecdf(dwellTimes);
                mu = sym('mu');

                obj_func = @(mu) norm(cdffun(mu(1),mu(2),k_ns_1_est,k_ns_2_est,mu(3),d) - f);
                LB = [0; eps*10; eps*10];
                UB = [1; 1; 1];
                mu = fmincon(obj_func, [0.5 k_ns_1_est min(1,max(1./mean(dwellTimes),1e-6))],[],[],[],[],LB,UB);  
                theta_est = mu;
                theta_conf = 0.*mu;
            else            
                a = matlabFunction(diff(cdffun,x));

                b = @(x,h_1,h_2,k_s) max(a(max(min(h_1,1),0),max(min(h_2,1),0),k_ns_1_est,k_ns_2_est,k_s,x),1e-2);

                if nargin < 2 || isempty(initialValueMLE)
                    k_s_initial =  min(1,max(1./mean(dwellTimes),1e-6))
                else
                    h_initial_1 = initialValueMLE(3);
                    h_initial_2 = initialValueMLE(4);
                    k_s_initial = initialValueMLE(5);
                end

                [theta_est,conf] = mle(dwellTimes,'pdf',b,'start',[0.5 0.5 k_s_initial]);
                if any(isnan(conf))
                    warningMessage=true;
                end
            end
            h_1_est = max(min(theta_est(1),1),0);
            h_1_est_conf = conf(:,1);

            h_2_est =max(min(theta_est(2),1),0);
            h_2_est_conf = conf(:,2);

            k_s_est = theta_est(3);
            k_s_est_conf = conf(:,3);
        end
    end
       
	function obj = guiDarkCorrector(h,gadObj,gdrObj)
        obj.status=0;
        obj.gadObj = gadObj;
        obj.gdrObj = gdrObj;

        htab1 = uitab(h, 'Title', 'NS-Corrector');
        obj.Figure = htab1;
        htabgroup = uitabgroup(htab1);
        obj.handels.onRateTab = uitab(htabgroup, 'Title', 'On-Rate');
        obj.handels.offRateTab = uitab(htabgroup, 'Title', 'Off-Rate');
            
        obj.handels.parentMenu = uimenu('Label','Dark Corrector');
        obj.handels.childMenu(1) = uimenu('Parent',obj.handels.parentMenu,'Label','Dark Corrector Settings','Separator','on','Callback',@obj.setDarkCorrectorSettings)
        uimenu(obj.handels.parentMenu,'Label','Estimate On-Rate','Callback',@obj.getCorrectedOnrate);
        uimenu(obj.handels.parentMenu,'Label','Estimate Off-Rate','Callback',@obj.getCorrectedOffrate);
    end
    
    function obj = setDarkCorrectorSettings(obj,src,~)
            Title = 'Rastergram Settings';

            %%%% SETTING DIALOG OPTIONS
            Options.Resize = 'on';
            Options.Interpreter = 'tex';
            Options.CancelButton = 'on';
            Options.ApplyButton = 'off';
            
            Option.Dim = 1; % Horizontal dimension in fields
            
            Prompt = {};
            Formats = {};
                      
            Prompt(1,:) = {'Fit type','List',[]};
            Formats(1,1).type = 'list';
            Formats(1,1).style = 'popupmenu';
            Formats(1,1).items = {'Probablity density function (Pdf)','Cumulative density function (Cdf)'};
            DefAns.List =obj.fitcfd+1;
            [result,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);

            if ~Cancelled
                obj.fitcfd = (result.List == 2);
            end                
    end
    
    
    function [p21DrVec, p21AdVec,pDTDrVec, pDTAdVec, pDwellTime,pTime21st]  = testDwellTimes(obj,src,~)     

        idx = [0 cumsum(obj.gdrObj.AddedTrackingInfo.spots)];
        for dataSetIndex=1:length(obj.gdrObj.AddedTrackingInfo.spots)
            spots(dataSetIndex) = sum(obj.gdrObj.spotsIncluded(1+idx(dataSetIndex):idx(dataSetIndex+1)));
        end
        idxspots = [0 spots];

        for idxDataSet = 1:length(obj.gdrObj.AddedTrackingInfo.spots)
            time21stDr{idxDataSet} = obj.gdrObj.getTimes2FirstEvent(true,1+idxspots(dataSetIndex):idxspots(dataSetIndex+1));
            dwellTimesDr{idxDataSet} =  obj.gdrObj.getDwellTimes(false,1+idxspots(dataSetIndex):idxspots(dataSetIndex+1));
        end
        
        idx = [0 cumsum(obj.gadObj.AddedTrackingInfo.spots)];
        for dataSetIndex=1:length(obj.gadObj.AddedTrackingInfo.spots)
            spots(dataSetIndex) = sum(obj.gadObj.spotsIncluded(1+idx(dataSetIndex):idx(dataSetIndex+1)));
        end
        idxspots = [0 spots];

        
        for idxDataSet = 1:length(obj.gadObj.AddedTrackingInfo.spots)
            time21stAd{idxDataSet} = obj.gadObj.getTimes2FirstEvent(true,1+idxspots(dataSetIndex):idxspots(dataSetIndex+1));
            dwellTimesAd{idxDataSet} = obj.gadObj.getDwellTimes(false,1+idxspots(dataSetIndex):idxspots(dataSetIndex+1));
        end
        
        if length(time21stDr) > 1
            p21Dr = zeros(length(time21stDr),length(time21stDr));
            
            for i=1:length(time21stDr)
                for j=i+1:length(time21stDr)
                    p21Dr(i,j) = ranksum(time21stDr{i},time21stDr{j});
                end
            end
            p21DrVec = nonzeros(triu(p21Dr)');
            if length(p21DrVec) > 1
                p21DrVec = fdr_bh(p21DrVec);
            end
        end
        
        if length(dwellTimesDr) > 1
            pDTDr = zeros(length(dwellTimesDr),length(dwellTimesDr));
            
            for i=1:length(dwellTimesDr)
                for j=i+1:length(dwellTimesDr)
                    pDTDr(i,j) = ranksum(dwellTimesDr{i},dwellTimesDr{j});
                end
            end
            
            pDTDrVec = nonzeros(triu(pDTDr)');
            if length(pDTDrVec) > 1
                pDTDrVec = fdr_bh(pDTDrVec);
            end
        end

        if length(time21stAd) > 1
            p21Ad = zeros(length(time21stAd),length(time21stAd));
            for i=1:length(time21stAd)
                for j=i+1:length(time21stAd)
                    p21Ad(i,j) = ranksum(time21stAd{i},time21stAd{j});
                    
                end
            end
            p21AdVec = nonzeros(triu(p21Ad)');
            if length(p21AdVec) > 1
                p21AdVec = fdr_bh(p21AdVec);
            end    
        end
        
        if length(dwellTimesAd) > 1
            pDTAd = zeros(length(dwellTimesAd),length(dwellTimesAd));
            for i=1:length(dwellTimesAd)
                for j=i+1:length(dwellTimesAd)
                    pDTAd(i,j) = ranksum(dwellTimesAd{i},dwellTimesAd{j});     
                end
            end
            pDTAdVec = nonzeros(triu(pDTAd)');
            if length(pDTAdVec) > 1
                pDTAdVec = fdr_bh(pDTAdVec);
            end    
        end
        
        allTime21stAd = [];
        allDwellTimesAd = [];
        for idxDataSet = 1:length(obj.gadObj.AddedTrackingInfo.spots)
            allTime21stAd = cat(1,allTime21stAd,time21stAd{idxDataSet});
            allDwellTimesAd= cat(1,allDwellTimesAd,dwellTimesAd{idxDataSet});
        end     
        
        allTime21stDr = [];
        allDwellTimesDr = [];
        for idxDataSet = 1:length(obj.gdrObj.AddedTrackingInfo.spots)
            allTime21stDr = cat(1,allTime21stDr,time21stDr{idxDataSet});
            allDwellTimesDr= cat(1,allDwellTimesDr,dwellTimesDr{idxDataSet});
        end     
        
        pTime21st = ranksum(allTime21stDr, allTime21stAd);
        pDwellTime =ranksum(allDwellTimesAd, allDwellTimesDr);
        
        if isa(src,'matlab.ui.container.Menu')
            msg{1} = 'ANALYSIS';
            msg{end+1} = sprintf('Firstevent    p-value');
            idx = 1;
            for i=1:length(time21stAd)
                for j=i+1:length(time21stAd)
                    msg{end+1} = sprintf('%d/%d              %0.2g',i,j,p21AdVec(idx));
                    idx=idx+1;
                end
            end
            msg{end+1} = sprintf('Dwelltime     p-value');
            idx = 1;
            for i=1:length(time21stAd)
                for j=i+1:length(time21stAd)
                    msg{end+1} = sprintf('%d/%d              %0.2g',i,j,pDTAdVec(idx));
                    idx=idx+1;
                end
            end
            msg{end+1} = [];
            msg{end+1} = 'DARK';
            msg{end+1} = sprintf('Firstevent    p-value');
            idx = 1;
            for i=1:length(time21stDr)
                for j=i+1:length(time21stDr)
                    msg{end+1} = sprintf('%d/%d              %0.2g',i,j,p21DrVec(idx));
                    idx=idx+1;
                end
            end
            msg{end+1} = sprintf('Dwelltime     p-value');
            idx = 1;
            for i=1:length(time21stDr)
                for j=i+1:length(time21stDr)
                    msg{end+1} = sprintf('%d/%d               %0.2g',i,j,pDTDrVec(idx));
                    idx=idx+1;
                end
            end

            msg{end+1} = [];
            msg{end+1} = sprintf('ALL             p-value');
            msg{end+1} = sprintf('Firstevent    %0.2g',pDwellTime);
            msg{end+1} = sprintf('Dwelltime    %0.2g',pTime21st);
            msg{end+1} = [];
            msg{end+1} = [];
            msg{end+1} = sprintf('STATS INFO');
            msg{end+1} = sprintf('Test = rank-sum');
            msg{end+1} = sprintf('P-values are corrected for multiple comparison');
            msgbox(msg);
        end

    end
    
    function getCorrectedOffrate(obj,src,~)
        
        if ischar(src)
            switch upper(src)
                case 'PDF'
                    obj.fitcfd=0;
                case 'CDF'
                    obj.fitcfd=1;
            end
        end
        
        warningMessage= false;
        
        if isfield(obj.handels,'histOffDarkPlot') && ~isempty(obj.handels.histOffDarkPlot)
            obj.deletetry(obj.handels.histOffDarkPlot)
        end
        if isfield(obj.handels,'histOffPlot') && ~isempty(obj.handels.histOffPlot)
            obj.deletetry(obj.handels.histOffPlot)
        end
        
        [k_s_est,k_s_est_conf,...
        k_ns_1_est,k_ns_1_est_conf,...
        k_ns_2_est,k_ns_2_est_conf,...
        h_1_est,h_1_est_conf,...
        h_2_est,h_2_est_conf,h_est,h_est_conf] = obj.getCorrectedOffrateFun; 
    

        obj.handels.histOffDarkPlot = subplot(2,2,[1:2]','Parent', obj.handels.offRateTab);
        dwellTimesDark = obj.gdrObj.getDwellTimes(false);
        obj.handels.histOffDark  = histogram(dwellTimesDark/1000,100,'Normalization','pdf','parent',obj.handels.histOffDarkPlot);
        axis tight
        hold on
        
        x = sym('x');
        k_ns_1 = sym('k_ns_1');
        k_ns_2 = sym('k_ns_2');
        h = sym('h');
        cdffun = 1-(1-h)*exp(-k_ns_1*x)-h*exp(-k_ns_2*x);
        a = matlabFunction(diff(cdffun,x));
        ntitle(sprintf('Dark off time events \n Off-rate dark k_{ns_1} = %0.2g [s^{-1}], k_{ns_2} = %0.2g [s^{-1}] with fraction h = %0.2g\n%s',k_ns_1_est*1000,k_ns_2_est*1000,h_est,sym2str(cdffun)))
        
        x = linspace(0,max(dwellTimesDark),100);
        pdfexp = a(h_est,k_ns_1_est,k_ns_2_est,x);
        pdfexp = pdfexp*max(obj.handels.histOffDark.Values)./max(pdfexp);
        plot(x/1000,pdfexp,'--r','linewidth',2,'parent',obj.handels.histOffDarkPlot)
        xlabel('Time [s]')
        
        obj.handels.histOffPlot = subplot(2,2,[3:4]','Parent', obj.handels.offRateTab);
        dwellTimesSpecific = obj.gadObj.getDwellTimes(false);
        obj.handels.histOff = histogram(dwellTimesSpecific/1000,100,'Normalization','pdf','parent',obj.handels.histOffPlot);
        axis tight
        hold on

        x = sym('x');
        k_ns_1 = sym('k_ns_1');
        k_ns_2 = sym('k_ns_2');
        k_s = sym('k_s');
        h_1 = sym('h_1');
        h_2 = sym('h_2');

        cdffun = h_1*(1-exp(-k_ns_1*x))+h_2*(1-exp(-k_ns_2*x))+(1-h_1-h_2)*(1-exp(-k_s*x));
        ntitle(sprintf('Off time events \n Off-rate k_{s} = %0.2g [s^{-1}] with fraction 1-h_1-h_2 = %0.2g\n %s',k_s_est*1000,1-h_1_est-h_2_est,sym2str(cdffun)))
        a = matlabFunction(diff(cdffun,x));

        xt = linspace(0,max(dwellTimesSpecific),10000);
        temp = a(h_1_est,h_2_est,k_ns_1_est,k_ns_2_est,k_s_est,xt);
        temp = max(obj.handels.histOff.Values)./max(temp)*temp;
        plot(xt/1000,temp,'-b','linewidth',2,'parent',obj.handels.histOffPlot)
 
        temp = a(h_1_est./(h_1_est+h_2_est),h_2_est./(h_1_est+h_2_est),k_ns_1_est,k_ns_2_est,0,xt);
        temp = max(obj.handels.histOff.Values)./max(temp)*temp;
        plot(xt/1000,temp,'--r','linewidth',2,'parent',obj.handels.histOffPlot)
        
        temp = a(h_1_est./(h_1_est+h_2_est),h_2_est./(h_1_est+h_2_est),0,0,k_s_est,xt);
        temp = max(obj.handels.histOff.Values)./max(temp)*temp;
        plot(xt/1000,temp,'--g','linewidth',2,'parent',obj.handels.histOffPlot)
        legend('Histogram','Complete distibution','Non-specific','Specific')
        xlabel('Time [s]')
        if any(isnan([k_s_est_conf(1) k_ns_1_est_conf(1) k_ns_1_est_conf(1) h_1_est_conf(1) h_2_est_conf(1) h_est_conf(1)]))
            warningMessage=true;
        end
     
        if warningMessage
            warndlg('Not enough data to find reliable estimates!')
        end

    end
    
    function getCorrectedOnrate(obj,src,~)
        
        if ischar(src)
            switch upper(src)
                case 'PDF'
                    obj.fitcfd=0;
                case 'CDF'
                    obj.fitcfd=1;
            end
        end
        
        [k_s_est,k_s_est_conf,k_ns_est, h_est, N_ns] = obj.getCorrectedOnrateFun;

        if isfield(obj.handels,'histOnDarkPlot') && ~isempty(obj.handels.histOnDarkPlot)
            obj.deletetry(obj.handels.histOnDarkPlot)
        end
        if isfield(obj.handels,'histOnPlot') && ~isempty(obj.handels.histOnPlot)
            obj.deletetry(obj.handels.histOnPlot)
        end
        
        obj.handels.histOnDarkPlot = subplot(2,2,[1:2]','Parent', obj.handels.onRateTab);
        time_dark_t21 = obj.gdrObj.getTimes2FirstEvent(false);
        obj.handels.histOnDark = histogram(time_dark_t21/1000,100,'Normalization','pdf','parent',obj.handels.histOnDarkPlot);
        axis tight
        hold on
        ntitle(sprintf('Dark on time events \n On-rate dark %0.2g [s^{-1}]',k_ns_est*1000))
        
        x = linspace(0,max(time_dark_t21),100);
        pdfexp = k_ns_est*exp(-k_ns_est*x);
        pdfexp = pdfexp*max(obj.handels.histOnDark.Values)./max(pdfexp);
        plot(x/1000,pdfexp,'--r','linewidth',2,'parent',obj.handels.histOnDarkPlot)
        xlabel('Time [s]')
        
        obj.handels.histOnPlot = subplot(2,2,[3:4]','Parent', obj.handels.onRateTab);
        time_t21 = obj.gadObj.getTimes2FirstEvent(false);
        obj.handels.histOn = histogram(time_t21/1000,100,'Normalization','pdf','parent',obj.handels.histOnPlot);
        axis tight
        ntitle(sprintf('On time events \n On-rate events %0.2g [s^{-1}] with fraction %0.2g ',k_s_est*1000,1-h_est))
        hold on
        x = sym('x');
        k_s = sym('k_s');
        k_ns = sym('k_ns');
        h = sym('h');

        cdffun = 1-(1-h)*exp(-k_s*x)-h*exp(-(k_s+k_ns)*x);
        pdffun = matlabFunction(diff(cdffun,x));
        
        xt = linspace(0,max(time_t21),10000);
        temp = pdffun(h_est,k_ns_est,k_s_est,xt);
        temp = max(obj.handels.histOn.Values)./max(temp)*temp;
        plot(xt/1000,temp,'-b','linewidth',2,'parent',obj.handels.histOnPlot)
 
        temp = pdffun(h_est,k_ns_est,0,xt);
        temp = max(obj.handels.histOn.Values)./max(temp)*temp;
        plot(xt/1000,temp,'--r','linewidth',2,'parent',obj.handels.histOnPlot)
        
        temp = pdffun(h_est,0,k_s_est,xt);
        temp = max(obj.handels.histOn.Values)./max(temp)*temp;
        plot(xt/1000,temp,'--g','linewidth',2,'parent',obj.handels.histOnPlot)
        legend('Histogram','Complete distibution','Non-specific','Specific')
        xlabel('Time [s]')
        
        if any(isnan([k_s_est_conf(1)]))
            warndlg('Not enough data to find reliable estimates!')
        end       
    end
    
      function [k_s_est,k_s_est_conf,k_ns_est, h_est, N_ns] = getCorrectedOnrateFun(obj,initialValueMLE)
               
        k_ns_est = 1./mean(obj.gdrObj.getTimes2FirstEvent);
        rastStructDrObj = obj.gdrObj.rasterGramOuput;
        % C_NS(t) = h*N_NS*(1-exp(-k_ns*t));
        N_ns = size(rastStructDrObj.cev,1);
        h_est = sum(sum(rastStructDrObj.cev,2)>0)./N_ns;

        x = sym('x');
        k_s = sym('k_s');
        k_ns = sym('k_ns');
        h = sym('h');

        cdffun = 1-(1-h)*exp(-k_s*x)-h*exp(-(k_s+k_ns)*x);
        dwellTimes = obj.gadObj.getTimes2FirstEvent;
        if nargin < 2 || isempty(initialValueMLE)
            initial =1./(1./k_ns_est+mean(dwellTimes));
            initial(isnan(initial))=0;
            initial = min(1,max(initial,1e-6));
        else
            initial = initialValueMLE;
        end
        if obj.fitcfd
            cdffun =  matlabFunction(cdffun);
            [f,d] = ecdf(dwellTimes);
            mu = sym('mu');

            obj_func = @(mu) norm(cdffun(h_est,k_ns_est,mu(1),d) - f);
            LB = [eps*10];
            UB = [1];
            mu = fmincon(obj_func, initial,[],[],[],[],LB,UB);  
            k_s_est = mu(1);
            k_s_est_conf = 0.*mu;            
        else        
            a = matlabFunction(diff(cdffun,x));
            b = @(x,k_s) max(a(h_est,k_ns_est,k_s,x),1e-10);
            
            [k_s_est,k_s_est_conf] = mle(dwellTimes,'pdf',b,'start',initial);

            if any(isnan(k_s_est_conf))
                warndlg('Not enough data to find reliable estimate!')
            end
        end
        
        k_s_est = max(k_s_est,0);
        k_s_est_conf = max(k_s_est_conf,0);
      end  
      
      function deletetry(obj,handle)
            try
                delete(handle)
            catch
            end
        end 
   end
end