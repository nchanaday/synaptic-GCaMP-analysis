for nfile=1:6
    if nfile==1
        in_filename='file1.xlsx';
        out='file1_Result';
    elseif nfile==2
        in_filename='file2.xlsx';
        out='file2_Result';
    elseif nfile==3
        in_filename='file3.xlsx';
        out='file3_Result';
    elseif nfile==4
        in_filename='file4.xlsx';
        out='file4_Result';
    elseif nfile==5
        in_filename='file5.xlsx';
        out='file5_Result';
    elseif nfile==6
        in_filename='file6.xlsx';
        out='file6_Result';    
    end

    [status,sheets]=xlsfinfo(in_filename);
    NSheets=numel(sheets);

    for sheet=1:NSheets
        out_filename=[out,num2str(sheet)];

        Data=xlsread(in_filename,sheet);
        t=Data(:,2);
        ROI=Data(:,5:end);
        t_APs=Data(1,1); % single AP (low frequency) stimualtion
        t_20Hz=Data(2,1); % high frequency stimulation (HFS)
        t_KCl=Data(3,1); % perfusion of 90mM KCl

        i_APs=find(t>t_APs,1);
        i_20Hz_start=find(t>(t_20Hz),1);
        i_20Hz_end=find(t>(t_20Hz+5),1);  % the +5 is to find the end of the peak, since the stimulation lasts 5 seconds
        i_KCl=find(t>t_KCl,1);

        acq_rate=t(end)/size(t,1); % acquisition rate (inverse of exposure/interval time)
        APs_interval=round(5/acq_rate); % interval between single AP stimulations (inverse of stimulation frequency)
        APs_number=20; % number of APs delivered - 10 for Syb2 and 20 for PSD95, same freq.

        % Amplitude of the reaponse to 90mM KCl:
        Ampli_90K=max(ROI(i_KCl:end,:),[],1)-mean(ROI(i_KCl-50:i_KCl,:),1);

        % this part creates just matrices filled with NaNs (to make code faster)
        Amplitude_singleAPs=NaN(APs_number,size(ROI,2));
        Tau=NaN(APs_number,size(ROI,2));
        Probability=NaN(1,size(ROI,2));
        sp=1;
        single_peaks=NaN(200,100);
        Tau_perROI=NaN(1,size(ROI,2));
        Ampli_perROI=NaN(1,size(ROI,2));
        Ampli_HFS=NaN(1,size(ROI,2));
        HFS_during_slope=NaN(1,size(ROI,2));
        HFS_after_Tau=NaN(1,size(ROI,2));
        HFS_after_base=NaN(1,size(ROI,2));
        HFS_peaks=NaN(size(i_20Hz_start-50:i_20Hz_end+600,2),size(ROI,2));

        %% ---> DENOISING algorithm

        % the values of these parameters were set after many tests (good denoising with little peak deformation)
        N = [4 8 16 32];                                                           % Sampling window length (size of the forward and backward mean for the predictors)
        K = length(N);                                                             % Number of forward and backward predictors
        M = 25;                                                                    % Analysis window to compare the predictors
        P = 2;                                                                     % Weighting factor

        % --- Forward-Backward non-linear Algorithm

        for y=1:size(ROI,2)
            testROI(:,1)=ROI(:,y);
            ltime=length(t);
            I_avg_f = zeros(ltime,K);
            I_avg_b = zeros(ltime,K);

            for g = 1:ltime
                for k = 1:K
                    % Average forward predictor
                    window = N(k);
                    if g == 1
                        I_avg_f(g,k) = testROI(1,1);
                    elseif g - window - 1 < 0
                        I_avg_f(g,k) = sum(testROI(1:g-1,1))/g;
                    else
                        epoint = g - window;
                        spoint = g - 1;
                        I_avg_f(g,k) = sum(testROI(epoint:spoint,1))/window;
                    end
                    % Average backward predictor
                    if g == ltime
                        I_avg_b(g,k) = testROI(g,1);
                    elseif g + window > ltime
                        sw = ltime - g;
                        I_avg_b(g,k) = sum(testROI(g+1:ltime,1))/sw;
                    else
                        epoint = g + window;
                        spoint = g + 1;
                        I_avg_b(g,k) = sum(testROI(spoint:epoint,1))/window;
                    end
                end
            end

            % Non-normalized forward and backward weights:
            f = zeros(ltime,K);
            b = zeros(ltime,K);
            for i = 1:ltime
                for k = 1:K
                    Mstore_f = zeros(M,1);
                    Mstore_b = zeros(M,1);
                    Pi=1/(2*K); % Natali added this (see paper by Chung and Kennedy)
                    for j = 0:M-1
                        t_f = i - j;
                        t_b = i + j;
                        if t_f < 1
                            Mstore_f(j+1,1) = (testROI(i,1) - I_avg_f(i,k))^2;
                        else
                            Mstore_f(j+1,1) = (testROI(t_f,1) - I_avg_f(t_f,k))^2;
                        end
                        % eqn. (4) and (5) in paper by Chung ang Kennedy:
                        if t_b > ltime
                            Mstore_b(j+1,1) = (testROI(i,1) - I_avg_b(i,k))^2;
                        else
                            Mstore_b(j+1,1) = (testROI(t_b,1) - I_avg_b(t_b,k))^2;
                        end
                    end
                    f(i,k) = Pi*(sum(Mstore_f)^(-P));
                    b(i,k) = Pi*(sum(Mstore_b)^(-P));
                end
            end

            % Vector of normalization factors for the weights:
            C = zeros(ltime,1);
            for i = 1:ltime
                Kstore = zeros(K,1);
                for k = 1:K
                    Kstore(k,1) = f(i,k) + b(i,k);
                end
                C(i,1) = 1/sum(Kstore);
            end

            % Putting parameters together and solving for intensities:
            ROIclean = zeros(ltime,1);
            for i = 1:ltime
                TempSum = zeros(K,1);
                for k = 1:K
                    TempSum(k,1) = f(i,k)*C(i,1)*I_avg_f(i,k) + b(i,k)*C(i,1)*I_avg_b(i,k);
                    % summatory over K of eqn. (2) in paper by Chung and Kennedy
                end
                ROIclean(i,1) = sum(TempSum);
            end

            ROI(2:ltime-1,y) = ROIclean(2:ltime-1,1);

        end

        ROI(ROI==Inf)=500;
        ROI(isnan(ROI))=0;

        ROI=(ROI-mean(ROI(1:100,:),1))./Ampli_90K;


        %% ---> Finding and fitting single AP PEAKs + Analysis of HFS

        for j=1:size(ROI,2) % loop through the ROIs (columns)
            if Ampli_90K(j)<100
                continue
            end
            N=0; sp_first=sp;

            for i=1:APs_number % loop through the rows (just the time points where single APs were delivered)
                point=i_APs+APs_interval*(i-1);
                [Max,Ipoint]=max(ROI((point-1:point+9),j));
                point=point+Ipoint-1;
                peak=mean(ROI((point:point+4),j),1);
                base=mean(ROI((point-30:point-5),j),1);
                noise=2*std(ROI((point-30:point-5),j),1);
                base2=mean(ROI((point-60:point-31),j),1);

                % conditional loop to detect peaks and calculate their amplitude, all peaks are saved in the "single_peaks" matrix
                % used a threshold of 0.04 for Syb2 and 0.06 for PSD95
                if (peak>base+noise) && (base2<=base+0.5*noise) && (peak-base>0.06) 
                    N=N+1;
                    Amplitude_singleAPs(i,j)=peak-base;
                    single_peaks(:,sp)=ROI(point-50:point+149,j)-base;
                    sp=sp+1;
                end

            end

            % probability of detection of single AP peaks (aka failure analysis)
            Probability(1,j)=N/APs_number;
            % mean Amplitude of all single AP peaks from each ROI
            Ampli_perROI=mean(Amplitude_singleAPs,1,'omitnan');

            % to fit the decay of the average of all single AP peaks of each ROI
            if sp-sp_first>3
                tfit=t(1:151,1);
                ROIfit=mean(single_peaks(50:200,sp_first:sp-1),2);
                options=optimoptions(@lsqcurvefit,'MaxFunEvals',1000,'MaxIter',1000,'Algorithm','levenberg-marquardt'); lb=[]; ub=[];
                expdecay=@(x,time)ROIfit(1)*exp(-x*time); x0=1;
                APs_expfit=lsqcurvefit(expdecay,x0,tfit,ROIfit,lb,ub,options);
                Tau_perROI(j)=1/APs_expfit;
            end

            % analysis if high frequency peak (HFS): amplitude, linear fit of peak during stimualtion and exponential fot of decay after the stimulation
            Ampli_HFS(j)=(max(ROI((i_20Hz_start:i_20Hz_end+10),j))-mean(ROI((i_20Hz_start-50:i_20Hz_start-5),j)));

            tfit_A=t(i_20Hz_start+50:i_20Hz_end-10,1);
            ROIfit_A=ROI(i_20Hz_start+50:i_20Hz_end-10,j)-mean(ROI((i_20Hz_start-50:i_20Hz_start-5),j));
            options=optimoptions(@lsqcurvefit,'MaxFunEvals',1000,'MaxIter',1000,'Algorithm','levenberg-marquardt'); lb=[]; ub=[];
            linfit=@(x,time)x(1)*time+x(2); x0=[-1,500];
            HFS_linfit=lsqcurvefit(linfit,x0,tfit_A,ROIfit_A,lb,ub,options);
            HFS_during_slope(j)=HFS_linfit(1);

            tfit_B=t(i_20Hz_end:i_20Hz_end+599,1);
            ROIfit_B=ROI(i_20Hz_end:i_20Hz_end+599,j)-mean(ROI((i_20Hz_start-50:i_20Hz_start-5),j));
            options=optimoptions(@lsqcurvefit,'MaxFunEvals',1000,'MaxIter',1000,'Algorithm','levenberg-marquardt'); lb=[]; ub=[];
            expdecay=@(x,time)ROIfit_B(1)*exp(-x*time); x0=0.0001;
            HFS_expfit=lsqcurvefit(expdecay,x0,tfit_B,ROIfit_B,lb,ub,options);
            HFS_after_Tau(j)=1/HFS_expfit;
            HFS_after_base(j)=mean(ROIfit_B(500:600));

            % this just saves all the HFS peaks
            HFS_peaks(:,j)=ROI(i_20Hz_start-50:i_20Hz_end+600,j)-mean(ROI((i_20Hz_start-50:i_20Hz_start-5),j));

        end

        % to fit the decay of the average of all single AP peaks (whole sheet - all ROIs)
        if ~isnan(mean(mean(single_peaks,2,'omitnan')))==1
        tfit=t(1:151,1);
        ROIfit=mean(single_peaks(50:200,1:sp-1),2,'omitnan');
        options=optimoptions(@lsqcurvefit,'MaxFunEvals',1000,'MaxIter',1000,'Algorithm','levenberg-marquardt'); lb=[]; ub=[];
        expdecay=@(x,time)ROIfit(1)*exp(-x*time); x0=1;
        APs_expfit=lsqcurvefit(expdecay,x0,tfit,ROIfit,lb,ub,options);
        Tau_singleAP=1/APs_expfit;
        else
          Tau_singleAP=NaN;  
        end

        if ~isnan(mean(mean(HFS_peaks,2,'omitnan')))==1
        % to fit the decay of the average of all HFS peaks (whole sheet - all ROIs)
        tfit_C=t(1:601,1);
        ROIfit_C=mean(HFS_peaks(end-600:end,:),2,'omitnan');
        options=optimoptions(@lsqcurvefit,'MaxFunEvals',1000,'MaxIter',1000,'Algorithm','levenberg-marquardt'); lb=[]; ub=[];
        expdecay=@(x,time)ROIfit_C(1)*exp(-x*time); x0=0.0001;
        HFS_expfit2=lsqcurvefit(expdecay,x0,tfit_C,ROIfit_C,lb,ub,options);
        Tau_meanHFS=1/HFS_expfit2;
        else
          Tau_meanHFS=NaN;
        end

        % to save the results (as MAT file)
        save(out_filename,"Amplitude_singleAPs","Probability","single_peaks","Tau_perROI","Ampli_perROI", ...
            "Ampli_HFS","HFS_during_slope","HFS_after_Tau","Tau_singleAP","Ampli_90K", ...
            "HFS_peaks","Tau_meanHFS","HFS_after_base");

        clearvars -except in_filename out status sheets NSheets sheet nfile
    end

end
clear