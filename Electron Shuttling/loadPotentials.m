function [ sparams, xxq, zzq ] = loadPotentials( sparams )
%LOADPOTENTIALS Function to either generate the potential profiles
%automatically or load them from an external file
    
%     g1 = [ones(1,60)*1.0 linspace(1.0,0.7,60) ones(1,120)*0.7];
%     g2 = [linspace(0.7,1.0,60) ones(1,120)*1.0 linspace(1.0,0.7,60)];
%     g3 = [ones(1,120)*0.7 linspace(0.7,1.0,60) ones(1,60)*1.0];
    g1 = sort([0.7:0.02:1.02, 1.01, 0.99]);
    g2 = g1;
    g3 = g2;

    % Make the waitbar to show run time
    h = waitbar(0,sprintf('On potential %0.3g, %0.3g, %0.3g',g1(1),g2(1),g3(1)),...
        'Name',sprintf('Loading potentials...'),...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    flag = 0;
    
    nn = 1;
    for ii = 1:length(g1)
        for jj = 1:length(g2)
            for kk = 1:length(g3)
                % Check for cancel button click
                if getappdata(h,'canceling')
                    flag = 1;
                    break;
                end
                
                if g1(ii) ~= 1
                    g1String = ['V_L1_' sprintf('%.3g',g1(ii))];
                else
                    g1String = 'V_L1_1.0';
                end
                if g2(jj) ~= 1
                    g2String = ['_V_C_' sprintf('%.3g',g2(jj))];
                else
                    g2String = '_V_C_1.0';
                end
                if g3(kk) ~= 1
                    g3String = ['_V_R1_' sprintf('%.3g',g3(kk))];
                else
                    g3String = '_V_R1_1.0';
                end

                currPotString = [g1String g2String g3String '.csv'];

                % Update waitbar every 3 loaded potentials
                if mod(nn,3) == 0
                    waitbar(nn/(length(g1)*length(g2)*length(g3)),...
                        h, sprintf('On potential %0.3g, %0.3g, %0.3g',g1(ii),g2(jj),g3(kk)));
                end

                
                % Load the figure
                data = dlmread([sparams.potDir currPotString]);
                [rows,cols] = size(data);
                zdata = data(2:rows,2:cols);

                % Then we are on the first potential
                if nn == 1
                    xdata = data(1,2:cols);
                    ydata = data(2:rows,1);

                    % The data may not be uniform in sampling, so we need to fix that for
                    % the fourier transforms in the main code for speed up.
                    % Find next highest power of 2 to the length of xx
                    desiredGridX = 2^(nextpow2(length(xdata)));
                    desiredGridZ = round(2.5*length(ydata));

                    % Make linearly spaced grid of points to interpolate the
                    % potential at
                    xxq = linspace(min(xdata),max(xdata),desiredGridX);
                    zzq = linspace(min(ydata),max(ydata),desiredGridZ);

                    [XX,ZZ] = meshgrid(xdata,ydata);
                    [XXq,ZZq] = meshgrid(xxq,zzq);

    %                 pots = zeros(1,desiredGridZ,desiredGridX);
                end
                sparams.potentials(nn).pot2D = interp2(XX,ZZ,zdata,XXq,ZZq);
                sparams.potentials(nn).pot2D = -sparams.potentials(nn).pot2D*sparams.ee; % Convert to J
                sparams.potentials(nn).gateValues = [g1(ii) g2(jj) g3(kk)];
                nn = nn + 1;
            end
            if flag
                break;
            end
        end
        if flag
            break;
        end
    end
    delete(h);
    xxq = xxq*1E-9; % Convert to m
    zzq = zzq*1E-9; 
end

