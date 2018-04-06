function[peak] = FeaturesLibWindow(W, nWind, NpWind, wFreq, df, f)
    %rmsa = calc_rmsa(W, nWind)
     %----------------- LET OP: Veel overeenkomsten met elkaar------------------
    %ampvar = calc_ampvar(W, nWind, NpWind)
    %ste = calc_ste(W, nWind)      
    %-----------------------------------------------------------------------------------------
    %zcr =  calc_zcr(W, nWind, NpWind)
    peak = calc_peak(wFreq, nWind, NpWind, df);
    centroid = calc_centroid(wFreq, nWind, NpWind, f);
    %spread = calc_spread(wFreq, nWind, NpWind, centroid)
    %flatness = calc_flatness(wFreq, nWind, NpWind)
    %     kurtosis = calc_kurtosis()
    feat_vals_window = [calc_rmsa(W, nWind); calc_ampvar(W, nWind, NpWind); calc_ste(W, nWind); calc_zcr(W, nWind, NpWind); calc_peak(wFreq, nWind, NpWind, df); calc_centroid(wFreq, nWind, NpWind, f); calc_spread(wFreq, nWind, NpWind, centroid); calc_flatness(wFreq, nWind, NpWind)];
end

function E_rmsa = calc_rmsa(W, numWindows)
    W_squared = W.^2;
    E_rmsa = zeros(1,numWindows);
    for k = 1:numWindows
        E_rmsa(k) = sqrt(sum(W_squared(:,k))/length(W_squared));
    end
end

function E_ampvar = calc_ampvar(W, numWindows, N_perWindow)
    E_ampvar = zeros(1,numWindows);
    for k = 1:numWindows
        window = W(:,k);
        mu = mean(window); %mean value of the window samples
        window_var = zeros(N_perWindow,1);
        for q = 1:N_perWindow
            window_var(q) = (window(q) - mu)^2;
        end
        window_var = sum(window_var);
        E_ampvar(k) = window_var;
    end
end

function E_ste = calc_ste(W, numWindows)
    E_ste = zeros(1,numWindows);
    for i = 1:numWindows
        E_ste(i) = sum(W(:,i).^2);
    end
end

function E_zcr = calc_zcr(W, numWindows, N_perWindow)
    E_zcr = zeros(1,numWindows);
    for i = 1:numWindows
        count = 0;
        for j = 1:N_perWindow
            if j < N_perWindow
                if W(j,i)>= 0 & W(j+1,i)< 0 || W(j,i)< 0 & W(j+1,i)>= 0
                count  = count + 1;
                end
            end
        end
        E_zcr(i) = count;
    end
end

function E_peak = calc_peak(W_freq, numWindows, N_perWindow, df)
    E_peak = zeros(1,numWindows);
    for i= 1:numWindows
        current_window = abs(W_freq(1:N_perWindow/2,i)); 
        E_peak(i) = find(current_window==max(current_window)) * df;
    end
end

function E_centroid = calc_centroid(W_freq, numWindows, N_perWindow, f)
    E_centroid = zeros(1,numWindows);
    for i= 1:numWindows
        step1 = 0;
        for j = 1:N_perWindow/2
            step1 = step1 + (f(j)*W_freq(j,i));
        end
        E_centroid(i) = step1 / sum(W_freq(:,i));
    end
end

function E_spread = calc_spread(W_freq, numWindows, N_perWindow, E_centroid)
    E_spread = zeros(1,numWindows);
    for i = 1: numWindows
        step1 = 0;
        for j = 1:N_perWindow/2
            step1 = step1 + (j-E_centroid(i))^2 * W_freq(j,i)^2;
        end
        E_spread(i) = sqrt( step1 / sum(W_freq(:,i).^2) );
    end
end

function E_flatness = calc_flatness(W_freq, numWindows, N_perWindow)
    E_flatness = zeros(1,numWindows);
    for i = 1: numWindows
        step1 = 0;
        for j = 1:N_perWindow/2
            step1 = step1 + log(W_freq(j,i));
        end
        step2 = exp(step1/(N_perWindow/2));
        E_flatness(i) = step2 / (sum(W_freq(:,i)) / (N_perWindow/2));
    end
end
