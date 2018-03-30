function [feat_vals] = FeaturesLib(y, Freq, f)
    %rmsa = calc_feat_rmsa(y)
    %----------------- LET OP: Veel overeenkomsten met elkaar------------------
    %ampvar = calc_feat_ampvar(y)
    %ste = calc_feat_ste(y)
     %-----------------------------------------------------------------------------------------
    %zcr =  calc_feat_zcr(y)
    peak = calc_feat_peak(Freq)
    centroid = calc_feat_centroid(Freq, f)
    %spread = calc_feat_spread(Freq, centroid)
    %flatness = calc_feat_flatness(Freq)
    %kurtosis = calc_feat_kurtosis(Freq, peak)
    feat_vals = [calc_feat_rmsa(y); calc_feat_ampvar(y); calc_feat_ste(y); calc_feat_zcr(y); calc_feat_peak(Freq); calc_feat_centroid(Freq, f); calc_feat_spread(Freq, centroid); calc_feat_flatness(Freq); calc_feat_kurtosis(Freq, peak)];
end

function E_rmsa = calc_feat_rmsa(y)
Y_squared = y.^2;
E_rmsa = sqrt(sum(Y_squared) / length(Y_squared));
end
    
function E_ampvar = calc_feat_ampvar(y)
    mu = mean(y);
    in_sum = zeros(size(y));
    for i = 1:size(y)
        in_sum(i) = (y(i) - mu)^2;
    end
    E_ampvar = sum(in_sum);
end

function E_ste = calc_feat_ste(y)
    E_ste = sum(y.^2);
end

function E_zcr = calc_feat_zcr(y)
    E_zcr = 0;
    n_Y= size(y);
    for i = 1:n_Y(1)-1
        if y(i) >= 0 && y(i+1) < 0 || y(i) < 0 && y(i+1) >=0
            E_zcr = E_zcr + 1;
        end
    end
end

function E_peak = calc_feat_peak(Freq)
    E_peak = max(Freq);
end

function E_centroid = calc_feat_centroid(Freq, f)
    step = 0;
    n_Freq = size(Freq);
    for i = 1:n_Freq(1)/2
        step = step + (f(i)*Freq(i));
    end
    E_centroid = step / sum(Freq);
end

function E_spread = calc_feat_spread(Freq, E_centroid)
    step = 0;
    n_Freq = size(Freq);
    for i = 1:n_Freq(1)/2
        step = step + (i - E_centroid)^2 * Freq(i)^2;
    end
    E_spread = sqrt(step / sum(Freq.^2));
end

function E_flatness = calc_feat_flatness(Freq)
    step1 = 0;
    n_Freq = size(Freq);
    for i = 1:n_Freq/2
        step1 = step1 + log(Freq(i));
    end
    step2 = exp(step1 / (n_Freq(1)/2));
    E_flatness = step2 ./ (sum(Freq) ./ (n_Freq(1)/2));
end

function E_kurtosis = calc_feat_kurtosis(Freq, E_peak)
     mu = mean(Freq);
     stdev = std(Freq);
     num_sum = 0;
     n_Freq = size(Freq);
     for i = 1:n_Freq(1)/2
        in_sum = (Freq(i) - mu)^4;
        num_sum = num_sum + in_sum;
     end
     num = 2 * num_sum;
     den = (E_peak * stdev)^4;
     E_kurtosis = num / den - 3;
end