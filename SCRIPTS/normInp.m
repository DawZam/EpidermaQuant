function [out, avgChn, stdDevChn] = normInp(X)
    % normalize channels independently (each channel persists as a column in X).
    avgChn = mean(X);
    stdDevChn = std(X);

    % EdgeCase Condition where standard Deviation is zero of any channel
    % Modify channel's stdDev=1 as the channel is irrelevant from clustering perspective.
    zeroLoc = stdDevChn==0;
    stdDevChn(zeroLoc) = 1;
    out = (X - avgChn)./stdDevChn;
end 