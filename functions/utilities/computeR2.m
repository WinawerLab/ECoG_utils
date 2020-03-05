function R2 = computeR2(DATA, MODEL)

% Compute measure of R2 between DATA and MODEL based on the SSE (sum of
% squares) of the residuals between model fit and data

R2 = 1 - sum((DATA-MODEL).^2) ./ sum((DATA-mean(DATA)).^2);

end
