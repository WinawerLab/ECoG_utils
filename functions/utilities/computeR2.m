function R2 = computeR2(DATA, MODEL)

R2 = 1-(sum((DATA-MODEL).^2) ./ sum((DATA-mean(DATA)).^2));

end
