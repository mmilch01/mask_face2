% Thresholding algorithm based on 
% T. W. Ridler, S. Calward. Picture Thresholding Using an Iterative
% Selection Method. IEEE Trans on Systems, Man and Cybernetics, 8:1978
function[t]=select_threshold(V)
[count intensity]=hist(V(:),1024);
T(1) = round(sum(count.*intensity) ./ sum(count));
delta = 1;	% initialisation before while loop
i=1;	%counter for the generations of T (threshold)
% the index of the threshold in the intensity list (T(i) is a threshold, not an index... it can be <0, for example.

while (delta ~= 0) && (i<15)
	% after the call to the "hist" function, the intensities are sorted
	% (ascending).
	T_indexes = find(intensity >= T(i));	
	T_i = T_indexes(1);	% finds the value (in "intensity") that is closest to the threshold.
	% calculate mean below current threshold: mbt
	mbt = sum(count(1:T_i) .* intensity(1:T_i) ) ./ sum(count(1:T_i));
	% calculate mean above current threshold: mat
	mat = sum(count(T_i:end) .* intensity(T_i:end) ) ./ sum(count(T_i:end));
	% the new threshold is the mean of mat and mbt
	i= i+1;
    T(i) = round(sqrt(mbt*mat));
%	T(i) = round( 0.75*mbt+0.25*mat );
	delta = T(i) - T(i-1);
end
t = T(i);

