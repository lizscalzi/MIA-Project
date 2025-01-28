function filled_array = sin_fill_gaps(logical_array,window,Fs)
% [xf] = ax_fill_gaps(x,samp_rate)
% 
% takes logical array x and sampling rate (Fs) and fills gaps between events of
% less than window seconds)  


window_samps=round(window*Fs);


% Find the indices of the true values
true_indices = find(logical_array);
    
% Initialize the filled array with the original logical array
filled_array = logical_array;
    
% Iterate through true indices to fill gaps
for i = 1:length(true_indices)-1
  % Check the gap size
  gap_size = true_indices(i+1) - true_indices(i) - 1;
        
  % If the gap is smaller than x, fill it
    if gap_size < window_samps
      filled_array(true_indices(i)+1:true_indices(i+1)-1) = true;
    end
  end
end
    
    


