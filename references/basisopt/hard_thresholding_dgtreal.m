function s = hard_thresholding_dgtreal(a, k)

%account for conj-symm extension when thresholding
b=a;
b(1,:)=a(1,:)/sqrt(2);
b(size(a,1),:)=a(size(a,1),:)/sqrt(2);

% sorting them
[d, ind] = sort(abs(b), 1,'descend');
% s = zeros(length(a),1);
s = zeros(size(a));
j= 1;

while j <= size(a,2)
    if k < size(a,1)
        s(ind(1:k,j),j) = a(ind(1:k,j),j);
    else
        s = a;
    end

    j= j+1;
end

end 

