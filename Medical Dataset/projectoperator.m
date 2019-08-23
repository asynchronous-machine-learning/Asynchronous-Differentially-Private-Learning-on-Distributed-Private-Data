function po = projectoperator(theta,bound)

i_out_v=[];
i_out_v=find(abs(theta)>bound);

po=theta;
for i = 1:length(i_out_v)
    po(i_out_v(i))=bound*sign(theta(i_out_v(i)));
end

