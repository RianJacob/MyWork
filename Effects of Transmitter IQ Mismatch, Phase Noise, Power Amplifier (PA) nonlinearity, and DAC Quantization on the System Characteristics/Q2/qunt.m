%%% function qunt

function y  = qunt(si,code)

% si is the signal to be quantized
% code is simply the quantization levels 
s= size(si);
y = zeros(s(1),s(2));
for i=1:s(1)
    for k=1:s(2)
        if (si(i,k)<0)
            sign = -1;
        else 
            sign = 1;
        end
        
        if (si(i,k)<code(1) || si(i,k)>code(end))
            y(i,k) = code(end)*sign;
        else
            [~,index] = min(abs(code - abs(si(i,k))));
            if (si(i,k)==0)
                y(i,k) = code(index+1);
            else
                y(i,k) = sign*code(index);
            end
        end
    end
end

                
            
        
   
   
    
    
    
    
    

    
   