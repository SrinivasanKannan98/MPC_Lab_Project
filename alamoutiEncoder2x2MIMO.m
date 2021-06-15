function encodedMatrix = alamoutiEncoder2x2MIMO(inp,timeSlot)

%Input 
%inp is expected either be a 1x2 vector(for encoding data symbols) or 2x2
%matrix (for encoding the channel matrix)
%TimeSlot defines the time slot in which the data symbols are being sent

%Output
%encodedMatrix is a 2x1 vector in case of data symbol encoding and a 2x2
%matrix in case of encoding of the channel matrix

    temp = zeros(size(inp));

    if(size(inp,1) == 1) %Encoding data

        %Encoding the data symbols from each antenna for different time
        %instances

        if(timeSlot == 2)
            temp(:,1) = -conj(inp(:,2));
            temp(:,2) = conj(inp(:,1));

        end

        encodedMatrix = temp(:);

    else
        %Encoding of the channel matrix in order to compute the composite
        %channel matrix
        if(timeSlot == 2)

            temp(:,1) = conj(inp(:,2));
            temp(:,2) = -conj(inp(:,1));

        end

        encodedMatrix = temp;

    end

end
    
