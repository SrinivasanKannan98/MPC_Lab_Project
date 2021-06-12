function encodedMatrix = alamoutiEncoder4x4MIMO(inp,timeSlot)

%Input 
%inp is expected either be a 1x4 vector(for encoding data symbols) or 4x4
%matrix (for encoding the channel matrix)
%TimeSlot defines the time slot in which the data symbols are being sent

%Output
%encodedMatrix is a 4x1 vector in case of data symbol encoding and a 4x4
%matrix in case of encoding of the channel matrix

    temp = zeros(size(inp));

    if(size(inp,1) == 1) %Encoding data

        %Encoding the data symbols from each antenna for different time
        %instances

        if(timeSlot == 2)
            temp(:,1) = conj(inp(:,2));
            temp(:,2) = -conj(inp(:,1));
            temp(:,3) = conj(inp(:,4));
            temp(:,4) = -conj(inp(:,3));
        elseif(timeSlot == 3)
            temp(:,1) = conj(inp(:,3));
            temp(:,2) = conj(inp(:,4));
            temp(:,3) = -conj(inp(:,1));
            temp(:,4) = -conj(inp(:,2));
        elseif(timeSlot == 4)
            temp(:,1) = inp(:,4);
            temp(:,2) = -inp(:,3);
            temp(:,3) = -inp(:,2);
            temp(:,4) = inp(:,1);
        end

        encodedMatrix = temp(:);

    else
        %Encoding of the channel matrix in order to compute the composite
        %channel matrix

        if(timeSlot == 2)
            temp(:,2) = conj(inp(:,1));
            temp(:,1) = -conj(inp(:,2));
            temp(:,4) = conj(inp(:,3));
            temp(:,3) = -conj(inp(:,4));
        elseif(timeSlot == 3)
            temp(:,3) = conj(inp(:,1));
            temp(:,4) = conj(inp(:,2));
            temp(:,1) = -conj(inp(:,3));
            temp(:,2) = -conj(inp(:,4));
        elseif(timeSlot == 4)
            temp(:,4) = inp(:,1);
            temp(:,3) = -inp(:,2);
            temp(:,2) = -inp(:,3);
            temp(:,1) = inp(:,4);
        end

        encodedMatrix = temp;

    end

end
    
    