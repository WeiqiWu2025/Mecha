function [demdData] = survey_Mecha_funcMLDetector(possiValue,estH,ofdmDerived,ofdmDemod,lenPayload,numTags)

demdData = zeros(lenPayload,numTags);
possiValue(possiValue==0)=-1;

len = size(possiValue,1);
for idx_1 = 1:lenPayload
    min = inf;
    y = ofdmDemod(:,idx_1);
    x = ofdmDerived(:,idx_1);
    indx = 0;
    for idx_2 = 1:len
        tmp_e = norm(y-(possiValue(idx_2,:)*estH.').'.*x);
        if tmp_e<min
            min = tmp_e;
            indx = idx_2;
        end
    end
    tmp_data = possiValue(indx,:);
    tmp_data(tmp_data==-1)=0;
    demdData(idx_1,:) = tmp_data;
            
end


end

