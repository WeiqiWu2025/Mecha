function [demodData] = survey_Mecha_funcDemd(estH,ofdmDerived,ofdmDemod)
global numTags;

lenPayload = size(ofdmDemod,2);
% demodData = zeros(lenPayload,numTags);

possiValue = survey_Mecha_funcPossibleValues(numTags);
demodData = survey_Mecha_funcMLDetector(possiValue,estH,ofdmDerived,ofdmDemod,lenPayload,numTags);

end

