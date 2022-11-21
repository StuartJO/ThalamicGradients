function [MouseOhParc,MouseBrain] = GetMouseVolData()

if exist('./data/ancillary/MouseOhParc.mat','file') == 0 || exist('./data/ancillary/MouseBrain.mat','file') == 0
    disp('MouseOhParc.mat and MouseBrain.mat don''t appear to exist')
    disp('Generating...')
    FormatMouseAtlas; 
end
load('./data/ancillary/MouseOhParc.mat','MouseOhParc')

load('./data/ancillary/MouseBrain.mat','MouseBrain')