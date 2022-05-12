function outcoords = MNI_mm2vox(coords,intype)

MNItransforms = load('MNI_vox2mm_transforms.mat');
m2v = MNItransforms.m2v;
outcoords = zeros(size(coords));
switch intype
    case 'mm'
        mm = coords;
        for i=1:size(mm,1)
            outcoords(i,1:3)=mm(i,:)*m2v(1:3,1:3) + m2v(1:3,4)';
        end    
    case 'vox'
        vox = coords;
        for i=1:size(vox,1)
            outcoords(i,:) = (vox(i,1:3) - m2v(1:3,4)')/m2v(1:3,1:3);
        end        
end




