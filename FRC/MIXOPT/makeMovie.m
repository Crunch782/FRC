% Creates a movie from the images in the figures subdirectory

param;

tag = append(num2str(ss),'_',num2str(TT),'_',num2str(ReRe));

vt = append('Data/',tag,'/Movies/vort_sca.avi');
vidfile = VideoWriter(vt);
open(vidfile)

for n = 1:nf+1
   it = append('Data/',tag,'/Figures/',num2str(n),'.jpg');
   im = imread(it);
   writeVideo(vidfile, im);
end
close(vidfile)

