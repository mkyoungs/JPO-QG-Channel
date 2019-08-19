function Cr(typ,nm,val)
%   Cr('filename.h','init.c');
%     opens files (to be included in main C program)
%     first is global defines, second is a reading routine init()
%     data is written to /tmp/init.dat
%     to be called at beginning of mail program
%   Cr('#define','varname',value);
%     constructs a cpp directive -- no semicolon at end
%   Cr('int','varname',value);
%   Cr('float','varname',value);
%   Cr('double','varname',value);
%     defines and initializes a scalar/vector/matrix
%     with the specified values
%   Cr()
%     closes files
%
global fidC fidCr fidCd;
if nargin==0
  fclose(fidC);
  fprintf(fidCr,'fclose(fid);\n}\n');
  fclose(fidCr);
  fclose(fidCd);
  clear fidC;
  clear fidCr;
  clear fidCd;
  return;
elseif nargin==2
  fidC=fopen(typ,'w');
  fidCr=fopen(nm,'w');
  fprintf(fidCr,'init()\n{\nFILE *fid=fopen("/tmp/init.dat","r");\n');
  fidCd=fopen('/tmp/init.dat','w');
  return;
end

if typ(1)=='#'
  fprintf(fidC,'%s %s %g\n',typ,nm,val);
  return;
end

sz=size(val);
nd=sum(sz>1);

if nd==0
  fprintf(fidC,'%s %s;\n',typ,nm);
  fprintf(fidCr,'fread(&%s,sizeof(%s),1,fid);\n',nm,typ);
  fwrite(fidCd,val,typ);
  return;
end

if nd==1
  l=length(val);
  fprintf(fidC,'%s %s[%d];\n',typ,nm,l);
  fprintf(fidCr,'fread(&%s[0],sizeof(%s),%d,fid);\n',nm,typ,l);
  fwrite(fidCd,val,typ);
  return;
end

% C array val[y][x]
%val=reshape(val',prod(sz),1);
% C array [x][y]
val=reshape(val,prod(sz),1);
switch (nd)
  case 2
    fprintf(fidC,'%s %s[%d][%d];\n',typ,nm,sz(2),sz(1));
    fprintf(fidCr,'fread(&%s[0][0],sizeof(%s),%d,fid);\n',nm,typ,length(val));
  case 3
    fprintf(fidC,'%s %s[%d][%d][%d];\n',typ,nm,sz(3),sz(2),sz(1));
    fprintf(fidCr,'fread(&%s[0][0][0],sizeof(%s),%d,fid);\n',nm,typ,length(val));
  case 4
    fprintf(fidC,'%s %s[%d][%d][%d][%d];\n',typ,nm,sz(4),sz(3),sz(2),sz(1));
    fprintf(fidCr,'fread(&%s[0][0][0][0],sizeof(%s),%d,fid);\n',nm,typ,length(val));
end
fwrite(fidCd,val,typ);
return;

end

