init()
{
FILE *fid=fopen("/tmp/init.dat","r");
fread(&r1,sizeof(float),1,fid);
fread(&r,sizeof(float),1,fid);
fread(&dtinv,sizeof(float),1,fid);
fread(&tmax,sizeof(float),1,fid);
fread(&kappa,sizeof(float),1,fid);
fread(&del,sizeof(float),1,fid);
fread(&rd,sizeof(float),1,fid);
fread(&width,sizeof(float),1,fid);
fread(&f0,sizeof(float),1,fid);
fread(&beta,sizeof(float),1,fid);
fread(&tpl,sizeof(int),1,fid);
fread(&tau[0],sizeof(float),129,fid);
fread(&Hy[0],sizeof(float),128,fid);
fread(&q[0][0][0],sizeof(float),66048,fid);
fread(&ubar[0][0],sizeof(float),258,fid);
fread(&topo[0][0],sizeof(float),33024,fid);
fclose(fid);
}
