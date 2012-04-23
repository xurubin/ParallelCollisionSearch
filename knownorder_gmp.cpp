#include "common_gmp.h"


void check_gmp(char *infile){
	FILE *fp;
	InputCell* data;
	CZZ order;
	int datac,i;

	fp=fopen(infile,"r");
    if (fp==NULL)
    {
        printf("file %s does not exist\n",infile);
        return;
    }

	cout<<"Input Order:";
	innum(order, stdin);

	datac = 0;
	cout<<"Allocating buffer..."<<flush;
	data = new InputCell[MAXSIZE];
	cout<<"done."<<endl<<flush;
	while (!feof(fp)) {
		if (datac >= MAXSIZE)
		{
			printf("InputCell buffer too small!\n");
			exit(1);
		}
		data[datac].g = CZZ(0);
		data[datac].k = CZZ(0);
		data[datac].x = CZZ(0);
		data[datac].y  = 0;
		innum(data[datac].x,fp);
		fscanf(fp,"%d\r\n",&(data[datac].y) );
		if ((data[datac].y != 0)&&(data[datac].y != 1))
			printf("There is an error in the datpoint file!\n");
		innum(data[datac].g ,fp);
		innum(data[datac].k ,fp);
		datac++;
	}
	fclose(fp);

	printf("Total points in data file:%d\r\n",datac);
	qsort(data,datac,sizeof(InputCell),(int (*)(const void *,const void *))cellcompare);

	CZZ a1,a2,b1,b2,t1,t2,t3;
	a1.SetFromInt(0);a2.SetFromInt(0);
	b1.SetFromInt(0);b2.SetFromInt(0);
	for(i=1;i<datac;i++) {
		//otnum(data[i-1].x,stdout);
		if (data[i-1].x == data[i].x) {
			a1 = data[i-1].g;
			b1 = data[i-1].k;
			if (data[i-1].y==data[i].y) {
				a2 = data[i].g;
				b2 = data[i].k;
			} else {
				a2 = order - data[i].g;
				b2 = order - data[i].k;
			}
			t1.SetFromInt(0);
			t2.SetFromInt(0);
			t3.SetFromInt(0);
		
			char sa1[1024],sb1[1024],sa2[1024],sb2[1024],sorder[1024];;
			a1.toString(sa1);
			a2.toString(sa2);
			b1.toString(sb1);
			b2.toString(sb2);
			order.toString(sorder);
			printf("%s+X*%s = %s+X*%s (MOD %s)\r\n",sa1,sb1,sa2,sb2,sorder);

			t1 = a1 - a2;
			t2 = b2 - b1;
			
			//if (exsign(t2)==-1) neg
			//xgcd(t2,order,0,0,0);
			CZZn::setModulus(sorder);

			t2.toString(sa1);
			(CZZn(1) / CZZn(sa1)).toString(sa1);
			t3.SetFromStr(sa1); //multi_inverse(1,&t2,order,&t3);

			t2 = (t1 * t3) % order;
			cout<<"Found privatekey:";
			otnum(t2,stdout);
		}
	}
	delete[] data;		

}

void generate_gmp(char *infile){
	FILE *fp,*fout;
	int bits;
	CECp curve;

	CZZ a,b,p,order,gx,gy,kx,ky,privatekey;

	CZZ cg[NUM_MULTIPOINTS],ck[NUM_MULTIPOINTS];

	char seed[1024];
	char buf[1024];
	fp=fopen(infile,"rt");
	if (fp==NULL)
	{
		printf("file %s does not exist\n",infile);
		return;
	}
	fscanf(fp,"%d\r\n",&bits); 

	a.SetFromInt(0);
	b.SetFromInt(0);
	p.SetFromInt(0);
	order.SetFromInt(0);
	gx.SetFromInt(0);
	gy.SetFromInt(0);
	kx.SetFromInt(0);
	ky.SetFromInt(0);
	privatekey.SetFromInt(0);

	innum(p,fp);
	cout <<"P=";otnum(p,stdout);
	innum(a,fp);
	cout <<"a=";otnum(a,stdout);
	innum(b,fp);
	cout <<"b=";otnum(b,stdout);
	innum(order,fp);
	cout <<"order=";otnum(order,stdout);
	innum(gx,fp);
	cout <<"Gx=";otnum(gx,stdout);
	innum(gy,fp);
	cout <<"Gy=";otnum(gy,stdout);
	innum(kx,fp);
	cout <<"Kx=";otnum(kx,stdout);
	innum(ky,fp);
	cout <<"Ky=";otnum(ky,stdout);
	fscanf(fp, "%s", seed); //NOTICE: \r\n is not read!
	cout<<"Seed = "<<seed<<endl;

	fclose(fp);

	CZZn::setModulus(p.toString());
	curve = CECp(a.toString(), b.toString()); /* initialise curve */

	CECPoint g(&curve, gx.toString(), gy.toString());
	if (!g.isValid()) /* initialise point of order q */
	{
		printf("Problem - point G(x,y) is not on the curve\n");
		exit(0);
	}

	CECPoint k(&curve, kx.toString(), ky.toString());
	if (!k.isValid()) /* initialise point of order q */
	{
		printf("Problem - point K(x,y) is not on the curve\n");
		exit(0);
	}

	CECPoint w(&curve);
	w = g * order;
	if (!w.isInfinity())
	{
		printf("Problem - point G(x,y) is not of the given order \n");
	}
	w = k * order;
	if (!w.isInfinity())
	{
		printf("Problem - point K(x,y) is not of the given order \n");
	}

	int i;
	Cell walking[16];
	for(i=0;i<16;i++) {
		MD5_String(seed, seed);
		walking[i].gi = CZZ(seed) % order;
		MD5_String(seed, seed);
		walking[i].ki = CZZ(seed) % order;

		walking[i].p = g * walking[i].gi + k * walking[i].ki;
		//tmppoint = g * (walking[i].gi + CZZ("123456789") * walking[i].ki);
		//if ((tmppoint.x != walking[i].p.x) || (tmppoint.y != walking[i].p.y))
		//	exit(1);
	}

	sprintf(buf, "%llu", GetTickCount());
	CZZ::RandInit(buf);

#ifdef WIN32
	SetThreadPriority(GetCurrentThread(),THREAD_PRIORITY_IDLE);
#endif
	CECMultiPoints cp(&curve);
	CECPoint* addpoints[NUM_MULTIPOINTS];

	fout = fopen("out.ecs","at");
	int index,dp=0;
	cout<<"Speed(ops /sec):"<<setw(20)<<0;

	//Initialise new trial
	for(int i=0;i<NUM_MULTIPOINTS;i++)
	{
		cg[i].Rand(order);
		ck[i].Rand(order);
		cp.set(g*cg[i] + k*ck[i], i);
	}
	int counter = 0;
	long t1 = GetTickCount();
	while (1) {

		for(int i=0;i<NUM_MULTIPOINTS;i++)
		{

			unsigned int lastlong= cp.get(i).x.Last32Bits();

			if ((lastlong & DP_MASK)==0) {
				int lsb = cp.get(i).y.Last32Bits() & 0x1;
				//cout<<"Found a distince point:\r\n"<<endl;
//if (!(g*cg[i] + k*ck[i] == cp.get(i)))
//{
//	printf("Warning: Invalid DP point!\n");
//	fflush(stdout);
//}
//cout<<buf<<endl;
//cout<<cg[i].toString()<<endl;
//cout<<ck[i].toString()<<endl;

				cp.get(i).x.toString(buf);fprintf(fout, "%s\n", buf);
				fprintf(fout,"%d\n",lsb);
				otnum(cg[i],fout);
				otnum(ck[i],fout);
				
				fflush(fout);
				dp++;

				cout<<"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
				cout<<setw(7)<<(long)counter*1000/(GetTickCount()-t1+1)<<"("<<setw(2)<<i<<")("<<setw(8)<<dp<<")"<<flush;
				
				//Reinitialise this point
				cg[i].Rand(order);
				ck[i].Rand(order);
				cp.set(g*cg[i] + k*ck[i], i);
				lastlong= cp.get(i).x.Last32Bits();

				t1 = GetTickCount(); counter = 0;
			}

			index = lastlong & 0xF;

			cg[i] += walking[index].gi;
			if (cg[i] >= order) cg[i] -= order;

			ck[i] += walking[index].ki;
			if (ck[i] >= order) ck[i] -= order;

			addpoints[i] = &walking[index].p;
		}
		cp.MultipleAdd(addpoints);
		counter += NUM_MULTIPOINTS;

	}
	return;
}



void generate_gmp_old(char *infile){
    FILE *fp,*fout;
    int bits;
	CECp curve;

    CZZ a,b,p,order,gx,gy,kx,ky,privatekey;

	CZZ cg,ck;

    char seed[1024];
	char buf[1024];
	fp=fopen(infile,"rt");
    if (fp==NULL)
    {
        printf("file %s does not exist\n",infile);
        return;
    }
    fscanf(fp,"%d\r\n",&bits); 

    a.SetFromInt(0);
    b.SetFromInt(0);
    p.SetFromInt(0);
    order.SetFromInt(0);
    gx.SetFromInt(0);
    gy.SetFromInt(0);
    kx.SetFromInt(0);
    ky.SetFromInt(0);
    privatekey.SetFromInt(0);

    innum(p,fp);
	cout <<"P=";otnum(p,stdout);
    innum(a,fp);
	cout <<"a=";otnum(a,stdout);
    innum(b,fp);
	cout <<"b=";otnum(b,stdout);
    innum(order,fp);
	cout <<"order=";otnum(order,stdout);
    innum(gx,fp);
	cout <<"Gx=";otnum(gx,stdout);
    innum(gy,fp);
	cout <<"Gy=";otnum(gy,stdout);
    innum(kx,fp);
	cout <<"Kx=";otnum(kx,stdout);
    innum(ky,fp);
	cout <<"Ky=";otnum(ky,stdout);
	fscanf(fp, "%s", seed); //NOTICE: \r\n is not read!
	cout<<"Seed = "<<seed<<endl;
    
    fclose(fp);

	CZZn::setModulus(p.toString());
	curve = CECp(a.toString(), b.toString()); /* initialise curve */

	CECPoint g(&curve, gx.toString(), gy.toString());
    if (!g.isValid()) /* initialise point of order q */
    {
        printf("Problem - point G(x,y) is not on the curve\n");
        exit(0);
    }

	CECPoint k(&curve, kx.toString(), ky.toString());
    if (!k.isValid()) /* initialise point of order q */
    {
        printf("Problem - point K(x,y) is not on the curve\n");
        exit(0);
    }

    CECPoint w(&curve);
	w = g * order;
    if (!w.isInfinity())
    {
        printf("Problem - point G(x,y) is not of the given order \n");
    }
    w = k * order;
    if (!w.isInfinity())
    {
        printf("Problem - point K(x,y) is not of the given order \n");
    }

	int i;
	Cell walking[16];
	for(i=0;i<16;i++) {
		MD5_String(seed, seed);
		walking[i].gi = CZZ(seed) % order;
		MD5_String(seed, seed);
		walking[i].ki = CZZ(seed) % order;

		walking[i].p = g * walking[i].gi + k * walking[i].ki;
		//tmppoint = g * (walking[i].gi + CZZ("123456789") * walking[i].ki);
		//if ((tmppoint.x != walking[i].p.x) || (tmppoint.y != walking[i].p.y))
		//	exit(1);
	}

	sprintf(buf, "%ld", GetTickCount());
	CZZ::RandInit(buf);

#ifdef WIN32
	SetThreadPriority(GetCurrentThread(),THREAD_PRIORITY_IDLE);
#endif
	cg.SetFromInt(0);
	ck.SetFromInt(0);
	CECPoint cp(&curve);

	fout = fopen("out.ecs","at");
	int index,dp=0;
	long t1;
	bool notfound;
	cout<<"Speed(ops /sec):"<<setw(14)<<0;
	while (1) {
		//Initialise new trial
		cg.Rand(order);
		ck.Rand(order);
		cp = g*cg + k*ck;
		notfound = true;
		
		t1 = GetTickCount(); i = 0;
		while (notfound) {

			unsigned int lastlong= cp.x.Last32Bits();

			if ((lastlong & 0xFFFF)==0) {
				int lsb = cp.y.Last32Bits() & 0x1;
				//cout<<"Found a distince point:\r\n"<<endl;
				cp.x.toString(buf);fprintf(fout, "%s\n", buf);
				fprintf(fout,"%d\n",lsb);
				otnum(cg,fout);
				otnum(ck,fout);
				fflush(fout);
				notfound=false;
				dp++;
			} else {
				index = lastlong & 0xF;

				cg += walking[index].gi;
				if (cg >= order) cg -= order;

				ck += walking[index].ki;
				if (ck >= order) ck -= order;

				cp += walking[index].p;
				i++;
			}
		}
		cout<<"\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
		cout<<setw(7)<<(long)i*1000/(GetTickCount()-t1+1)<<"("<<setw(5)<<dp<<")"<<flush;
	}
	return;
}


