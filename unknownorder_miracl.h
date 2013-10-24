

void ucheck(char *infile){
    miracl *mip;
	FILE *fp;
	InputCell *data;
	big order;
	int datac,i;

    mip=mirsys(-1024,16);  /* Use Hex internally */
	mip->IOBASE = 16;

	fp=fopen(infile,"rt");
    if (fp==NULL)
    {
        printf("file %s does not exist\n",infile);
        return;
    }

	order = mirvar(0);

	datac = 0;
	data = (InputCell*)malloc(sizeof(InputCell)*MAXSIZE);
	while (!feof(fp)) {
		data[datac].g = mirvar(0);
		data[datac].k = mirvar(0);
		data[datac].x = mirvar(0);
		data[datac].y  = 0;
		innum(data[datac].x,fp);
		fscanf(fp,"%d\r\n",&(data[datac].y) );
		innum(data[datac].g ,fp);
		innum(data[datac].k ,fp);
		datac++;
	}
	fclose(fp);

	printf("Total points in data file:%d\r\n",datac);
	//printf("How many points is being processed:");
	//scanf("%d",&datac);

	qsort(data,datac,sizeof(InputCell),cellcompare);

	big a1,a2,b1,b2,A1,A2,B1,B2,t1,t2,t3;
	int foundcollision=0;
	a1 = mirvar(0);a2 = mirvar(0);
	b1 = mirvar(0);b2 = mirvar(0);
	A1 = mirvar(0);A2 = mirvar(0);
	B1 = mirvar(0);B2 = mirvar(0);
	for(i=1;i<datac;i++) {
		//otnum(data[i-1].x,stdout);
		if (compare(data[i-1].x,data[i].x)==0) {
				
			copy(data[i-1].g,a1);
			copy(data[i-1].k,b1);

			copy(data[i].g,a2);
			copy(data[i].k,b2);
			if (data[i-1].y!=data[i].y) {
				negify(a2,a2);
				negify(b2,b2);
			}
			t1 = mirvar(0);
			t2 = mirvar(0);
			t3 = mirvar(0);
			
			char sa1[1024],sb1[1024],sa2[1024],sb2[1024],sorder[1024];;
			otstr(a1,&sa1[0]);otstr(a2,&sa2[0]);otstr(b1,&sb1[0]);otstr(b2,&sb2[0]);otstr(order,&sorder[0]);
			printf("%s+X*%s = %s+X*%s (MOD ORDER)\r\n",&sa1,&sb1,&sa2,&sb2);

			if (foundcollision==0) {
				copy(a1,A1);
				copy(b1,B1);
				copy(a2,A2);
				copy(b2,B2);
			};
			foundcollision++;
			if (foundcollision!=2) continue;

			subtract(B2,B1,t1);
			subtract(a1,a2,t2);
			multiply(t1,t2,t3);
			subtract(b2,b1,t1);
			subtract(A1,A2,t2);
			multiply(t1,t2,t2);
			subtract(t3,t2,t1);
			printf("Use external tool to factorize this number to find the order:\r\n");otnum(t1,stdout);
			printf("Input Order:");
			innum(order,stdin);
			

			if (exsign(a2)==-1) add(a2,order,a2);
			if (exsign(b2)==-1) add(b2,order,b2);
			subtract(a1,a2,t1);
			subtract(b2,b1,t2);
			
			//if (exsign(t2)==-1) neg
			//xgcd(t2,order,0,0,0);
			multi_inverse(1,&t2,order,&t3);
			multiply(t3,t1,t2);
			divide(t2,order,order);
			if (exsign(t2)==-1) add(t2,order,t2);
			cout<<"Found privatekey:";
			otnum(t2,stdout);
			//exit(0);
		}
	}
	free(data);		

}

void ugenerate(char *infile){
    FILE *fp,*fout;
    int bits;
    epoint *g,*w,*k;

	Cell walking[16];

    big a,b,p,order,gx,gy,kx,ky,privatekey,t;

	epoint *cp;
	big cg,ck;

    long seed,orderbits;
    miracl *mip;

	fp=fopen(infile,"rt");
    if (fp==NULL)
    {
        printf("file %s does not exist\n",infile);

        return;
    }
    fscanf(fp,"%d\r\n",&bits); 
    mip=mirsys(-bits/4,16);  /* Use Hex internally */

    t=mirvar(0);
    a=mirvar(0);
    b=mirvar(0);
    p=mirvar(0);
    order=mirvar(0);
    gx=mirvar(0);
    gy=mirvar(0);
    kx=mirvar(0);
    ky=mirvar(0);
    privatekey=mirvar(0);

	mip->IOBASE = 16;

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
    fscanf(fp,"%d\r\n",&seed);
	cout<<"Seed = "<<seed<<endl;
    
    fclose(fp);

    ecurve_init(a,b,p,MR_AFFINE);  /* initialise curve */

    g=epoint_init();
    k=epoint_init();
    w=epoint_init();

    if (!epoint_set(gx,gy,0,g)) /* initialise point of order q */
    {
        printf("Problem - point G(x,y) is not on the curve\n");
        exit(0);
    }
    if (!epoint_set(kx,ky,0,k)) /* initialise point of order q */
    {
        printf("Problem - point K(x,y) is not on the curve\n");
        exit(0);
    }

    ecurve_mult(order,g,w);
    if (!point_at_infinity(w))
    {
        printf("Point G(x,y) is not of the given order,take this parameter as the bitlength\n");
		orderbits = remain(order,0x80000000);
		printf("Bitlength of order:%d\r\n",orderbits);
    } else
	{
        printf("Given order is correct, use -g to generate \n");
		exit(0);
	}

	irand(seed);
	int i;
	zero(t);incr(t,1,t);
	sftbit(t,(orderbits+1)/2,t);
	epoint *p1 = epoint_init();
	for(i=0;i<16;i++) {
		walking[i].gi = mirvar(0);
		walking[i].ki = mirvar(0);
		walking[i].p  = epoint_init();

		bigrand(t,walking[i].gi);
		bigrand(t,walking[i].ki);
		ecurve_mult(walking[i].gi,g,p1);
		walking[i].p = epoint_init();
		ecurve_mult(walking[i].ki,k,walking[i].p);
		ecurve_add(p1,walking[i].p);
	}

	zero(t);incr(t,1,t);
	sftbit(t,orderbits,t);
	irand(GetTickCount());

	SetThreadPriority(GetCurrentThread(),THREAD_PRIORITY_IDLE);
	cg = mirvar(0);
	ck = mirvar(0);
	cp = epoint_init();

	fout = fopen("out.ecs","at");
	int index,t1,dp=0;
	bool notfound;
	INPUT input;
	cout<<"Speed(ops /sec):"<<setw(14)<<0;
	while (1) {
		//Initialse new trial
		bigrand(t,cg);
		bigrand(t,ck);
		ecurve_mult(cg,g,p1);
		ecurve_mult(ck,k,cp);
		ecurve_add(p1,cp);
		notfound = true;
		
		t1 = GetTickCount(); i = 0;
		while (notfound) {
			unsigned int lastlong= remain(cp->X,0xFFFF);

			if ((lastlong & 0xFFFF)==0) {
				int lsb = epoint_get(cp,gx,gy);
				//cout<<"Found a distince point:\r\n"<<endl;
				/*
				otnum(gx,stdout);
				otnum(gy,stdout);
				otnum(cg,stdout);
				otnum(ck,stdout);
				*/

				otnum(gx,fout);
				fprintf(fout,"%d\n",lsb);
				otnum(cg,fout);
				otnum(ck,fout);
				fflush(fout);
				notfound=false;
				dp++;
				if(dp%50==1) {
					memset(&input,0,sizeof(INPUT));
					input.type = INPUT_MOUSE;
					
					input.mi.dx = dp%50 %4 -2;
					input.mi.dy = dp%50 %4 -2;
					input.mi.dwFlags = MOUSEEVENTF_MOVE;
					SendInput(1,&input,sizeof(INPUT));
				}
			} else {
				index = lastlong & 0xF;

				add(cg,walking[index].gi,cg);

				add(ck,walking[index].ki,ck);

				ecurve_add(walking[index].p, cp);
				i++;
			}
		}
		cout<<"\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
		cout<<setw(7)<<i*1000/(GetTickCount()-t1+1)<<"("<<setw(5)<<dp<<")";
	}
	return;
}


