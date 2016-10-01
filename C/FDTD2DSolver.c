#include <stdio.h>
#include <stdlib.h>

#define IE 200
#define JE 200
#define cia 79 
#define cib 119
#define cja 79 
#define cjb 119
#define IsSoftSource 1

int FDTD2DSolver()
{
	float ga[IE][JE],dz[IE][JE],ez[IE][JE],hx[IE][JE],hy[IE][JE];
	int l,n,i,j,ic,jc,nsteps,npml,nPosition;
	float ddx,dt,T,epsz,pi,epsilon,sigma,eaf;
	float xn,xxn,xnum,xd,curl_e;
	float t0,spread,pulse,sFreq;
	float gi2[IE],gi3[IE];
	float gj2[JE],gj3[JE];
	float fi1[IE],fi2[IE],fi3[IE];
	float fj1[JE],fj2[JE],fj3[JE];
	float ihx[IE][JE],ihy[IE][JE];
	FILE *fp,*fopen(),*fp1;
//------------------------------------------------------------------------------
// The Pulse source position
	ic=IE/2-5-1;
	jc=JE/2-5-1;
	ddx=.005;				// Cell size, Use m unit
	dt=ddx/3e8/2;		// Time steps, Use c=3e8 m/sec
	npml=10;
	sFreq=150.0e6;		// Central Frequency
	spread=1000.0;		// Gaussian Broaden
	epsilon=3.2*3.2;
	nPosition=3;
	nsteps=1;
	t0=40.0;
	printf("%f \n",sFreq);
//------------------------------------------------------------------------------
	T=0;
	epsz=8.8e-12;
	pi=3.14159;
//------------------------------------------------------------------------------
// initialize the arrays
	for (j=0;j<JE;j++) {
//		printf("%2d ",j);
		for (i=0;i<IE;i++) {
			dz[i][j]	=0.0;
			hx[i][j]	=0.0;
			hy[i][j]	=0.0;
			ihx[i][j]=0.0;
			ihy[i][j]=0.0;
			ga[i][j]	=1.0;
		}
//		printf (" \n");
	}
//------------------------------------------------------------------------------
	// For Microcavity (Rectangular)
	for (i=cia;i<=cib;i++) {
		for (j=cja;j<=cjb;j++) {
			ga[i][j]	=1/epsilon;
		}
	}
//------------------------------------------------------------------------------
// Calculate the PML parameters
	for (i=0;i<IE;i++) {
		gi2[i]=1.0;
		gi3[i]=1.0;
		fi1[i]=0.0;
		fi2[i]=1.0;
		fi3[i]=1.0;
	}
	for (j=0;j<JE;j++) {
		gj2[j]=1.0;	
		gj3[j]=1.0;
		fj1[j]=0.0;
		fj2[j]=1.0;
		fj3[j]=1.0;
	}
	printf("Number of PML cell--> %d \n",npml);

	for (i=0;i<=npml;i++) {
		xnum=npml-i;
		xd=npml;

		xxn=xnum/xd;
		xn=0.33*pow(xxn,3.0);
//	fi1 is different with the Program in the BOOK
		fi1[i]=xn;						fi1[IE-2-i]=xn;
		gi2[i]=1.0/(1.0+xn);			gi2[IE-1-i]=1.0/(1.0+xn);
		gi3[i]=(1.0-xn)/(1.0+xn);	gi3[IE-1-i]=(1.0-xn)/(1.0+xn);

		xxn=(xnum-.5)/xd;
		xn=0.33*pow(xxn,3.0);
//		fi1[i]=xn;						fi1[IE-2-i]=xn;
		fi2[i]=1.0/(1.0+xn);			fi2[IE-2-i]=1.0/(1.0+xn);
		fi3[i]=(1.0-xn)/(1.0+xn);	fi3[IE-2-i]=(1.0-xn)/(1.0+xn);
//		printf("%d %7.4f %7.4f  %7.4f \n",i,fi1[i],fi2[i],fi3[i]);
	}

	for (j=0;j<=npml;j++) {
		xnum=npml-j;
		xd=npml;

		xxn=xnum/xd;
		xn=0.33*pow(xxn,3.0);
//		printf("%d %7.4f %7.4f \n",j,xxn,xn);
		gj2[j]=1.0/(1.0+xn);			gj2[JE-1-j]=1.0/(1.0+xn);
		gj3[j]=(1.0-xn)/(1.0+xn);	gj3[JE-1-j]=(1.0-xn)/(1.0+xn);

		xxn=(xnum-.5)/xd;
		xn=0.33*pow(xxn,3.0);
		fj1[j]=xn;						fj1[JE-2-j]=xn;
		fj2[j]=1.0/(1.0+xn);			fj2[JE-2-j]=1.0/(1.0+xn);
		fj3[j]=(1.0-xn)/(1.0+xn);	fj3[JE-2-j]=(1.0-xn)/(1.0+xn);
	}

//	printf("gi + fi \n");
//	for (i=0;i<IE;i++) {
//		printf("%2d  %5.2f  %5.2f \n",i,gi2[i],gi3[i]);
//		printf("     %5.2f  %5.2f  %5.2f \n",fi1[i],fi2[i],fi3[i]);
//	}
//	for (j=0;j<JE;j++) {
//		printf("%2d  %5.2f  %5.2f \n",j,gj2[j],gj3[j]);
//		printf("     %5.2f  %5.2f  %5.2f \n",fj1[j],fj2[j],fj3[j]);
//	}
//	system("PAUSE");
	
	printf ("nsteps -->");
	scanf ("%d",&nsteps);
	printf ("%d \n",nsteps);
//------------------------------------------------------------------------------
	fp1=fopen("cEz2DIcJc.dat","w");
	fprintf(fp1,"IsSoftSource %d\n",IsSoftSource);
	fprintf(fp1,"Ic %d\n",ic);
	fprintf(fp1,"Jc %d\n",jc);
	fprintf(fp1,"cia %d\n",cia);
	fprintf(fp1,"cib %d\n",cib);
	fprintf(fp1,"cja %d\n",cja);
	fprintf(fp1,"cjb %d\n",cjb);
	fprintf(fp1,"npml %d\n",npml);
	fprintf(fp1,"dx %10.4e\n",ddx);
	fprintf(fp1,"sFreq %10.4e\n",sFreq);
	fprintf(fp1,"spread %8.4f\n",spread);
	fprintf(fp1,"nPos %3d\n",nPosition);
	fprintf(fp1,"Esi %8.4f\n",epsilon);
	fprintf(fp1,"nsteps %8d\n",nsteps);
//------------------------------------------------------------------------------
	while (nsteps>0) {
//------------------------------------------------------------------------------
		// Start of the main FDTD loop //
		for (n=1;n<=nsteps;n++) {
			T=T+1;
			// calculate the Dz field
			for (i=1;i<IE;i++) {
				for (j=1;j<JE;j++) {
					dz[i][j]=gi3[i]*gj3[j]*dz[i][j]
					+0.5*gi2[i]*gj2[j]*(hy[i][j]-hy[i-1][j]-hx[i][j]+hx[i][j-1]);
				}
			}

			// Sinisoidal Source
			pulse=sin(2*pi*sFreq*dt*T);
			pulse=pulse*exp(-pow(T/spread,2));
			dz[ic][jc]=dz[ic][jc]*IsSoftSource+pulse;

			// calculate the Ez field
			for (i=0;i<IE;i++) {
				for (j=0;j<JE;j++) {
					ez[i][j]=ga[i][j]*dz[i][j];
				}
			}
//			printf("%3f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f \n",T, ez[ic][jc],
//				dz[40][40],ez[40][40],hx[40][40],hy[40][40]);

			// set the ez edge to 0, as part of the PML
			for (j=0;j<JE-1;j++) {
				ez[0][j]=0.0;
				ez[IE-1][j]=0.0;
			}
			for (i=0;i<IE-1;i++) {
				ez[i][0]=0.0;
				ez[i][JE-1]=0.0;
			}

			// calculate the Hx field
			for (i=0;i<IE;i++) {
				for (j=0;j<JE-1;j++) {
					curl_e=ez[i][j]-ez[i][j+1];
					ihx[i][j]=ihx[i][j]+fi1[i]*curl_e;
					hx[i][j]=fj3[j]*hx[i][j]+fj2[j]*0.5*(curl_e+ihx[i][j]);
				}
			}

			// calculate the Hy field
			for (i=0;i<IE-1;i++) {
				for (j=0;j<JE;j++) {
					curl_e=ez[i+1][j]-ez[i][j];
					ihy[i][j]=ihy[i][j]+fj1[j]*curl_e;
					hy[i][j]=fi3[i]*hy[i][j]+fi2[i]*0.5*(curl_e+ihy[i][j]);
				}
			}
			// Save Ez to the file at each time step
			fprintf(fp1,"%8.4e \t %14.10e \t %14.10e \t %14.10e \t %14.10e \n"
				,dt*T,pulse,ez[99][99],ez[105][105],ez[105][102]);
		}
//------------------------------------------------------------------------------
// Output Results
		// Show the result on the monitor,
//		for (j=0;j<JE;j++) {
//			for (i=0;i<IE;i++) {
//				printf("%8.3e  ",ez[i][j]);
//			}
//			printf("\n");
//		}
		// Save to File
		fp=fopen("cEz2D.dat","w");
		for (j=0;j<JE;j++) {
			for (i=0;i<IE;i++) {
				fprintf(fp,"%8.3e  ",ez[i][j]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
//------------------------------------------------------------------------------
		// For Next Calculation
		printf ("nsteps -->");
		scanf ("%d",&nsteps);
		printf ("%d \n",nsteps);
//------------------------------------------------------------------------------
	}
	fclose(fp1);	// For Time Step
//------------------------------------------------------------------------------
	system("PAUSE");
	return 0;
}
