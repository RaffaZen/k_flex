//////////////////////////////////////////////////////////////////////
//
//   k-meros apilados y flexibles en equilibrio
//
// - 1D + 1D
// - interaccion con el sustrato V (primera capa)
// - interaccion lateral W (solo hacia abajo, capa 2 en adelante)
// - maneja lista de kuplas adsorbibles y kmeros desorbibles
//
// Marcelo Pasinetti Fecha: 12/2019 - Version: 4.XX
//////////////////////////////////////////////////////////////////////

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<limits.h>
#include<time.h>



typedef double flt;

#define ULONGMAX	4294967296.0				// unsigned long max + 1 (=2^32) para normalizar a 1 el generador
#define _for(a,b,c)	for(a=b; a<c; ++a)			// para simplificar
void erroraso(const char *msg) { printf("\nERROR: %s\n",msg); exit(EXIT_FAILURE); } // imprime un error y termina el programa

unsigned int seed[256];					        // cosas del generador pseudorandom
unsigned int r;
unsigned char irr;



#define TTMC1   5000   // numero de pasos de MC para equilibrar
#define TTMC2   5000   // numero de pasos de MC para promediar

// Rango de potencial quimico
#define UI      -15
#define UF       20
#define DU      0.2

// Rango de presion
#define pi      0.00001
#define pf      1
#define dp      0.01

// Tama�os
#define L       2000	// lado del sistema
#define HMAX    8 
#define K       2      	// tama�o del kmero

// Energias
#define c	    0.01   					// cociente q1/qi
#define W       0.0				     		// energia lateral entre particulas (solo hacia abajo)
#define V       (-1.0*log(c))/K   			// energia con el sustrato

#define namearch        "k2_c001_2C.dat"

#define busca_eventos   busca_eventos_flexible   // podes poner "busca_eventos_apilado" o "busca_eventos_flexible" segun el caso



flt U;				// potencial quimico
flt p;				// presion


int MA[L][HMAX+5];	// MATRIZ DE OCUPACION: 0:vacio, label>0:cabeza, label<0:cuerpo
int H[L];           // alturas de las columnas
int LABEL;			// label del proximo nuevo kmero que se adsorba

int N;              // numero instantaneo de sitios ocupados

flt Et;

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void inicializa_randommmm(void)
{
	int i;
// inicializo generador pseudorandom: uso semilla = segundos desde 1970, para que genera siempre una secuencia distinta
	srand((unsigned)time(0));
	irr=1;
	_for(i,0,256) seed[i]=rand();
	r=seed[0];
	_for(i,0,70000) r=seed[irr++]+=seed[r>>24];

}


flt randommmm(void)
{
//	esto genera un nuevo numero pseudoaleatorio uniformemente distribuido en el intervalo [0,1)
	r=seed[irr++]+=seed[r>>24];
	return ( (flt)r/ULONGMAX );
}



void inicializa1(void)
{
	int i,j;

	inicializa_randommmm();

	_for(i,0,L) _for(j,0,HMAX) MA[i][j]=0;  // matriz vacia

	_for(i,0,L) H[i]=0; // columnas vacias

    LABEL=1;    // contador de etiquetas de kmeros

    N=0;
    Et=0.0;


}




///////////////////////////////////////////////////////////////////////////////////////////////////
/// BUSCADOR DE EVENTOS ADS/DES (LENTO PERO SEGURO)
///////////////////////////////////////////////////////////////////////////////////////////////////

char f_ads[L],f_des[L];
int n_ads,n_des;

void busca_eventos_flexible(void)
{
    int i,j,a,b;

    _for(i,0,L) { f_ads[i]=1; f_des[i]=0; }
    n_ads=L;
    n_des=0;

    _for(i,0,L)
        if( abs( H[i] - H[(i-1+L) % L] ) > 1 )
            _for(j,1,K)
                if( f_ads[(i-j+L) % L] ) { f_ads[(i-j+L) % L]=0; --n_ads; }          // aca hay una kupla no-adsorbibles

    _for(i,0,L)
        if(f_ads[i])
            if(H[i]>0)
			{
                a=MA[i][H[i]-1];
                if(a>0)
                {
                    _for(j,1,K)
                    {
                    	b=(i+j) % L;
                    	if(MA[b][H[b]-1]!=-a) break;
                    }
                    if(j==K) { f_des[i]=1; ++n_des; }                                // aca hay un kmero desorbible
                }
			}

    // AQUI:
    // los unos de f_ads[] son kuplas adsorbibles, y n_ads es la cantidad
    // los unos de f_des[] son kmeros desorbibles, y n_des es la cantidad

}



void busca_eventos_apilado(void)
{
    int i,j,a,b;

    _for(i,0,L) f_ads[i]=f_des[i]=0;
    n_ads=0;
    n_des=0;

    _for(i,0,L) if(H[i]==0)		// busca kuplas adsorbibles en primera capa
    {
    	_for(j,1,K) if(H[(i+j) % L]>0) break;
    	if(j==K) { f_ads[i]=1; ++n_ads; }
    }
    else if(MA[i][H[i]-1]>0)						// busca kmeros y kuplas en multicapa
	{
		f_ads[i]=1,++n_ads;
		f_des[i]=1,++n_des;
	}

	// AQUI:
    // los unos de f_ads[] son kuplas adsorbibles, y n_ads es la cantidad
    // los unos de f_des[] son kmeros desorbibles, y n_des es la cantidad

}




///////////////////////////////////////////////////////////////////////////////////////////////////
/// Paso elemental de metropolis
///////////////////////////////////////////////////////////////////////////////////////////////////

void intento_adsorber(int x)
{

    int i,j,ii,jj,a,b,flag;
    flt E,q;

    E=0.0;
    _for(a,0,K)
    {
        i=(x+a) % L;
        j=H[i];

        if(j==0) E+=V; else E+=W;   // interacciona solo hacia abajo
    }

    q=exp(U-E);


    flag=0;

    if(randommmm()<q)
    {

        // ACA ADSORBO...
        _for(a,0,K)
        {
            i=(x+a) % L;
            j=H[i];
            if(a==0) MA[i][j]=LABEL; else MA[i][j]=-LABEL;
            ++H[i];

            if(H[i]>HMAX) flag=1;  //erroraso("te fuiste");

        }
		N+=K;
		Et+=E;
        ++LABEL;


        if(flag)    // vuelvo para atras
        {

            _for(a,0,K)
            {
                i=(x+a) % L;
                j=H[i];
                MA[i][j]=0;
                --H[i];
            }
            N-=K;
            Et-=E;
            --LABEL;


        }



    }




}


void intento_desorber(int x)
{
    int i,j,ii,jj,a,b;
    flt E,q;

    E=0.0;
    _for(a,0,K)
    {
        i=(x+a) % L;
        j=H[i]-1;
        if(j==0) E+=V; else E+=W;   // interacciona solo hacia abajo
    }

    q=exp(E-U);

    if(randommmm()<q)
    {

        // ACA DESORBO...
        _for(a,0,K)
        {
            i=(x+a) % L;
            j=H[i]-1;
            MA[i][j]=0;
            --H[i];
        }
		N-=K;
        Et-=E;
    }

}


void METROPOLIS(void)
{
    int i,j,a,b;


    busca_eventos();        // esto puede ser apilado o flexible


    a=(n_ads + n_des)*randommmm();	 // elije un evento al azar


	if(a<n_ads)		// intenta adsorcion
	{
		b=0;
		_for(i,0,L) // busca el a-esimo uno...
		{
			b+=f_ads[i];
			if(b>a) break;
		}
		intento_adsorber(i);
	}
	else			// intenta desorcion
	{
		a-=n_ads;
		b=0;
		_for(i,0,L) // busca el a-esimo uno...
		{
			b+=f_des[i];
			if(b>a) break;
		}
		intento_desorber(i);
	}

}







void voy_al_equilibrio(void)
{
    int i,j;

    _for(i,0,TTMC1)
    {
        _for(j,0,L) METROPOLIS();	// un MCS

    }

}


flt NACC;       // acumulador para promedio
flt ETACC;
int MUESTRAS;   // numero de muestras


void promedio_en_equilibrio(void)
{
    int i,j;

    NACC=0.0;
    ETACC=0.0;
    MUESTRAS=0;

    _for(i,0,TTMC2)
    {
        _for(j,0,L) METROPOLIS();	// un MCS

        NACC+=N;
        ETACC+=Et;
        ++MUESTRAS;


    }





}






///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/// EL MAIN
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
int main(void)
{
		int i,j,hmax;
        FILE *arch;


        inicializa1();


		printf("L=%d\nk=%d\nW=%g\nV=%g\n\n",L,K,W,V);


/***/
        arch=fopen(namearch,"wt");
        printf("U p CU \n");
        fprintf(arch,"U p CU \n");
        fclose(arch);


		for(U=UI; U<=UF; U+=DU)
//		for(p=pi; p<=pf; p+=dp)
        {

			p=exp(U);

            voy_al_equilibrio();

            promedio_en_equilibrio();

            printf("%g %g %g %g \n",U,p,(flt)NACC/MUESTRAS/(L*HMAX),(flt)ETACC/MUESTRAS/(L*HMAX));

            arch=fopen(namearch,"at");
            fprintf(arch,"%g %g %g %g \n",U,p,(flt)NACC/MUESTRAS/(L*HMAX),(flt)ETACC/MUESTRAS/(L*HMAX));
            fclose(arch);

        }

/***/


/***
  //      U=0.0;

        arch=fopen("Ncapas.dat","wt");
        _for(i,0,1)
        {
        	 _for(j,0,L) METROPOLIS();		// un MCS


//			 printf("LABEL=%d\n",LABEL);
//			 al_rest(0.1);

            if((i % 1000)==0)
            {
                    hmax=0;
                    _for(j,0,L) if(H[j]>hmax) hmax=H[j];
                    printf("%d %d\n",i,hmax);
                    fprintf(arch,"%d %d\n",i,hmax);
            }
            if(hmax>9000) break;
        }
        fclose(arch);

	//	_for(i,0,L) printf("%d",f_ads[i]); printf("\n%d\n",n_ads);
	//	_for(i,0,L) printf("%d",f_des[i]); printf("\n%d\n",n_des);

		for(;;) { }



		return 0;
***/
}

