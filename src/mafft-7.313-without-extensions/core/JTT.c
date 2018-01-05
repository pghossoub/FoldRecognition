#if 0
#include "mltaln.h"
#endif
#define DEFAULTGOP_J -1530
#define DEFAULTGEP_J   -00 
#define DEFAULTOFS_J  -123  /* +10 -- -50  teido ka ? */
#define DEFAULTPAMN  200

void JTTmtx( double **rsr, double *freq, unsigned char locamino[26], char locgrp[26], int isTM )
{
	int i, j;
	double r[20][20];
//	char locamino0[] = "ARNDCQEGHILKMFPSTWYVBZX.-U";
	char locamino0[] = "ARNDCQEGHILKMFPSTWYVBZX.-J";
	char locgrp0[] = 
	{
		0, 3, 2, 2, 5, 2, 2, 0, 3, 1, 1, 3, 1, 4, 0, 0, 0, 4, 4, 1, 2, 2,
		6, 6, 6, 1, 
	};

	double freq0[20] = 
	{
		0.077,
		0.051,
		0.043,
		0.052,
		0.020,
		0.041,
		0.062,
		0.074,
		0.023,
		0.052,
		0.091,
		0.059,
		0.024,
		0.040,
		0.051,
		0.069,
		0.059,
		0.014,
		0.032,
		0.066,
	};
	double freq0_TM[20] = 
	{
		 0.1051,
		 0.0157,
		 0.0185,
		 0.0089,
		 0.0219,
		 0.0141,
		 0.0097,
		 0.0758,
		 0.0168,
		 0.1188,
		 0.1635,
		 0.0112,
		 0.0333,
		 0.0777,
		 0.0260,
		 0.0568,
		 0.0523,
		 0.0223,
		 0.0324,
		 0.1195,
	};

    /* Lower triangular is JTT's Accepted point mutations */
    r[ 1][ 0]=  247;

    r[ 2][ 0]=  216; r[ 2][ 1]=  116;

    r[ 3][ 0]=  386; r[ 3][ 1]=   48; r[ 3][ 2]= 1433;

    r[ 4][ 0]=  106; r[ 4][ 1]=  125; r[ 4][ 2]=   32; r[ 4][ 3]=   13;

    r[ 5][ 0]=  208; r[ 5][ 1]=  750; r[ 5][ 2]=  159; r[ 5][ 3]=  130;
    r[ 5][ 4]=    9;

    r[ 6][ 0]=  600; r[ 6][ 1]=  119; r[ 6][ 2]=  180; r[ 6][ 3]= 2914;
    r[ 6][ 4]=    8; r[ 6][ 5]= 1027;

    r[ 7][ 0]= 1183; r[ 7][ 1]=  614; r[ 7][ 2]=  291; r[ 7][ 3]=  577;
    r[ 7][ 4]=   98; r[ 7][ 5]=   84; r[ 7][ 6]=  610;

    r[ 8][ 0]=   46; r[ 8][ 1]=  446; r[ 8][ 2]=  466; r[ 8][ 3]=  144;
    r[ 8][ 4]=   40; r[ 8][ 5]=  635; r[ 8][ 6]=   41; r[ 8][ 7]=   41;

    r[ 9][ 0]=  173; r[ 9][ 1]=   76; r[ 9][ 2]=  130; r[ 9][ 3]=   37;
    r[ 9][ 4]=   19; r[ 9][ 5]=   20; r[ 9][ 6]=   43; r[ 9][ 7]=   25;
    r[ 9][ 8]=   26;

    r[10][ 0]=  257; r[10][ 1]=  205; r[10][ 2]=   63; r[10][ 3]=   34;
    r[10][ 4]=   36; r[10][ 5]=  314; r[10][ 6]=   65; r[10][ 7]=   56;
    r[10][ 8]=  134; r[10][ 9]= 1324;

    r[11][ 0]=  200; r[11][ 1]= 2348; r[11][ 2]=  758; r[11][ 3]=  102;
    r[11][ 4]=    7; r[11][ 5]=  858; r[11][ 6]=  754; r[11][ 7]=  142;
    r[11][ 8]=   85; r[11][ 9]=   75; r[11][10]=   94;

    r[12][ 0]=  100; r[12][ 1]=   61; r[12][ 2]=   39; r[12][ 3]=   27;
    r[12][ 4]=   23; r[12][ 5]=   52; r[12][ 6]=   30; r[12][ 7]=   27;
    r[12][ 8]=   21; r[12][ 9]=  704; r[12][10]=  974; r[12][11]=  103;

    r[13][ 0]=   51; r[13][ 1]=   16; r[13][ 2]=   15; r[13][ 3]=    8;
    r[13][ 4]=   66; r[13][ 5]=    9; r[13][ 6]=   13; r[13][ 7]=   18;
    r[13][ 8]=   50; r[13][ 9]=  196; r[13][10]= 1093; r[13][11]=    7;
    r[13][12]=   49;

    r[14][ 0]=  901; r[14][ 1]=  217; r[14][ 2]=   31; r[14][ 3]=   39;
    r[14][ 4]=   15; r[14][ 5]=  395; r[14][ 6]=   71; r[14][ 7]=   93;
    r[14][ 8]=  157; r[14][ 9]=   31; r[14][10]=  578; r[14][11]=   77;
    r[14][12]=   23; r[14][13]=   36;

    r[15][ 0]= 2413; r[15][ 1]=  413; r[15][ 2]= 1738; r[15][ 3]=  244;
    r[15][ 4]=  353; r[15][ 5]=  182; r[15][ 6]=  156; r[15][ 7]= 1131;
    r[15][ 8]=  138; r[15][ 9]=  172; r[15][10]=  436; r[15][11]=  228;
    r[15][12]=   54; r[15][13]=  309; r[15][14]= 1138;

    r[16][ 0]= 2440; r[16][ 1]=  230; r[16][ 2]=  693; r[16][ 3]=  151;
    r[16][ 4]=   66; r[16][ 5]=  149; r[16][ 6]=  142; r[16][ 7]=  164;
    r[16][ 8]=   76; r[16][ 9]=  930; r[16][10]=  172; r[16][11]=  398;
    r[16][12]=  343; r[16][13]=   39; r[16][14]=  412; r[16][15]= 2258;

    r[17][ 0]=   11; r[17][ 1]=  109; r[17][ 2]=    2; r[17][ 3]=    5;
    r[17][ 4]=   38; r[17][ 5]=   12; r[17][ 6]=   12; r[17][ 7]=   69;
    r[17][ 8]=    5; r[17][ 9]=   12; r[17][10]=   82; r[17][11]=    9;
    r[17][12]=    8; r[17][13]=   37; r[17][14]=    6; r[17][15]=   36;
    r[17][16]=    8;

    r[18][ 0]=   41; r[18][ 1]=   46; r[18][ 2]=  114; r[18][ 3]=   89;
    r[18][ 4]=  164; r[18][ 5]=   40; r[18][ 6]=   15; r[18][ 7]=   15;
    r[18][ 8]=  514; r[18][ 9]=   61; r[18][10]=   84; r[18][11]=   20;
    r[18][12]=   17; r[18][13]=  850; r[18][14]=   22; r[18][15]=  164;
    r[18][16]=   45; r[18][17]=   41;

    r[19][ 0]= 1766; r[19][ 1]=   69; r[19][ 2]=   55; r[19][ 3]=  127;
    r[19][ 4]=   99; r[19][ 5]=   58; r[19][ 6]=  226; r[19][ 7]=  276;
    r[19][ 8]=   22; r[19][ 9]= 3938; r[19][10]= 1261; r[19][11]=   58;
    r[19][12]=  559; r[19][13]=  189; r[19][14]=   84; r[19][15]=  219;
    r[19][16]=  526; r[19][17]=   27; r[19][18]=   42;


    /* Upper triangular is JTT's Accepted point mutations for transmembrane */
 r[ 0][ 1]=   21; r[ 0][ 2]=    2; r[ 0][ 3]=    7; r[ 0][ 4]=   13;
 r[ 0][ 5]=    4; r[ 0][ 6]=    6; r[ 0][ 7]=  160; r[ 0][ 8]=    6;
 r[ 0][ 9]=   44; r[ 0][10]=   43; r[ 0][11]=    5; r[ 0][12]=   10;
 r[ 0][13]=   21; r[ 0][14]=   34; r[ 0][15]=  198; r[ 0][16]=  202;
 r[ 0][17]=    0; r[ 0][18]=    1; r[ 0][19]=  292; 
 
 r[ 1][ 2]=    0; r[ 1][ 3]=    1; r[ 1][ 4]=    2; r[ 1][ 5]=   21;
 r[ 1][ 6]=    3; r[ 1][ 7]=   22; r[ 1][ 8]=   21; r[ 1][ 9]=    4;
 r[ 1][10]=    8; r[ 1][11]=   53; r[ 1][12]=   19; r[ 1][13]=    0;
 r[ 1][14]=    1; r[ 1][15]=    5; r[ 1][16]=    5; r[ 1][17]=   28;
 r[ 1][18]=    0; r[ 1][19]=    0; 
 
 r[ 2][ 3]=   14; r[ 2][ 4]=    1; r[ 2][ 5]=    7; r[ 2][ 6]=    0;
 r[ 2][ 7]=    0; r[ 2][ 8]=    8; r[ 2][ 9]=    4; r[ 2][10]=    5;
 r[ 2][11]=   11; r[ 2][12]=    3; r[ 2][13]=    1; r[ 2][14]=    2;
 r[ 2][15]=   32; r[ 2][16]=   19; r[ 2][17]=    1; r[ 2][18]=    1;
 r[ 2][19]=    2; 
 
 r[ 3][ 4]=    0; r[ 3][ 5]=    0; r[ 3][ 6]=   12; r[ 3][ 7]=   15;
 r[ 3][ 8]=    4; r[ 3][ 9]=    1; r[ 3][10]=    0; r[ 3][11]=    2;
 r[ 3][12]=    1; r[ 3][13]=    0; r[ 3][14]=    1; r[ 3][15]=    0;
 r[ 3][16]=    6; r[ 3][17]=    0; r[ 3][18]=    1; r[ 3][19]=    4;
 
 r[ 4][ 5]=    0; r[ 4][ 6]=    0; r[ 4][ 7]=   13; r[ 4][ 8]=    2;
 r[ 4][ 9]=    4; r[ 4][10]=   11; r[ 4][11]=    0; r[ 4][12]=    1;
 r[ 4][13]=   34; r[ 4][14]=    0; r[ 4][15]=   48; r[ 4][16]=   13;
 r[ 4][17]=    8; r[ 4][18]=   23; r[ 4][19]=   47; 
 
 r[ 5][ 6]=   16; r[ 5][ 7]=    1; r[ 5][ 8]=   26; r[ 5][ 9]=    1;
 r[ 5][10]=   16; r[ 5][11]=    6; r[ 5][12]=    3; r[ 5][13]=    0;
 r[ 5][14]=    5; r[ 5][15]=    7; r[ 5][16]=    2; r[ 5][17]=    0;
 r[ 5][18]=    0; r[ 5][19]=    0; 
 
 r[ 6][ 7]=   21; r[ 6][ 8]=    0; r[ 6][ 9]=    0; r[ 6][10]=    0;
 r[ 6][11]=    0; r[ 6][12]=    0; r[ 6][13]=    0; r[ 6][14]=    0;
 r[ 6][15]=    4; r[ 6][16]=    2; r[ 6][17]=    0; r[ 6][18]=    0;
 r[ 6][19]=    7; 
 
 r[ 7][ 8]=    1; r[ 7][ 9]=   10; r[ 7][10]=    0; r[ 7][11]=    0;
 r[ 7][12]=    3; r[ 7][13]=    4; r[ 7][14]=    7; r[ 7][15]=   64;
 r[ 7][16]=   12; r[ 7][17]=    5; r[ 7][18]=    0; r[ 7][19]=   53;
 
 r[ 8][ 9]=    3; r[ 8][10]=    2; r[ 8][11]=    0; r[ 8][12]=    1;
 r[ 8][13]=    0; r[ 8][14]=    0; r[ 8][15]=    0; r[ 8][16]=    4;
 r[ 8][17]=    0; r[ 8][18]=   29; r[ 8][19]=    2;

 r[ 9][10]=  273; r[ 9][11]=    0; r[ 9][12]=  161; r[ 9][13]=   66;
 r[ 9][14]=    4; r[ 9][15]=   22; r[ 9][16]=  150; r[ 9][17]=    1;
 r[ 9][18]=    4; r[ 9][19]=  883;

 r[10][11]=    1; r[10][12]=  153; r[10][13]=  251; r[10][14]=   37;
 r[10][15]=   43; r[10][16]=   26; r[10][17]=   20; r[10][18]=    6;
 r[10][19]=  255;

 r[11][12]=    4; r[11][13]=    0; r[11][14]=    0; r[11][15]=    1;
 r[11][16]=    2; r[11][17]=    0; r[11][18]=    5; r[11][19]=    1;

 r[12][13]=    8; r[12][14]=    0; r[12][15]=    1; r[12][16]=   32;
 r[12][17]=    1; r[12][18]=    5; r[12][19]=   89;

 r[13][14]=    0; r[13][15]=   32; r[13][16]=    9; r[13][17]=    2;
 r[13][18]=   54; r[13][19]=   37;

 r[14][15]=    9; r[14][16]=   10; r[14][17]=    0; r[14][18]=    1;
 r[14][19]=    1;

 r[15][16]=  134; r[15][17]=    1; r[15][18]=   22; r[15][19]=   13;

 r[16][17]=    1; r[16][18]=    3; r[16][19]=   48;

 r[17][18]=    2; r[17][19]=   18;

 r[18][19]=    2;



	for (i = 0; i < 20; i++) r[i][i] = 0.0;
	if( isTM )
	{
		for (i = 1; i < 20; i++) for (j = 0; j < i; j++)
		{
			r[j][i] /= 400.0 * freq0_TM[i] * freq0_TM[j];
			r[i][j] = r[j][i];
		}
		for( i=0; i<20; i++ ) freq[i] = freq0_TM[i];
	}
	else
	{
		for (i = 1; i < 20; i++) for (j = 0; j < i; j++)
		{
			r[i][j] /= 400.0 * freq0[i] * freq0[j];
			r[j][i] = r[i][j];
		}
		for( i=0; i<20; i++ ) freq[i] = freq0[i];
	}

	for( i=0; i<26; i++ ) locamino[i] = locamino0[i];
	for( i=0; i<26; i++ ) locgrp[(int)locamino[i]] = locgrp0[i];
	for( i=0; i<20; i++ ) for( j=0; j<20; j++ ) rsr[i][j] = r[i][j];
}
