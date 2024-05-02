RequireVersion("2.5.55");

NORMALIZE_SEQUENCE_NAMES = TRUE;

DataSet msa = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filter = CreateFilter (msa, 1);


global AC = 1;
global AT = 1;
global CG = 1;
global CT = 1;
global GT = 1;

HarvestFrequencies (freqs, filter, 1, 1, 1);

global alpha = .5;
	alpha:>0.01;alpha:<100;
	category c = (4, EQUAL, MEAN, 
					GammaDist(_x_,alpha,alpha), 
					CGammaDist(_x_,alpha,alpha), 
					0 , 
			  	    1e25,
			  	    CGammaDist(_x_,alpha+1,alpha)
			  	 );

Q = {
    {*,     c*AC*t,   c*t,      c*AT*t}
    {c*AC*t,  *   ,  c*CG*t,    c*CT*t}
    {c*t,  c*CG*t   ,  *,       c*GT*t}
    {c*AT*t,  c*CT*t   , c*GT*t,  *}
};

Model GTR = (Q,freqs,1);

fscanf (PROMPT_FOR_FILE, "Tree", T);

LikelihoodFunction lf = (filter, T);

Optimize (res, lf);

fprintf (PROMPT_FOR_FILE, CLEAR_FILE, Format (T,0,1));


