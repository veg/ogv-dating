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

Q = {
    {*,     AC*t,   t,      AT*t}
    {AC*t,  *   ,  CG*t,    CT*t}
    {t,  CG*t   ,  *,       GT*t}
    {AT*t,  CT*t   , GT*t,  *}
};

Model GTR = (Q,freqs,1);

fscanf (PROMPT_FOR_FILE, "Tree", T);

LikelihoodFunction lf = (filter, T);
Optimize (res, lf);


fprintf (PROMPT_FOR_FILE, CLEAR_FILE, Format (T,0,1));


