RequireVersion ("2.4.0");


LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/convenience/regexp.bf");
LoadFunctionLibrary     ("libv3/IOFunctions.bf");

filter.analysis_description = {terms.io.info :
                            "
                            Map a protein MSA back onto nucleotide sequences
                            ",
                            terms.io.version :          "0.1",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "A protein MSA and the corresponding nucleotide alignment"
                          };
                          
io.DisplayAnalysisBanner (filter.analysis_description);

SetDialogPrompt ("Protein MSA");

filter.code_info = alignments.LoadGeneticCode ("Universal");


console.log ("Loading protein MSA");
filter.protein = alignments.ReadProteinDataSet ("filter.protein_data", None);

console.log ("Loading nucleotide unaligned sequences");
filter.nuc = alignments.ReadNucleotideDataSet ("filter.nuc_data", None);

alignments.GetSequenceByName ("filter.nuc_data", None);

filter.path = io.PromptUserForFilePath ("Save nucleotide in-frame MSA to ");
fprintf (filter.path, CLEAR_FILE, KEEP_OPEN);

filter.tags = {};

for (filter.id = 0; filter.id < filter.protein[terms.data.sequences]; filter.id += 1) {
    //filter.nuc_seq  = alignments.GetSequenceByName ("filter.nuc_data", _sequence_);
    filter.prot_seq = alignments.GetIthSequence ("filter.protein_data", filter.id);
    filter.nuc_seq = alignments.GetSequenceByName ("filter.nuc_data", filter.prot_seq[terms.id]);    
    
    if (/*regexp.Find(filter.prot_seq[terms.id], "QVOA")*/ TRUE) {
        filter.seq_tag = filter.prot_seq[terms.id];
    } else {    
        filter.seq_tags = regexp.Split(filter.prot_seq[terms.id], "_");
        filter.seq_tag = filter.seq_tags[0] + "_";
        for (k = 1; k < utility.Array1D (filter.seq_tags); k+=1){
            if ((filter.seq_tags[k] $ "WPI")[0] >= 0) {
                filter.seq_tag  += filter.seq_tags[k];
                break;
            }
        }
    
        filter.tags [filter.seq_tag] += 1;
        filter.seq_tag +=  "_" + filter.tags [filter.seq_tag];
    }
    //console.log (filter.prot_seq[terms.id]);
       
    fprintf (filter.path, ">", filter.seq_tag, "\n", alignment.MapCodonsToAA(filter.nuc_seq, filter.prot_seq[terms.data.sequence] , 1, filter.code_info[terms.code.mapping]), "\n");
}

fprintf (filter.path, CLOSE_FILE);

filter.aligned = alignments.ReadNucleotideDataSet ("filter.dataset", filter.path);
DataSetFilter filter.all = CreateFilter (filter.dataset, 1);

filter.names = alignments.GetSequenceNames ("filter.all");

filter.splits = regexp.PartitionByRegularExpressionsWithSub (filter.names, {"0" : "QVOA|OGV", "1" : "([0-9]+WPI)"});

fprintf (filter.path, CLEAR_FILE);

utility.SetEnvVariable ("DATA_FILE_PRINT_FORMAT",9);


utility.ForEach (filter.splits, "_sequences_", 'filter.compress ( _sequences_)');

filter.count = 0;


function filter.compress (s) {
    if (utility.Array1D (s)) {
        filter.selector  = utility.SwapKeysAndValues (s);
        DataSetFilter filter.subfilter = CreateFilter (filter.dataset, 1, "", filter.selector / filter.names[speciesIndex]);
        io.ReportProgressMessage ("UNIQUE SEQUENCES", "Retained `alignments.CompressDuplicateSequences ('filter.subfilter','filter.subfilter.compressed', TRUE)` unique  sequences");
        if (filter.count > 0) {
            fprintf (filter.path, "\n", filter.subfilter.compressed);
        } else {
            fprintf (filter.path, filter.subfilter.compressed);
        }
        filter.count += 1;
    }
    
}

