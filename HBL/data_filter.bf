RequireVersion ("2.3.12");


LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/convenience/regexp.bf");
LoadFunctionLibrary     ("libv3/IOFunctions.bf");
LoadFunctionLibrary     ("lib/igscueal.bf");


utility.SetEnvVariable  ("NORMALIZE_SEQUENCE_NAMES", TRUE);
utility.SetEnvVariable  ("ACCEPT_ROOTED_TREES", TRUE);


filter.analysis_description = {terms.io.info :
                            "Perform initial filtering on HIV-1 sequencing reads.
                            (1) Split RNA and QVOA reads into separate files based on reg exp
                            (2) Normalize RNA sequence reads names
                            (3) Find longest reads, perform corrective alignment to maintain frame
                            (4) Translate to A.A for alignment

                            ",
                            terms.io.version :          "0.1",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "A sequence alignment"
                          };
                          
io.DisplayAnalysisBanner (filter.analysis_description);

filter.code_info = alignments.LoadGeneticCode ("Universal");

KeywordArgument ("input", "A FASTA file with RNA and QVOA sequences (at least some must be full length ORF in one of three frames)");

filter.nuc_data = alignments.ReadNucleotideDataSet ("filter.raw_data", None);

io.ReportProgressMessage ("Data QC", "Loaded `filter.nuc_data[terms.data.sequences]` sequences on `filter.nuc_data[terms.data.sites]` sites from **`filter.nuc_data[terms.data.file]`**");

/* split off QVOA reads */

filter.partitionining_expressions   = {"QVOA" : "(QVOA)|(OGV)"};
filter.partitoned_sequence_names    = regexp.PartitionByRegularExpressions(alignments.GetSequenceNames ("filter.raw_data"), filter.partitionining_expressions); 

//io.CheckAssertion ("Abs(filter.partitoned_sequence_names['(QVOA)|(OGV)'])>0", "There were no sequences marked as QVOA in the input file");

KeywordArgument ("qvoa", "Write QVOA sequences here ", filter.nuc_data[terms.data.file] + "_qvoa.fas");
filter.qvoa_path = io.PromptUserForFilePath ("Save QVOA sequences to");
io.ReportProgressMessage ("Data QC", "Found `Abs(filter.partitoned_sequence_names['(QVOA)|(OGV)'])` QVOA sequences in the file; writing them to **`filter.qvoa_path`**");
alignments.GetSequenceByName ("filter.raw_data", None);
fprintf (filter.qvoa_path, CLEAR_FILE, KEEP_OPEN);

filter.qvoa = {};

utility.ForEach (filter.partitoned_sequence_names["(QVOA)|(OGV)"], "_seq_record_", 
'
    sanitized_name = ((_seq_record_ ^ {{"_xxxx_000WPI"}{""}}) ^ {{"OGV_"}{""}}) + "_QVOA";
    filter.qvoa [sanitized_name] = alignments.StripGaps (alignments.GetSequenceByName ("filter.raw_data", _seq_record_));
    

    fprintf (filter.qvoa_path, ">", sanitized_name , "\n", 
        filter.qvoa [sanitized_name], 
        "\n");
');
fprintf (filter.qvoa_path, CLOSE_FILE);


/* convert sequence names from RNA reads */

filter.RNA_reads     = {};
filter.lookup_cache  = {};
filter.seq_count     = 1;
filter.sequence_info = {};
filter.longest_seq   = "";
filter.longest_seq_L = 0;
filter.clean_seqs   = {};
filter.frameshifted = {};
filter.unique       = {};


KeywordArgument ("protein", "Translated RNA sequences", filter.nuc_data[terms.data.file] + "_protein.fas");
filter.protein_path = io.PromptUserForFilePath ("Save translated RNA sequences file to");

KeywordArgument ("rna", "Reduced RNA sequences", filter.nuc_data[terms.data.file] + "_nuc.fas");
filter.nuc_path = io.PromptUserForFilePath ("Save reduced RNA sequences file to");


KeywordArgument ("combined-protein", "Translated RNA+QVOA sequences", filter.nuc_data[terms.data.file] + "_combined_protein.fas");
filter.combined_protein_path = io.PromptUserForFilePath ("Save translated RNA+QVOA sequences file to");

KeywordArgument ("combined-rna", "Reduced RNA+QVOA sequences", filter.nuc_data[terms.data.file] + "_combined_nuc.fas");
filter.combined_nuc_path = io.PromptUserForFilePath ("Save reduced RNA sequences file to");


io.ReportProgressMessage ("Data QC", "Will write unaligned unique protein sequences for MSA to **`filter.protein_path`**, and the corresponding nucleotide sequences to **`filter.nuc_path`**");
io.ReportProgressMessage ("Data QC", "Will write unaligned unique protein sequences + QVOA for MSA to **`filter.combined_protein_path`**, and the corresponding nucleotide sequences to **`filter.combined_nuc_path`**");

fprintf (filter.protein_path, CLEAR_FILE, KEEP_OPEN);
fprintf (filter.nuc_path, CLEAR_FILE, KEEP_OPEN);

fprintf (filter.combined_protein_path, CLEAR_FILE, KEEP_OPEN);
fprintf (filter.combined_nuc_path, CLEAR_FILE, KEEP_OPEN);

filter.sequences_with_copies = {};

utility.ForEach (filter.partitoned_sequence_names[""], "_seq_record_", 
'
    io.ReportProgressBar ("filter","Filtering RNA sequence " + filter.seq_count);
    filter.read_to_check = alignments.Strip (alignments.StripGaps (alignments.GetSequenceByName ("filter.raw_data", _seq_record_)));
    
    if (filter.unique / filter.read_to_check == FALSE) {    
        
        filter.RNA_reads[_seq_record_] = filter.read_to_check;
        
        filter.sequences_with_copies [filter.read_to_check] = {"0" : _seq_record_};
        filter.unique [filter.read_to_check] = _seq_record_;
    
        filter.sequence_info[_seq_record_] = alignments.TranslateCodonsToAminoAcidsWithAmbigsAllFrames (filter.RNA_reads[_seq_record_],
                               filter.code_info, filter.lookup_cache);
                               
                                                                                                 
        for (frame = 0; frame < 3; frame += 1) {

            if (((filter.sequence_info[_seq_record_])[frame])[terms.stop_codons] == 0) {
                if (((filter.sequence_info[_seq_record_])[frame])[terms.sense_codons] > filter.longest_seq_L) {
                    filter.longest_seq_L = ((filter.sequence_info[_seq_record_])[frame])[terms.sense_codons];
                    filter.longest_seq = ( filter.RNA_reads[_seq_record_]) [frame][ Abs(filter.RNA_reads[_seq_record_]) - 1];
                    filter.longest_seq_NL = Abs(filter.longest_seq);
                    if ( filter.longest_seq_NL % 3) {
                        filter.longest_seq_NL = filter.longest_seq_NL$3*3;
                        filter.longest_seq = filter.longest_seq[0][filter.longest_seq_NL-1];
                    }
                }
                filter.clean_seqs [_seq_record_] = ((filter.sequence_info[_seq_record_])[frame])[terms.data.sequence];
                break;
            }
        }   
                           
        if (frame == 3) {
            filter.frameshifted [_seq_record_] = 1;
        }  
    } else {
        filter.sequences_with_copies [filter.read_to_check] + _seq_record_;
    }
                           
    filter.seq_count += 1;
');


io.ClearProgressBar ();
io.ReportProgressMessage ("Data QC", "Found `Abs(filter.clean_seqs)` unique RNA sequences that were in frame");
io.CheckAssertion ("filter.longest_seq_L>0", "There were no RNA sequences that were in frame and had no stop codons");


filter.ref_seq = {"REF" : {'stripped' : filter.longest_seq}};
filter.options = IgSCUEAL.define_alignment_settings (filter.code_info);   
filter.options["code"] = filter.code_info;


if (Abs(filter.frameshifted)) {
    io.ReportProgressMessage ("Data QC", "Correcting frame-shifts in the remaining reads");
    filter.seq_count = 1;
 
    utility.ForEachPair (filter.frameshifted, "_sequence_", "_value_",
    '
        io.ReportProgressBar ("filter","Processing sequence " + filter.seq_count);
        filter.cleaned = IgSCUEAL.align_sequence_to_reference_set (filter.RNA_reads[_sequence_], filter.ref_seq, filter.options);
        
        //if (_sequence_ == "CAP206_159wpi_ENV_2_all_aligned_loops_trimmed_98_2") {
        //    console.log (filter.cleaned);
        //}
 
        filtered.aa_seq = alignments.StripGaps(filter.cleaned["AA"]);
        filtered.na_seq = IgSCUEAL.strip_in_frame_indels(filter.cleaned["QRY"]);
        
        (filter.sequences_with_copies[filter.RNA_reads[_sequence_]])["_write_to_file"][""];
        filter.seq_count += 1;
    ');

    io.ClearProgressBar ();
    
} 

io.ReportProgressMessage ("Data QC", "Checking for frame-preserving indels in other reads");
filter.seq_count = 1;

utility.ForEachPair (filter.clean_seqs, "_sequence_", "_value_",
'
    io.ReportProgressBar ("filter","Processing sequence " + filter.seq_count);
    filter.cleaned = IgSCUEAL.align_sequence_to_reference_set (filter.RNA_reads[_sequence_], filter.ref_seq, filter.options);

    //if (_sequence_ == "CAP206_159wpi_ENV_2_all_aligned_loops_trimmed_98_2") {
    //    console.log (filter.cleaned);
    //}

    filtered.aa_seq = alignments.StripGaps(filter.cleaned["AA"]);
    filtered.na_seq = IgSCUEAL.strip_in_frame_indels(filter.cleaned["QRY"]);
    (filter.sequences_with_copies[filter.RNA_reads[_sequence_]])["_write_to_file"][""];
    filter.seq_count += 1;
');

io.ClearProgressBar ();


io.ReportProgressMessage ("Data QC", "Performing QA on QVOA reads");

utility.ForEachPair (filter.qvoa, "_sequence_", "_value_",
'
    filter.cleaned = IgSCUEAL.align_sequence_to_reference_set (_value_, filter.ref_seq, filter.options);
    filtered.aa_seq = alignments.StripGaps(filter.cleaned["AA"]);
    filtered.na_seq = IgSCUEAL.strip_in_frame_indels(filter.cleaned["QRY"]);
    
    fprintf (filter.combined_protein_path, ">", _sequence_, "\n",  filtered.aa_seq, "\n");
    fprintf (filter.combined_nuc_path, ">", _sequence_, "\n", filtered.na_seq , "\n");
');
    


fprintf (filter.protein_path, CLOSE_FILE);
fprintf (filter.nuc_path, CLOSE_FILE);
fprintf (filter.combined_protein_path, CLOSE_FILE);
fprintf (filter.combined_nuc_path, CLOSE_FILE);

io.ReportProgressMessage ("Next steps", "Please run **`filter.protein_path`** through an MSA program, and then run data_filter-2.bf on the output and **`filter.nuc_path`** to recover the nucleotide MSA");


//console.log (filter.clean_seqs);

function _write_to_file (key, value) {
    fprintf (filter.protein_path, ">", value, "\n",  filtered.aa_seq, "\n");
    fprintf (filter.nuc_path, ">", value, "\n", filtered.na_seq , "\n");
    fprintf (filter.combined_protein_path, ">", value, "\n",  filtered.aa_seq, "\n");
    fprintf (filter.combined_nuc_path, ">", value, "\n", filtered.na_seq , "\n");
 }




