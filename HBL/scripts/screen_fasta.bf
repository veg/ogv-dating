LoadFunctionLibrary ("libv3/all-terms.bf");
LoadFunctionLibrary ("libv3/IOFunctions.bf");
LoadFunctionLibrary ("libv3/UtilityFunctions.bf");
LoadFunctionLibrary ("libv3/tasks/mpi.bf");
LoadFunctionLibrary ("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("../lib/igscueal.bf");

// SETTINGS
/*
igh_human.settings = {
    "segments" : {
        "0" : {{"V", PATH_TO_CURRENT_BF + "../data/human/IGH.nex","MRCA|01$"}},
        "1" : {{"J", PATH_TO_CURRENT_BF + "../data/human/IGJ.nex","MRCA|01$"}},
        "2" : {{"C", PATH_TO_CURRENT_BF + "../data/human/IGC.nex","MRCA|01$"}}
    },
    "d" : PATH_TO_CURRENT_BF+"../data/human/IGD.seq"
};
*/


reads_path = io.PromptUserForString ("Path to RNA reads file");

igh_human.settings = {"segments" : 
    {
        "0" : {{"HIV", 
                reads_path, 
                "."}}
    }
};



report_headers = {{
/*0*/ "Index",
/*1*/ "Name",
/*2*/ "Best Rearrangement",
/*3*/ "Support",
/*4*/ "Sequence",
/*5*/ "FW1",
/*6*/ "FW1_AA",
/*7*/ "CDR1",
/*8*/ "CDR1_AA",
/*9*/ "FW2",
/*10*/ "FW2_AA",
/*11*/ "CDR2",
/*12*/ "CDR2_AA",
/*13*/ "FW3",
/*14*/ "FW3_AA",
/*15*/ "CDR3",
/*16*/ "CDR3_AA",
/*17*/ "JUNCTION",
/*18*/ "JUNCTION_AA",
/*19*/ "D_ALLELE",
/*20*/ "J",
/*21*/ "J_AA",
/*22*/ "CH",
/*23*/ "CH_AA",
/*24*/ "Mapped Read",
/*25*/ "V-length",
/*26*/ "J-length",
/*27*/ "C-length",
/*28*/ "V-divergence",
/*29*/ "J-divergence",
/*30*/ "C-divergence"
}};

//------------------------------------------------------------------------------------------------

io.DisplayAnalysisBanner({
    "info": "IgSCUEAL (Immunoglobulin subtype classification using evolutionary algorithms)
    uses a phylogenetic placement algorithm to assign a rearrangement to an Ig* (human) read
    V, J, and CH1 are classified this way; D region is assigned based on best nucleotide
    homology.",
    "version": "2.00",
    "reference": "Frost SDW, et al (2015) *Assigning and visualizing germline genes in antibody repertoires* (Phil. Trans. R. Soc. B 2015 370)",
    "authors": "Sergei L Kosakovsky Pond",
    "contact": "spond@temple.edu",
    "requirements": "a FASTA file with (human) Ig reads"
});

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", FALSE);

IgSCUEAL.set_up (igh_human.settings);

SetDialogPrompt ("Select an alignment to screen");
read_set = alignments.ReadNucleotideDataSet ("reads", None);
io.ReportProgressMessageMD("IGSCUEAL", "Load", "Loaded " + read_set["sequences"] + " sequences from `read_set['file']`");

SetDialogPrompt ("Save main screening results to");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN);
IgSCUEAL.result.csv = utility.GetEnvVariable ("LAST_FILE_PATH");
//console.log (IgSCUEAL.result.csv);

SetDialogPrompt ("Save alternative rearrangements to:");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN, "Name\tRearrangement\tSupport\n");
IgSCUEAL.result.rearr.tsv = utility.GetEnvVariable ("LAST_FILE_PATH");


screening.queue = mpi.CreateQueue ({"Headers" : utility.GetListOfLoadedModules ("libv3/"),
                                    "Variables" : {{"igh_human.settings", "NORMALIZE_SEQUENCE_NAMES"}}});



fprintf (IgSCUEAL.result.csv, Join ("\t", report_headers), "\n");

igscueal.start_time = Time (0);
IgSCUEAL.batch_size = Min (25, read_set["sequences"] $ Max (1, utility.GetEnvVariable ("MPI_NODE_COUNT")));
igscueal.finished = 0;

for (seq_id = 0; seq_id < read_set["sequences"]; ) {


    send.indices = {};
    send.ids = {};
    send.sequences = {};

    send.size = Min (IgSCUEAL.batch_size, read_set["sequences"] - seq_id);

    for (batch_id = 0; batch_id < send.size; batch_id += 1) {
        this_read = alignments.GetIthSequence ("reads", seq_id + batch_id);
        send.indices + (seq_id + batch_id + 1);
        send.ids + this_read[terms.id];
        send.sequences + this_read[terms.data.sequence];
    }


    mpi.QueueJob (screening.queue, "igh_human_process_read", {"0" : send.ids,
                                                                "1" : send.indices,
                                                                "2" : send.sequences
                                                              },
                                                              "igh_human_report_read_batch");

    seq_id += send.size;
}

mpi.QueueComplete (screening.queue);
fprintf (IgSCUEAL.result.csv, CLOSE_FILE);
fprintf (IgSCUEAL.result.rearr.tsv, CLOSE_FILE);

// ---------------------------------------------------------------------------------------------------------

lfunction igh_human_report_read_batch (node, batch_results, arguments) {
    //console.log (batch_results);

    size = utility.Array1D (batch_results);
    for (k = 0; k < size; k += 1) {
        igh_human_report_read (node, batch_results[k], arguments);
    }


    ^"igscueal.finished" += size;
    time_so_far = Max (1,Time(0) - ^"igscueal.start_time");

    io.ReportProgressBar("", "\tScreened " + ^"igscueal.finished" + "/" + (^"read_set")["sequences"] + " sequences" +
            ". Elapsed time : " + io.FormatTimeInterval (time_so_far ) +
            ". ETA : " + io.FormatTimeInterval (time_so_far * (((^"read_set")["sequences"]-^"igscueal.finished")/(^"igscueal.finished"))) +
            ". " + Format (^"igscueal.finished"/time_so_far, 5,2) + " sequences/ second" + ".");

}

lfunction igh_human_report_read (node, read_result, arguments) {


    col_count = utility.Array1D (^"report_headers");

    result    = {col_count,1};

    for (i = 0; i < col_count; i+=1) {
        result[i] = "";
    }

    //console.log (read_result);

    result[0] = "" + read_result ["index"];
    result[1] = read_result["ID"];

    //console.log (result);

    if (Type (read_result["assignment"]) != "AssociativeList") {
        result [2] = "Error: alignment failed";
    } else {
        result [2] = (((read_result["assignment"])["REARRANGEMENTS"])["best"])["type"];
        result [3] = "" + (((read_result["assignment"])["REARRANGEMENTS"])["best"])["support"];
        result [4] = (read_result["assignment"])["input-sequence"];
        for (i = 5; i < 24; i+=1) {
            if ((read_result["assignment"])["ANNOTATION"] / (^"report_headers")[i]) {
                result[i] = ((read_result["assignment"])["ANNOTATION"])[(^"report_headers")[i]];
            }
        }
        result [24] = ((read_result["assignment"])["three-way"])[2];
        if ((read_result["assignment"])["SEGMENTS"] / "0") { // has V
            result [25] = "" + (((read_result["assignment"])["SEGMENTS"])["0"])["span"];
        }
        if ((read_result["assignment"])["PHYLO-PLACEMENT"] / "0") { // has V
            result [28] = "" + (((read_result["assignment"])["PHYLO-PLACEMENT"])["0"])["divergence"];
        }
        if ((read_result["assignment"])["SEGMENTS"] / "1") { // has J
            result [26] = "" + (((read_result["assignment"])["SEGMENTS"])["1"])["span"];
        }
        if ((read_result["assignment"])["PHYLO-PLACEMENT"] / "1") { // has V
            result [29] = "" + (((read_result["assignment"])["PHYLO-PLACEMENT"])["1"])["divergence"];
        }
        if ((read_result["assignment"])["SEGMENTS"] / "2") { // has C
            result [27] = "" + (((read_result["assignment"])["SEGMENTS"])["2"])["span"];
        }
        if ((read_result["assignment"])["PHYLO-PLACEMENT"] / "2") { // has V
            result [30] = "" + (((read_result["assignment"])["PHYLO-PLACEMENT"])["2"])["divergence"];
        }
   }

   utility.ForEachPair (((read_result["assignment"])["REARRANGEMENTS"])["credible"], "_id_", "_rearr_",
                         'fprintf (IgSCUEAL.result.rearr.tsv, (`&read_result`)["ID"], "\t", _id_, "\t", _rearr_, "\n")');

   fprintf (^"IgSCUEAL.result.csv", Join ("\t", result), "\n");
}

function igh_human_process_read (sequence_names, indices, sequence_data) {
    if (!igh_human_process_read.setup ) {
        igh_human_process_read.setup = TRUE;
        IgSCUEAL.set_up (igh_human.settings);
    }

    batch_size = utility.Array1D (sequence_names);
    igh_human_process_read.result = {};

    for (igh_human_process_read.k = 0; igh_human_process_read.k < batch_size; igh_human_process_read.k += 1) {
        igh_human_process_read.result + {"index" : indices[igh_human_process_read.k], "ID": sequence_names[igh_human_process_read.k], "assignment" : IgSCUEAL.screen_a_read (sequence_data[igh_human_process_read.k])}
    }

    return igh_human_process_read.result;
}
