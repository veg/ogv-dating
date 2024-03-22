import os, argparse, json, re, sys, csv
from   Bio import SeqIO



if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='scan the directory of annotation files and summarize them'
    )
    
    parser.add_argument(
        '-d', '--directory',
        type=str,
        help='the directory to scan',
        required=True,
    )

    
    parser.add_argument(
        '-j', '--json',
        type=argparse.FileType('r'),
        help='the JSON annotation file',
        required=True,
    )
    
    parser.add_argument(
        '-o', '--output',
        type=argparse.FileType('w'),
        help='the resulting CSV file',
        required=False,
        default = sys.stdout
    )

    parser.add_argument(
        '-c', '--counts',
        action = 'store_true'
    )

    parser.add_argument(
        '-s', '--sequences',
        type=str,
        help='the directory of sequences to scan for QVOA files',
        required=False,
    )

    args = parser.parse_args()
    
    
    prefix_extractor = re.compile ('^([^_]+)_([^_]+_[^_]+)(.+)fa(.*)_classification')
    has_counts = args.counts
    
    summary = {}

    sequence_attributes = {}    
    if args.sequences:
       msa_extractor = re.compile ('combined_nuc')
       for root, dirs, files in os.walk(args.sequences):
            for each_file in files:
                name, ext = os.path.splitext(each_file)
                if len(ext) > 0 and ext in ['.fas']:
                    m = msa_extractor.search (name)
                    if m:
                        with open (os.path.join (root, each_file)) as fh:
                            for seq_record in SeqIO.parse(fh, "fasta"):
                                seq_id   = seq_record.name
                                seq = str (seq_record.seq).upper()
                                sequence_attributes[seq_id] = len ([k for k in seq if k=='N'])/len (seq)
                                
    
    annotations = json.load (args.json)
    
    output = csv.writer (args.output, delimiter = ',')
    headers = ['Subject','Region','Sequence','Combined','Partistic','Clade','Placement','Regression','Partistic_support','Clade_support','Placement_support','Partistic_confidence','Clade_confidence','Root_to_tip_r2']

    if args.sequences:
        headers.append ('N_fraction')
        
    output.writerow (headers)

    #id,partistic,partistic_support,clade,clade_support,regression,placement,placement_support

    def date_value (v, s):
        if v != 'N/A' and v != '':
            return s - float (v)
        return 'N/A'

    for root, dirs, files in os.walk(args.directory):
        for each_file in files:
            name, ext = os.path.splitext(each_file)
            if len(ext) > 0 and ext in ['.csv']:
                m = prefix_extractor.search (name)
                if m:
                   subject = m.group (1)
                   gene = m.group (2)
                   with open (os.path.join (root, each_file)) as fh:
                        reader = csv.reader (fh)
                        try:
                            next (fh)
                        except StopIteration as e:
                            print ("EMPTY file %s" % name, file = sys.stderr)
                            
                        
                        for line in reader:
                           combined = 0.
                           weights  = 0.
                           record = [subject, gene]
                           record.append (line[0])
                           record.append (0)
                           start2art = annotations[subject]['Start2ART'] if subject in annotations else 0
                           #print (line[0], subject, start2art, file = sys.stderr)
                           for c in [2,7,12]:
                                dv = date_value(line[c],start2art)
                                record.append (dv)
                                if dv != 'N/A':
                                    wt = 1/float (line[c+4])
                                    weights  += wt
                                    combined += dv * wt
                                   
                           record.append (start2art - float(line[-1]))
                           for c in [3,8,13]:
                                record.append (line[c])
                           for c in [4,9,14]:
                                record.append (line[c])
                           #record.append (line[1])
                           if combined > 0.:
                            record[3] = combined / weights
                           else:
                            record[3] = 'N/A'
                           if args.sequences:
                                record.append (sequence_attributes[record[2]])
                           output.writerow (record)
                        #print (record)
                

    #with open (args.output, "w") as oh:
    #     json.dump (summary, oh, indent = 1)

    

    sys.exit(0)
