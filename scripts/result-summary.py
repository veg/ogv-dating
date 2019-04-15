import os, argparse, json, re, sys, csv


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

    args = parser.parse_args()
    
    prefix_extractor = re.compile ('^([^_]+)_([^_]+_[^_]+)(.+)fasta_classification')
    
    summary = {}
    
    annotations = json.load (args.json)
    
    output = csv.writer (args.output, delimiter = ',')
    output.writerow (['Subject','Region','Sequence','Copy_number','Partistic','Clade','Placement','Regression','Partistic_support','Clade_support','Placement_support','Partistic_confidence','Clade_confidence','Root_to_tip_r2'])

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
                        next (fh)
                        for line in reader:
                           record = [subject, gene]
                           record.append (line[0])
                           record.append (int(line[0].split ('_')[-1]))
                           for c in [2,5,9]:
                                record.append (date_value(line[c], annotations[subject]['Start2ART']))
                           record.append (annotations[subject]['Start2ART'] - float(line[8]))
                           for c in [3,6,10]:
                                record.append (line[c])
                           for c in [4,7]:
                                record.append (line[c])
                           record.append (line[1])
                           output.writerow (record)
                        #print (record)

    #with open (args.output, "w") as oh:
    #     json.dump (summary, oh, indent = 1)

    

    sys.exit(0)
