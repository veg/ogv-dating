#!/usr/bin/env node

const fs = require('fs');
const commander = require('commander');
const pd = require("pretty-data").pd;
const csvParse = require("d3").csvParse;
const _ = require('underscore');
//const stats = require("stats-lite");
const csv   = require("fast-csv");
const phylotree = require('phylotree');

commander.description('Read a Newick tree, partition leaves into IN/OUT groups and compute the closest OUT neighbor for each leaf in the IN').
  version('0.0.1', '-v, --version').
  option(
    '-n --newick <newick>',
    'Input newick file'
  )
  .option(
    '-p --placement [placement]',
    'CSV file with phylogenetic placement information')
  .option(
    '-b --bootstrap [value]',
    'Minimim bootstrap support [0-1] to consider a node supported',
    0.9
  ).option(
    '-r --regex [regex] [otherRegexes...]',
    'Name,Regex pairs to match (comma delimited), i.e. Blood,BL[0-9]+'
  ).option(
    '-d --date [date]',
    'A regexp with at least one subexpression used to extract date information from sequence names (default ([0-9]+)WPI)',
    '([0-9]+)WPI'
  ).parse(process.argv);

if (!commander.regex){
  throw 'ERROR: Regex (-r) input is required';
}

if (!commander.placement){
  throw 'ERROR: Placement (-p) input is required';
}

const date_regexp = new RegExp(commander.date);
const date_extract = (n) => 
	{
	    const m = date_regexp.exec(n);
        if (m) {
            return m[1];
        }
        return null;
    };


var placement_stream  = fs.createReadStream(commander.placement);
let placement_support = {};
let partistic_support  = {};
let clade_support = {};


var csvStream = csv.parse({delimiter : "\t", headers : true})
    .on("data", (data) => {
            const id = data.Name;
    
            if (!(id in placement_support)) {
                placement_support[id] = {};
            }
   
            const date_tag = date_extract (data.Rearrangement);
            if (date_tag) {
                placement_support [id][date_tag] = (placement_support [id][date_tag] ? placement_support [id][date_tag] : 0) + (+data.Support);
            }
        })
    .on("end", function(){
        
         fs.readFile(commander.newick, (err, newick_data) => {
         const tree = new phylotree.phylotree(newick_data.toString());
  
          var  partitions = [commander.regex].concat (commander.args).map(regex_str => {
                const split = regex_str.split(',');
                return {
                  name: split[0],
                  regex: new RegExp(split[1]),
                  count: 0,
                  node_indices : new Object
                };
              });
      
    
            var group_by_index = {};
    
            phylotree.extract_dates (tree, (n) => {
                return date_extract (n.data.name);
            }, (d) => {
                // convert WPI into years post-infection
                var date = new Date (2000 + (+d), 0, 1);
                return date;
            }
            );
    
     
    
            tree.traverse_and_compute(function(node) {
              partitions.forEach((c,i) => {
                if(node.data.name.match(c.regex)) {
                  if (node.decimal_date_value) {
                    node.decimal_date_value -= 2000;
                  }
                }
              });
            });
    
            const rtt = phylotree.fit_root_to_tip (tree);
            tree.reroot (rtt.root);
            phylotree.pairwise_distances (tree);
            phylotree.root_to_tip (tree);
               
            tree.traverse_and_compute(function(node) {
              partitions.forEach((c,i) => {
                if(node.data.name.match(c.regex)) {
                  node.data.trait = c.name;
                  group_by_index[node.cot_leaf_index]  = [i, node];
                  c.node_indices [node.cot_leaf_index] = 1; 
                  c.count++;
                }
              });
            });
    
            const bs_cut = +commander.bootstrap; 
    
            const summarize_set = (s) => {
                if (s) {
                    const size = _.reduce (s, (i,n) => i + n, 0);
                    const s2 = _.map (s, (i,k) => [k,i/size]);
                    return _.max (s2, (s) => s[1]);
                } 
                return s;
            };
            
            
     
            var csvStream = csv
                .createWriteStream({headers: true})
                .transform(function (row) { 
                        return  {
                           id: row.id, 
                           r2: row.r2,
                           partistic: (row.partistic ? row.partistic[0] : 'N/A'),
                           partistic_support: row.partistic ? row.partistic[1] : 'N/A',
                           partistic_tree_confidence: row.partistic && (row.partistic[0] in partistic_support) ? partistic_support[row.partistic[0]][1] / partistic_support[row.partistic[0]][0] : 'N/A',
                           clade: row.clade ? row.clade[0] : 'N/A',
                           clade_support: row.clade ? row.clade[1] : 'N/A',
                           clade_tree_cofidence: row.clade && (row.clade[0] in  clade_support)? clade_support[row.clade[0]][1] / clade_support[row.clade[0]][0] : 'N/A',
                           regression: row.regression,
                           placement: row.placement[0],
                           placement_support: row.placement[1]
                        };
                    }
                );
        
            csvStream.pipe (process.stdout);
            
            let qvoa_reports      = [];

            tree.traverse_and_compute(function(node) {
              if (tree.is_leafnode (node)) {
                var neighbor_nodes = _.map (node.cot_path_to_leaves_above, (v,k) => [v + node.cot_computed_length,+k]).filter (v=>group_by_index[v[1]][0] != 0).sort ((a,b) => a[0] - b[0]);
                const min_distance   = neighbor_nodes[0][0]; 
        
                /*const mean_distance = stats.mean (neighbor_nodes.map (v=>v[0]));
                const sd_distance = stats.stdev (neighbor_nodes.map (v=>v[0]));
                const z_distance = neighbor_nodes.map (v => [Math.abs (v[0]-mean_distance) / sd_distance,v[1]]);
                console.log (z_distance);
                */

                //if (node.data.trait == partitions[0].name)
                //    console.log ("\n", node.data.name, neighbor_nodes, min_distance);
        
                neighbor_nodes = _.filter (neighbor_nodes, (v) => v[0] <= 2*min_distance);
               
                let nn_count = neighbor_nodes.length;
                
                 neighbor_nodes = _.countBy (neighbor_nodes.map (v => [v[0],group_by_index[v[1]][1].decimal_date_value]), (v) => v[1]);
                
        
                
                let supported_parent = node.parent;
                while (supported_parent.parent && (+supported_parent.data.bootstrap_values) < bs_cut) {
                    supported_parent = supported_parent.parent;
                }
                var spp = _.map (_.filter (supported_parent.descendants(), (n) => tree.is_leafnode (n) && n.data.trait !=  partitions[0].name), (n) => n.decimal_date_value);
                while (supported_parent.parent && _.size (spp) == 0 && (+supported_parent.data.bootstrap_values) < bs_cut) {
                    supported_parent = supported_parent.parent;
                    spp = _.map (_.filter (supported_parent.descendants(), (n) => tree.is_leafnode (n) && n.data.trait !=  partitions[0].name), (n) => n.decimal_date_value);
                }
                
                if (node.decimal_date_value) {
                    if (!(node.decimal_date_value in partistic_support)) {
                        partistic_support [node.decimal_date_value] = [0,0]
                    } 
                    partistic_support [node.decimal_date_value][0] ++;
                    partistic_support [node.decimal_date_value][1] += summarize_set(neighbor_nodes)[0] == node.decimal_date_value;
                    
                    if (supported_parent.parent) {
                        if (!(node.decimal_date_value in clade_support)) {
                            clade_support [node.decimal_date_value] = [0,0]
                        } 
                        clade_support [node.decimal_date_value][0] ++;
                        clade_support [node.decimal_date_value][1] += summarize_set(_.countBy (spp))[0] == node.decimal_date_value;
                    }
                }
                
                if (node.data.trait == partitions[0].name) {            
                    let name_wihout_counts = node.data.name.split ('_');
                    name_wihout_counts.splice (-1,1);
                    name_wihout_counts = name_wihout_counts.join ('_');
        
                    qvoa_reports.push ({"id" : node.data.name,
                        "r2" : rtt.fit.r2,
                        partistic : summarize_set(neighbor_nodes),
                        regression : rtt.fit.slope * node.root_to_tip + rtt.fit.intercept,
                        clade : summarize_set(supported_parent.parent ? _.countBy (spp) : null),
                        placement: _.max(_.map (placement_support[name_wihout_counts], (i,k) => [k,i]), (s) => s[1])
                        
                    });
                }
             }
            });
    
    
            qvoa_reports.forEach (r => csvStream.write (r));
            csvStream.end();  
            
            
 });

});

placement_stream.pipe (csvStream);
     
