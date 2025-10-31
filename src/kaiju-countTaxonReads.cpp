/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. 
 * Added by Sean Golez */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <string>
#include <list>
#include <utility>
#include <stdexcept>
#include <deque>

#include "util.hpp"

void usage(char *progname);

int main(int argc, char** argv) {

	std::unordered_map<uint64_t,uint64_t> nodes;
	std::unordered_map<uint64_t, std::string> node2name;
	std::unordered_map<uint64_t, std::string> node2rank;
	std::unordered_map<uint64_t, std::string> node2path;

	std::string nodes_filename = "";
	std::string names_filename = "";
	std::string in_filename = "";
	std::string out_filename;

	bool verbose = false;

    std::unordered_map<uint64_t,uint64_t> taxid_counts;


	// --------------------- START ------------------------------------------------------------------
	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "hvn:t:i:o:")) != -1) {
		switch(c)  {
			case 'h':
				usage(argv[0]);
			case 'v':
				verbose = true; break;
			case 'o':
				out_filename = optarg; break;
			case 'n':
				names_filename = optarg; break;
			case 't':
				nodes_filename = optarg; break;
			case 'i':
				in_filename = optarg; break;
			default:
				usage(argv[0]);
		}
	}
	if(names_filename.length() == 0) { error("Please specify the location of the names.dmp file with the -n option."); usage(argv[0]); }
	if(nodes_filename.length() == 0) { error("Please specify the location of the nodes.dmp file, using the -t option."); usage(argv[0]); }
	if(in_filename.length() == 0) { error("Please specify the location of the input file, using the -i option."); usage(argv[0]); }

	/* read nodes.dmp */
	std::ifstream nodes_file;
	nodes_file.open(nodes_filename);
	if(!nodes_file.is_open()) { std::cerr << "Error: Could not open file " << nodes_filename << std::endl; usage(argv[0]); }
	if(verbose) std::cerr << "Reading taxonomic tree from file " << nodes_filename << std::endl;
	parseNodesDmpWithRank(nodes,node2rank,nodes_file);
	nodes_file.close();

	/* read names.dmp */
	std::ifstream names_file;
	names_file.open(names_filename);
	if(!names_file.is_open()) { std::cerr << "Error: Could not open file " << names_filename << std::endl; usage(argv[0]); }
	if(verbose) std::cerr << "Reading taxon names from file " << names_filename << std::endl;
	parseNamesDmp(node2name,names_file);
	names_file.close();

	std::ifstream in_file;
	in_file.open(in_filename);
	if(!in_file.is_open()) {  std::cerr << "Could not open file " << in_filename << std::endl; exit(EXIT_FAILURE); }

	std::ostream * out_stream;
	if(out_filename.length()>0) {
		if(verbose) std::cerr << "Output file: " << out_filename << std::endl;
		std::ofstream * output_file = new std::ofstream();
		output_file->open(out_filename);
		if(!output_file->is_open()) {  std::cerr << "Could not open file " << out_filename << " for writing" << std::endl; exit(EXIT_FAILURE); }
		out_stream = output_file;
	}
	else {
		out_stream = &std::cout;
	}

	if(verbose) std::cerr << "Processing " << in_filename <<"..." << "\n";

	// read file and count reads
	std::string line;
	while(getline(in_file,line)) {
		if(line.length() == 0) { continue; }

		size_t found = line.find('\t');
		found = line.find('\t',found+1);
		size_t end = line.find_first_not_of("0123456789",found+1);
		uint64_t taxonid;
		try {
			taxonid = stoul(line.substr(found,end-found));
		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Error: Found bad taxon id in line: " << line << std::endl;
			continue;
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Error: Found bad taxon id (out of range error) in line: " << line << std::endl;
			continue;
		}

		if(taxonid != 0 && nodes.count(taxonid)==0) {
			std::cerr << "Warning: Taxon ID " << taxonid << " in output file is not contained in taxonomic tree file "<< nodes_filename << ".\n";
			continue;
		}
		if(taxonid != 0 && node2name.count(taxonid)==0) {
			std::cerr << "Warning: Taxon ID " << taxonid << " in output file is not found in file "<< names_filename << ".\n";
			continue;
		}

        taxid_counts[taxonid] += 1;
	}  // end while getline

	// sort by count descending
	std::vector<std::pair<uint64_t,uint64_t>> sorted_counts(taxid_counts.begin(), taxid_counts.end());
	std::sort(sorted_counts.begin(), sorted_counts.end(), 
		[](const std::pair<uint64_t, int>& a, const std::pair<uint64_t, int>& b) {
			if( a.second == b.second ) {
				return a.first < b.first;
			}
			return a.second > b.second;
		});

	// output to file
	for(auto& count_pair : sorted_counts) {
		*out_stream << count_pair.first << '\t' << count_pair.second << "\n";
	}

	if(in_file.is_open()) {
		in_file.close();
	}
	out_stream->flush();
	if(out_filename.length()>0) {
		((std::ofstream*)out_stream)->close();
		delete ((std::ofstream*)out_stream);
	}

	return 0;

}

void usage(char *progname) {
	print_usage_header();
	fprintf(stderr, "Usage:\n   %s -t nodes.dmp -n names.dmp -i kaiju.out -o kaiju-counts.out\n", progname);
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -i FILENAME   Name of input file\n");
	fprintf(stderr, "   -o FILENAME   Name of output file. If not specified, output will be printed to STDOUT.\n");
	fprintf(stderr, "   -t FILENAME   Name of nodes.dmp file\n");
	fprintf(stderr, "   -n FILENAME   Name of names.dmp file.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -v            Enable verbose output.\n");
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
}