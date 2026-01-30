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
#include <cmath>

#include "util.hpp"

void usage(char *progname);

int main(int argc, char** argv) {
    
	std::string in_filename = "";
	std::string out_filename;

	bool verbose = false;


	// --------------------- START ------------------------------------------------------------------
	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "hvi:o:")) != -1) {
		switch(c)  {
			case 'h':
				usage(argv[0]);
			case 'v':
				verbose = true; break;
			case 'o':
				out_filename = optarg; break;
			case 'i':
				in_filename = optarg; break;
			default:
				usage(argv[0]);
		}
	}

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

    // write header
    *out_stream << "readID" << "\t" << "taxID" << "\t" << "score" << "\t" << "queryLength" << "\n";

	// read file and count reads
	std::string line;
    size_t start;
    size_t end;
	while(getline(in_file,line)) {
		if(line.length() == 0) { continue; }

		start = 0;
		end = line.find('\t',start);
        char classified_code;
		try {
			classified_code = line.substr(start,end-start)[0];
		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Error: Found bad classified/unclassifed indicator in line: " << line << std::endl;
			continue;
		}

		start = end+1;
		end = line.find('\t',start);
        std::string readID = line.substr(start,end-start);
		
        start = end+1;
		end = line.find('\t',start);
		uint64_t taxID;
		try {
			taxID = stoul(line.substr(start,end-start));
		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Error: Found bad taxon id in line: " << line << std::endl;
			continue;
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Error: Found bad taxon id (out of range error) in line: " << line << std::endl;
			continue;
		}

		start = end+1;
		end = line.find('\t',start);
		uint64_t score;
		try {
			score = stoul(line.substr(start,end-start)) * 3; // multiply by 3 for NT length
		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Error: Found bad score in line: " << line << std::endl;
			continue;
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Error: Found bad score (out of range error) in line: " << line << std::endl;
			continue;
		}

		start = end+1;
		end = line.find('\t',start);
		double queryLength;
		try {
			queryLength = stod(line.substr(start,end-start)) * 3; // multiply by 3 for NT length
			queryLength = static_cast<int>(std::round(queryLength));

		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Error: Found bad query length in line: " << line << std::endl;
			continue;
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Error: Found bad query length (out of range error) in line: " << line << std::endl;
			continue;
		}

        *out_stream << readID << "\t" << taxID << "\t" << score << "\t" << queryLength << "\n";
        
	}  // end while getline

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
	fprintf(stderr, "Usage:\n   %s -i kaiju.out -o kaiju-recentrifugeformat.out\n", progname);
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -i FILENAME   Name of input file (kaiju output)\n");
	fprintf(stderr, "   -o FILENAME   Name of output file. If not specified, output will be printed to STDOUT.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -v            Enable verbose output.\n");
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
}