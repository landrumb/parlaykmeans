#include <iostream>
#include <fstream>
#include <string>

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/io.h"
#include "parlay/primitives.h"
#include "parlay/internal/get_time.h"

/* 
code for parsing the NIH open-i text files to fbin files.

The NIH open-i text files are in the following format:
    - every even line is an image url
    - every odd line is a comma separated list of floats describing a vector embedding
        - for resnet, d = 512
    - each float is rounded to 5 decimal places, has a leading 0, and can be negative

The fbin files are in the following format:
    - the first 4 bytes are the number of vectors in the file
    - the next 4 bytes are the dimension of each vector
    - the rest of the file is the vectors, one after another
        - each vector is a sequence of 4 byte floats
 */

int main(int argc, char* argv[]){
    // get the input and output filenames
    std::string input = argv[1];
    std::string output = argv[2];

    // open the input file
    std::ifstream infile(input);

    // get the number of lines in the file (unfortunate but need to start with number of vectors)
    // the smarter thing to do would be leave space at the beginning, record the number of vectors as we read the file, and then write it to the beginning of the output file at the end
    int32_t num_lines = 0;
    std::string line;
    while (std::getline(infile, line)){
        num_lines++;
    }
    std::cout << "Detected " << num_lines << " vectors" << std::endl;

    // open the input file again
    infile = std::ifstream(input);

    // for completeness, determine the dimension of the first vector (the second line)
    std::getline(infile, line);
    line.clear();
    std::getline(infile, line);
    int32_t d = 1;
    for (size_t i = 0; i < line.size(); i++){
        if (line[i] == ','){
            d++;
        }
    }

    std::cout << "Detected dimension " << d << std::endl;

    infile = std::ifstream(input); // open the input file again

    // open the output file
    std::ofstream outfile(output, std::ios::binary);
    outfile << num_lines;
    outfile << d;

    parlay::internal::timer t;
    t.start();
    // read the input file line by line
    for (int i = 0; i < num_lines; i++){
        if (i % 100000 == 0 && i != 0) {
            double elapsed_time = t.total_time();
            double projected_time = (elapsed_time * num_lines / i) - elapsed_time;

            long elapsed_m = (long) elapsed_time / 60;
            long elapsed_s = (long) elapsed_time % 60;
            long projected_m = (long) projected_time / 60;
            long projected_s = (long) projected_time % 60;
            std::printf("Processed %d/%d vectors in %ldm%lds, projected %ldm%lds remaining\n", i, num_lines, elapsed_m, elapsed_s, projected_m, projected_s);
        }

        line.clear();
        std::getline(infile, line);
        
        // traverse the line, writing each float to the output file
        size_t start = line.find('\t');
        size_t end;
        for (int j = 1; j < d; j++){ // j = 1 because the last float doesn't have a comma after it
            // find the next comma
            end = line.find(',', start);

            // if (i == 0 && j == 1) std::cout << line.substr(start, end - start) << std::endl;

            // get the float
            float f = std::stof(line.substr(start, end - start));

            // write the float to the output file
            outfile.write((char*)&f, sizeof(float));

            // update the start index
            start = end + 1;
        }

        // write the last float to the output file
        float f = std::stof(line.substr(start, line.size() - start));
        outfile.write((char*)&f, sizeof(float));
    }

    // close the files
    infile.close();
    outfile.close();

    return 0;
}
