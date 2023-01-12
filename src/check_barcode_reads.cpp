#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "Trie.h"
#include "ResizeArray.h"
using namespace std;


#define BLOCKSIZE 10000000
#define TERMINATOR '@'
#define MAX_LL 8192
struct a_barcode {
    string sequence;
    int original_pos;
};


a_barcode **barcodes;
//trie_node *barcode_single_trie_head;

////////////// Functions for reading barcodes from input files, and storing in struct arrays
int Get_Lines_In_File(string filename)
{
    /*
  Iterates over the file and counts the number of lines in the file.
  */
    fstream file;
    file.open(filename, ios::in);
    
    string line;
    int N = 0;

    while (getline(file, line)) {
        N++;
    }
    file.close();
    return N;
}

int Read_In_Barcodes(string filename)
{
    /*
  Read in the barcodes from a given textfile, which needs to contain newline char seperated barcode entries.
  All barcodes must be the same length.
  Barcodes get read into an a_barcode struct, along with other information,
  and the struct is then stored in the barcodes array.
  filename: a string of the file name to read from
  return: the number of barcodes read in.
  post-condition: barcodes contains structs of all the barcodes in the given file.
  */
    fstream file;
    //size_t len = 1000;
    file.open(filename, ios::in);

    int num_barcode = Get_Lines_In_File(filename);
    barcodes = new a_barcode *[num_barcode];
    string line;
    a_barcode *new_barcode;                          /// create a new_barcode struct variable
    

    int count = 0;
    //char *token;
    while (getline(file, line))
    {
        new_barcode = new a_barcode; // new makes a pointer
        new_barcode->sequence = line;
        new_barcode->original_pos = count;
        
        barcodes[count] = new_barcode;
        count++;
    }
    file.close();

    return num_barcode;
}


void
Build_Trie_Barcodes(Trie *this_trie, int num_barcode)
{
    /*
  Build a trie using the barcodes array, with parameters allowing support for building paired read, or dual indexing barcodes.
  For every barcode in the barcodes array, add it to the trie by:
    Starting at the head node, add a node for the current char if none exists, or follow to the relevant node
    Repeat for every char in the barcode, finally inserting a terminator character at the end of the trie path
  
  is_paired: boolean value indicating if we should create a paired read trie, storing all sequences in barcodes->sequenceRev
  is_dualindex: boolean value indicating if we should insert barcodes->sequence2 into the created trie.
    Both of these booleans cannot be true. is_paired will be prioritised over is_dualindex
  return: a pointer to the head of the trie, which contains the empty character
  */
    //trie_node *head = Initialise_Node('\0');
    //trie_node *current_node;
    string cur_seq;

    // For every barcode in the barcodes array, add it to the trie by iteratively following the trie
    // through each character in the barcode, adding the characters that don't exist yet.
    // The original position of the barcodes are recorded at the terminator character
    for (int bc_i = 0; bc_i < num_barcode; bc_i++)
    {
        cur_seq = barcodes[bc_i]->sequence;

        this_trie->Add_String(cur_seq, barcodes[bc_i]->original_pos, bc_i);
    }
}


// barcode sorting
int
Base_to_Int(char* base) {
  /*
  Determine the position of a base in the array used for count sort
  base: the base to convert to an int
  return: an int in the range 0-4 inclusive
  */
  switch (*base) {
        case 'A':
          return 1;
          // positions[1]++ usingthis above returns the pre-incremented value of positions
        case 'C':
          return 2;
        case 'G':
          return 3;
        case 'T':
          return 4;
        case TERMINATOR:
        default:
          return 0;
      }
}

void
Count_Sort_Barcodes(long index, a_barcode** input_barcodes, a_barcode** sorted_barcodes,
                    int num_barcode, int barcode_length) {
  /* 
  Implements Count Sort, which stable-y sorts the input_barcodes, making use of the intermediate
  sorted_barcodes array given, to store the sorted barcodes as we find their positions.
  
  index: the index of the barcode to sort based on
  input_barcodes: a pointer to the unsorted array of barcodes
  sorted_barcodes: a pointer to a allocated array of a_barcode structs, which may or may not contain pointers.
                  These will be overridden
  */
  long counts[5]; // counts for empty, A, C, G, T
  long positions[5]; // positions for empty A, C, G, T
  int arr_pos;
  // intialise the counts array to 0
  for (arr_pos = 0; arr_pos < 5; arr_pos++) {
    counts[arr_pos] = 0;
  }

  // Count each appearance of the bases in the array
  char base;
  for (arr_pos = 0; arr_pos < num_barcode; arr_pos++) {
      // count the number of occurances of each base
      base = input_barcodes[arr_pos]->sequence[index];
      counts[Base_to_Int(&base)]++;
  }

  // positions holds [whitespace, A, C, G, T]
  positions[0] = 0;
  // determine the positions of the first appearance of each base in the array.
  // position[i] = position[i-1] + count[i-1]
  for (arr_pos = 1; arr_pos < 5; arr_pos++) {
    positions[arr_pos] = positions[arr_pos-1] + counts[arr_pos-1];
  }

  // construct the sorted barcodes array from our stored positions
  for (arr_pos = 0; arr_pos < num_barcode; arr_pos++) {
      base = input_barcodes[arr_pos]->sequence[index];
      // using the position of this base, insert this barcode into the sorted array, and increment that position
      sorted_barcodes[positions[Base_to_Int(&base)]++] = input_barcodes[arr_pos];
  }

  // now, sorted_barcodes should now contain all values from input_barcodes but sorted
  // make input_barcodes reflect the values of sorted_barcodes
  int j;
  for (j= 0; j < num_barcode; j++) {
    input_barcodes[j] = sorted_barcodes[j];
  }
}

void 
Sort_Barcodes(int num_barcode, int barcode_length) {
  /* 
  Implements Radix Sort, which performs count sort on an array
  of barcodes, on each subsequent base from right to left.
  At end, barcodes array contains the same structs but sorted lexographically
  Makes use of an intermediate sort temporary array, which holds references to existing structs.
  These should not be freed on function termination, as doing so will destroy our barcode data.
    We only need to free the intermediate array pointers.
  */
  // Create our storage for the intermediate steps of count sort.
  a_barcode** temporary_barcodes = new a_barcode *[num_barcode];

  long i;
  // run radix_sort on this local array of the barcodes, freeing memory as we go
  for (i = barcode_length; i >= 0; i--) {
    // run count sort, which will sort based on the given index
    Count_Sort_Barcodes(i, barcodes, temporary_barcodes, num_barcode, barcode_length);
    // the barcodes array is now stable sorted based on the ith index
  }
  // free the array we created.
  delete [] temporary_barcodes;

}


void Search_Barcodes_At_Index(Trie* this_trie, string filename, int index, int barcode_length, 
                                int num_reads_search, long *barcodes_found, long *barcodes_not_found)
{
    // this function will only check for barcodes at the given index
    // search each line of the file (up to a certain point) and record matches at the required position.
    /*
    fstream file;
    file.open(filename, ios::in);

    string line;
    */
   
    gzFile fq_gz;
    char c_line[MAX_LL];

    fq_gz = gzopen(filename.c_str(), "r");

    int barcode_index = -1;
    long line_count = 0, found = 0, not_found = 0;
    string line;
    //while (getline(file, line) && (line_count / 4 < num_reads_search))
    gzgets(fq_gz, c_line, MAX_LL);
    while (!gzeof(fq_gz) && (line_count / 4 < num_reads_search))
    {
        line_count++;

        if (line_count % 4 == 2)
        {
            // line_count % 4 == 2 is the read line we're after
            line = string(c_line);
            // given that we are in the actual read of the fastq file (the second line of every group of 4), check for a barcode at the given position
            barcode_index = 
                this_trie->Locate_Seq_At_Pos(line, index, barcode_length);
            if (barcode_index != -1) {
                found++;
            } else {
                not_found++;
            }
        }
        // read the next line before looping
        gzgets(fq_gz, c_line, MAX_LL);
    }

    gzclose(fq_gz);

    *barcodes_found = found;
    *barcodes_not_found = not_found;
}

ResizeArray *Search_Barcodes_Section_Read(Trie* this_trie, string filename, int index, int barcode_length, long num_to_check,
                                long *out_found, long *out_not_found) {
    // should this function search through the entire read or only a small subset like the above function
    /*
    fstream file;
    file.open(filename, ios::in);

    string line;
    */
    gzFile fq_gz;
    char c_line[MAX_LL];

    fq_gz = gzopen(filename.c_str(), "r");
    int barcode_index = -1, found_position = -1;
    long line_count = 0, found = 0, not_found = 0;
    string line;

    //barcodes_found b_found = {0, 0, 0, ResizeArray(100)};
    ResizeArray *positions = new ResizeArray(100);
    // search each line and record the position barcode was found in.
    gzgets(fq_gz, c_line, MAX_LL);
    while (!gzeof(fq_gz) && (line_count / 4 < num_to_check)) {
        line_count++;

        if (line_count % 4 == 2) {
            line = string(c_line);
            // scan through the given section of the read and check for barcode matches
            barcode_index = this_trie->Locate_Seq_Subsection(line, 0, index+10, &found_position);

            if (barcode_index != -1) {
                positions->Increment(found_position);
                found++;
            } else not_found++;
        }

        gzgets(fq_gz, c_line, MAX_LL);
    }

    gzclose(fq_gz);

    *out_found = found;
    *out_not_found = not_found; 
    return positions;
}


void Clean_Up(int num_barcode)
{
    /*
  Deallocate all space for arrays created
  */
    // free the barcode array
    for (int i = 0; i < num_barcode; i++) {
        delete barcodes[i];
    }
    delete [] barcodes;
}

// [[Rcpp::export]]
bool check_barcode_reads(Rcpp::String fastq, Rcpp::String barcodeseqs, Rcpp::String barcodeRealname,
                int barcode_start, int barcode_length,
                int lines_to_search, double threshold) {

    int num_barcode;
    bool finish_program = false;
    Trie* barcodes_trie = new Trie;
    try {
        string barcode_file = barcodeseqs.get_cstring();
        // First we need to read in all the barcodes
        num_barcode = Read_In_Barcodes(barcode_file);

        // then we need to sort the barcodes
        // Sort_Barcodes(num_barcode, barcode_length); // only need to sort if we are allowing mismatches in the check
        //barcode_single_trie_head = Build_Trie_Barcodes(num_barcode, barcode_length); // are the barcodes going to be paired or not paired?? should I remove the dual indexed option?
        Build_Trie_Barcodes(barcodes_trie, num_barcode);
        //long lines_to_search = 100000; 
        string fastq_file = fastq.get_cstring();

        long found, not_found;
        double search_result, new_search_result;
        // then we need to actually search the file lines.
        // quick search of every read line (to a point)
        Search_Barcodes_At_Index(barcodes_trie, fastq_file, barcode_start, barcode_length, lines_to_search, &found, &not_found);

        // once we've check the given position for barcodes, if we have a certain hit amount, proceed with program
        // if we don't reach the threshold, we should search the whole read for barcode locations and give a break down
        // of the best spots
        search_result = (double) found / (double) (found + not_found);
        if (search_result >= threshold) {
            Rcpp::Rcout << "Successful; continuing with program.\n";
            finish_program = true;
        } else {
            // we now need to search through a subset of the read to locate a better position for the barcode start
            ResizeArray *positions_found = 
                Search_Barcodes_Section_Read(barcodes_trie, fastq_file, barcode_start, barcode_length, lines_to_search, &found, &not_found);
            long max_value;
            int max_position = positions_found->Max(&max_value);
            new_search_result = (double) found / (double) (found + not_found);

            if (new_search_result >= .5) {
                Rcpp::Rcout << "Invalid barcode start index given, with only " <<  search_result * 100
                            << " percent of reads containing a barcode match. However, a better barcode start location is " << max_position
                            << " , where " << new_search_result * 100 
                            << " percent of barcodes were found.\n";
            } else {
                Rcpp::Rcout << "Unsuccessful. No location was found with a high number of barcode matches. Did both " << barcodeRealname.get_cstring()
                            << " and " << fastq_file.c_str() 
                            << " come from the same provider?\n";
            }
            positions_found->Delete();
            delete positions_found;
        }

        Clean_Up(num_barcode);
        barcodes_trie->Clear_Trie();
        delete barcodes_trie;
    } catch (std::exception &e) {
        Rcpp::stop(e.what());
    }
    return finish_program;
}
