#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

#define BLOCKSIZE 10000000
#define TERMINATOR '@'
// below here to // [[Rccp::export]] is C code, not C++. As it is converted and tested, move this down
/*typedef struct
{
   char *sequence;
    char *sequence2;   // dual indexed reads
    char *sequenceRev; // paired reads
    int original_pos;
} a_barcode;
typedef struct
{
    int original_seq_index;
    int current_seq_index;
} end_node;
typedef struct trie_node trie_node;
struct trie_node
{
    char base;
    long count;
    // links is [@, A, C, G, T]
    trie_node *links[5];
    end_node *end;
};
*/

struct a_barcode {
    string sequence;
    int original_pos;
};

struct end_node {
    int original_seq_index;
    int current_seq_index;
};

struct trie_node {
    char base;
    long count;
    trie_node *links[5]; // array of pointers to trie nodes
    end_node *end;
};

a_barcode **barcodes;
trie_node *barcode_single_trie_head;
trie_node *barcode_paired_trie_head;



////////////// Functions for reading barcodes and hairpins from input files, and storing in struct arrays
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






void Read_In_Barcodes(string filename, int num_barcode)
{
    /*
  Read in the barcodes from a given textfile, which needs to contain newline char seperated barcode entries.
  All barcodes must be the same length.
  Barcodes get read into an a_barcode struct, along with other information,
  and the struct is then stored in the barcodes array.
  filename: a string of the file name to read from
  return: void.
  post-condition: barcodes contains structs of all the barcodes in the given file.
  */
    Rprintf("\tRead Function, file: %s\n", filename.c_str());
    fstream file;
    //size_t len = 1000;
    file.open(filename, ios::in);

    barcodes = new a_barcode *[num_barcode];
    string line; /// allocate space for the line
    a_barcode *new_barcode;                          /// create a new_barcode struct variable
    
    int count = 0;
    //char *token;
    while (getline(file, line))
    {
        //new_barcode = (a_barcode *)malloc(sizeof(a_barcode));
        new_barcode = new a_barcode; // new makes a pointer
        new_barcode->sequence = line;
        new_barcode->original_pos = count;
        /*
        if (is_PairedReads > 0)
        {
            // strtok returns the first token from the string up to the first sep \t. Subsequent calls return the next token
            token = strtok(line, "\t");
            token = strtok(NULL, "\t");
            new_barcode->sequenceRev = (char *)malloc(barcode_length_rev * sizeof(char));
            strncpy(new_barcode->sequenceRev, token, barcode_length_rev);
        }
        else if (is_DualIndexingReads > 0)
        {
            token = strtok(line, "\t");
            token = strtok(NULL, "\t");
            new_barcode->sequence2 = (char *)malloc(barcode_length_rev * sizeof(char));
            strncpy(new_barcode->sequence2, token, barcode2_length);
        }
        else
        {
        new_barcode->sequenceRev = NULL;
        new_barcode->sequence2 = NULL;
        };
        */
        
        barcodes[count] = new_barcode;
        count++;
    }
    file.close();
    //delete(line);

    for (int i = 0; i < 2; i++) {
        Rprintf("%s\n", barcodes[i]->sequence.c_str());
    }
    Rprintf(" -- Number of Barcodes : %d\n", num_barcode);
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
          break;
        case 'C':
          return 2;
          break;
        case 'G':
          return 3;
          break;
        case 'T':
          return 4;
          break;
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


void Print_Barcodes(int num_barcode) {
    Rprintf("Barcode Array:\n");
    for (int i = 0; i < num_barcode; i++) {
        Rprintf("\tBarcode # %d: %s\n", barcodes[i]->original_pos, barcodes[i]->sequence.c_str());
    }
    Rprintf("End Barcode Array.\n");
}

// [[Rcpp::export]]
void check_barcode_reads(String fastq, String barcodeseqs, 
                int barcodeStart, int barcodeLength) {

    int num_barcode;
    int barcode_length;
    long num_read;
    Rprintf("Fastq given: %s\n", fastq.get_cstring());
    Rprintf("Barcodes given: %s\n", barcodeseqs.get_cstring());

    string barcode_file = barcodeseqs.get_cstring();
    num_barcode = Get_Lines_In_File(barcode_file);
    // First we need to read in all the barcodes
    Read_In_Barcodes(barcode_file, num_barcode);

    Print_Barcodes(num_barcode);
    // then we need to sort the barcodes
    // IF BARCODE SORTING NEEDS REVERSE READS AND STUFF THEN WE SHOULD USE OLD SORTING METHOD. NOT COUNTING SORT
    // Sort_Barcodes(num_barcode, barcode_length); // THIS NEEDS TO BE PROPERLY FIGURED OUT. SHOULD BE ABLE TO USE COUNT SORT FROM HAIRPINR SORTING.
    // and the build the barcode trie for searching (each barcode node contains the current position in the barcode array, as well as the original position in the barcode array,
    // which is why we need to sort them before building the trie, otherwise the current and original will be the same but the barcodes will be in a different position in the array)
    barcode_single_trie_head = Build_Trie_Barcodes(false, false); // are the barcodes going to be paired or not paired?? should I remove the dual indexed option?

    int lines_to_search = 500;
    int barcodes_found = 0;
    int barcodes_not_found = 0;
    // then we need to actually search the file lines.
    // quick search of every read line (to a point)
    //Search_Barcodes_At_Index(*file, *barcodeStart, &barcodes_found, &barcodes_not_found);

    Rprintf("Number of Barcodes Found: %d\n", barcodes_found);
    Rprintf("Number of reads with barcodes not found: %d\n", barcodes_not_found);
}

//NumericVector timesTwo(NumericVector x) {
//  return x * 2;
//}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
